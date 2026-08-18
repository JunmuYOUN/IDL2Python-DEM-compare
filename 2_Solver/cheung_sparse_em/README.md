# cheung_sparse_em — 적은 수의 성분으로 DEM 만들기

이 코드는 여러 후보 DEM 모양 중에서 **관측값을 설명하는 데 필요한 성분만 골라** 합친다.
계수는 음수가 될 수 없으며, 가능한 한 적은 수의 성분을 쓰도록 LP로 계산한다.

## 먼저 읽기

`aia_sparse_em_init`은 좁은 봉우리와 폭이 다른 Gaussian basis를 모아 후보 목록을 만든다.
`aia_sparse_em_solve`는 각 픽셀에서 관측 오차 범위를 만족하면서 계수 합이 가장 작은 조합을
찾는다. 여기서 “희소하다”는 말은 **0이 아닌 계수가 몇 개뿐**이라는 뜻이다.

결과 DEM은 매끈한 곡선보다 몇 개의 뾰족한 성분으로 보일 수 있다. 이것은 계산 오류가 아니라
방법 자체의 특징이다. `tolfac`을 키우면 관측값을 허용하는 범위가 넓어지고, 보통 더 적은 성분을
선택한다. 상태값이 1 또는 10인 픽셀은 계수가 모두 0이므로 결과에서 제외해야 한다.

IDL의 `simplex` 대신 Python에서는 SciPy의 `linprog`를 쓴다. 답이 여러 개인 경우 계수 배열은
달라도 다시 계산한 채널 밝기는 같을 수 있다. 변환본은 두 시험 세트에서 각각 15개 비교 항목과
회귀 테스트 4개를 통과했다.

- IDL 원본: `idl/`
- Python 변환본: `python/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 방법 요약 — LP 기반 sparse inversion
2. 왜 L1인가 — sparsity의 의미
3. `aia_sparse_em_init` — 사전(후보 모양 모음) 만들기
4. `aia_sparse_em_solve` — 픽셀별 LP
5. 변환의 핵심 — IDL simplex를 scipy로 바꾸다
6. 변환 및 검증 범위

## 1. 방법 요약 — LP 기반 sparse inversion

Cheung 방법은 least-squares 대신 LP로 DEM inversion을 풉니다.

**DEM 복원을 선형계획법(Linear Programming)으로 푼다.**

**최소화:** 기저 계수의 합 (= L1 norm, 계수가 모두 ≥0이므로)

**제약:** 재구성한 DN이 관측 DN의 허용 밴드 `[y−tol, y+tol]` 안에 있을 것

**제약:** 모든 계수 ≥ 0

이 형태를 **basis pursuit**이라 부르며, compressed sensing 계열의 표준 문제입니다.

|  | 최소제곱 계열 | sparse-EM |
|---|---|---|
| 목적함수 | Σ(잔차)² 최소화 | Σ(계수) 최소화 — 잔차가 아님 |
| 데이터 적합 | 목적함수 안에 포함 (χ²) | 부등식 제약으로 분리 |
| 양수성 | log 치환 / 잘라내기 / 없음 | LP의 변수 하한 ( x ≥ 0 ) |
| 해의 성격 | 매끄러운 곡선 | 희소(sparse) — 대부분의 계수가 정확히 0 |
| 풀이 도구 | Cholesky / SVD / 반복 | simplex (LP 솔버) |

## 2. 왜 L1인가 — sparsity의 의미

계수의 합을 최소화하는 L1 regularization은 다음과 같이 작동합니다.

제약조건(데이터 밴드)을 만족하는 해는 무수히 많습니다. 그중에서 **계수 합이 가장 작은 것**을 고르면, 수학적으로 **대부분의 계수가 정확히 0인 해**가 선택됩니다. 이것이 L1 최소화의 잘 알려진 성질입니다.

**물리적 해석:** "이 픽셀의 플라즈마를 설명하는 데 필요한 **최소한의 온도 성분**만 쓴다."

Tikhonov 계열이 **"매끄러운 DEM"**을 선호한다면, sparse-EM은 **"성분 개수가 적은 DEM"**을 선호합니다. 코로나 플라즈마가 몇 개의 구별되는 온도 구조로 이뤄져 있다는 관점입니다.

실제 검증에서 이 성질이 확인됩니다 — 기저가 **3×nt개**나 되는데 실제로 0이 아닌 계수는 **set1에서 3개, set2에서 10개**뿐이었습니다.

## 3. aia_sparse_em_init — 사전(후보 모양 모음) 만들기

```
def aia_sparse_em_init(tresp, lgtaxis, bases_sigmas=(0.0, 0.1, 0.2), dictfac=1e26):
    nchannels, ntemp = tresp.shape
    nsigmas = len(bases_sigmas)
    basis_funcs = np.zeros((nsigmas * ntemp, ntemp))
    for s, sig in enumerate(bases_sigmas):
        if sig == 0.0:
            for m in range(ntemp):
                basis_funcs[ntemp*s + m, m] = 1.0                      # ① 델타 기저
        else:
            for m in range(ntemp):
                line = np.exp(-((lgtaxis - lgtaxis[m]) / sig) ** 2.0)  # ② 가우시안 기저
                line[line < 0.04] = 0.0                                # ③ 꼬리 절단
                basis_funcs[ntemp*s + m, :] = line
    Dict = (tresp @ basis_funcs.T) * dictfac                           # ④ 사전
    return Dict, basis_funcs
```

### ① 델타 기저 (σ = 0)

온도격자 한 점에만 1이고 나머지는 0인 기저입니다. 결과적으로 **단위행렬**이 되죠. "폭이 0인 아주 좁은 온도 성분"을 표현합니다.

### ② Gaussian 기저 (σ > 0)

```
line = exp( −((logT − logT_m) / σ)² )
```

> **⚠ 지수의 형태에 주의.** 통상의 정규분포는 `exp(−x²/2σ²)` 인데, 여기는 `exp(−(x/σ)²)` 입니다. 즉 이 `σ`는 표준편차가 아니라 **폭 파라미터**이고, 표준편차로 환산하면 `σ/√2` 입니다. 값을 물리적으로 해석할 때 이 차이를 감안해야 합니다.

### ③ 꼬리 절단

```
line[line < 0.04] = 0.0
```

진폭이 4% 미만인 꼬리를 잘라 **기저를 국소화**합니다. 기저가 온도축 전체에 퍼져 있으면 LP가 무의미하게 커지고, sparsity의 의미도 흐려집니다.

### ④ 사전 Dict

```
Dict = (tresp @ basis_funcs.T) * dictfac        # (nchan, nbasis), dictfac = 1e26
```

`Dict[c, k]` = "기저 `k`가 계수 1일 때 채널 `c`가 받는 DN". 순방향 모델을 미리 계산해 표로 만든 것으로, compressed sensing 용어로 **후보 모양 모음**입니다.

`dictfac = 1e26`은 **단위 환산 상수**입니다. EM 값이 10²⁶~10²⁸ 규모라 그대로 두면 LP 솔버의 수치 조건이 나빠지므로, 계수가 O(1) 근처에서 놀도록 스케일을 맞춥니다.

기저 개수는 `nsigmas × ntemp` = **3 × nt**입니다. 같은 온도에 대해 폭이 다른 3가지 버전을 모두 후보로 두고, LP가 알아서 고르게 합니다.

## 4. aia_sparse_em_solve — 픽셀별 LP

### LP 구성 50–57

```
c    = np.ones(nbasis)                   # 목적: minimize Σx  (= L1)
A_ub = np.vstack([Dict, -Dict])          # 부등식 제약 행렬 (루프 밖, 1회)
for 픽셀:
    y    = np.maximum(image[i,j,:], 0.0)                          # ① 음수 관측 → 0
    tol  = tolfac * np.sqrt(y)                                    # ② 허용 밴드 폭
    b_ub = np.concatenate([y + tol,
                           -np.maximum(y - symmbuff*tol, 0.0)])   # ③ 상·하한
    res  = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=(0, None), method="highs")
```

**② 허용 밴드 — χ² 대신 쓰는 것**

```
tol = tolfac × √y          (tolfac = 1.4)
```

`√y`는 Poisson 잡음의 크기입니다. 즉 **"관측값에서 잡음의 1.4배 이내면 맞는 것으로 친다"**는 뜻이고, 다른 solver의 χ²≈1 조건과 **역할이 같습니다.** 다만 여기서는 목적함수가 아니라 **제약**으로 들어갑니다.

**③ 부등식 두 줄로 밴드 만들기**

```
Dict·x ≤ y + tol                    ← 상한
```

`Dict·x ≥ y − tol` ⟺ `−Dict·x ≤ −(y − tol)   ← 하한 (부호 뒤집어 상한 형태로)`

`linprog`는 `≤` 형태만 받으므로, 하한을 표현하려고 `−Dict`를 아래에 쌓았습니다. 그래서 `A_ub`가 `[Dict;` `−Dict]` 입니다.

하한에는 `np.maximum(..., 0.0)`이 한 번 더 걸려 있습니다 — `y − tol`이 음수가 되면 **하한을 0으로 완화**합니다. 어두운 픽셀에서 제약이 모순되는 것을 막는 장치입니다.

### 결과 처리 58–69

```
if res.success:  x = res.x;  s = 0
else:            x = zeros;  s = 1                    # LP 실패
if x.min() < 0.0:                                     # 음수 해 → 전체를 0으로
    x = np.zeros(nbasis);  s = 10
coeffs[i,j,:] = x
zmax[i,j]     = -np.sum(x)                            # IDL r[0] 규약
status[i,j]   = s
```

> **⚠** `status = 10` **분기.** LP의 변수 하한이 이미 `x ≥ 0` 이라 이론상 음수 해는 나올 수 없습니다. 하지만 IDL `simplex`는 수치오차로 아주 작은 음수를 낼 수 있고, 원본은 그럴 때 **해 전체를 버리고 status 10을** **기록**합니다. 변환본은 **실제로 발생하지 않더라도 이 규약을 그대로 재현**했습니다 — IDL과 같은 입력에서 같은 status 값이 나와야 결과 일치가 성립하기 때문입니다.

`zmax = −Σx`의 음수 부호도 IDL 관용구입니다. IDL `simplex`는 **최대화** 문제로 정의되므로, 최소화하려면 목적함수의 부호를 뒤집어 넣고 결과도 음수로 나옵니다.

### EM 곡선 복원 73–76

```
oem[i,j,:] = coeffs[i,j,:] @ basis_funcs
```

계수에 기저함수를 곱해 **온도축 위의 실제 EM 곡선**을 만듭니다. `coeffs`는 (3·nt)개, `oem`은 nt개입니다.

## 5. 변환의 핵심 — IDL simplex를 scipy로 바꾸다

이 패키지에서 가장 위험했던 부분입니다.

```
IDL     내장 simplex  (Numerical Recipes의 simplx 구현)
        ↓
Python  scipy.optimize.linprog(method='highs')
```

수학적으로 **같은 문제**를 풀지만, **알고리즘 구현이 완전히 다릅니다.**

### 왜 위험한가 — 퇴화(degeneracy)

sparse L1 LP에서는 목적값(Σx)이 같지만 계수 조합이 서로 다른 **degenerate optimum**이 생길 수 있습니다.

simplex 계열 알고리즘이 어느 vertex를 반환하는지는 구현에 따라 달라질 수 있습니다. 따라서
IDL과 SciPy의 coefficient가 달라도 objective와 reconstructed DN이 같으면 동등한 해일 수 있습니다.

### 검증 결과

| 세트 | 성격 | 결과 |
|---|---|---|
| set1 | 단봉 Gaussian, nonzero 3개 | 동일 vertex, coeffs max_rel ≈ 3e-7 |
| set2 | 광폭 이중봉 — 퇴화가 잘 일어나는 케이 스, nonzero 10개 | 동일 vertex, coeffs max_rel ≈ 3e-7 |

두 케이스 모두 scipy HiGHS와 IDL simplex가 **같은 sparse vertex**를 float32 정밀도로 반환했습니다. 특히 set2는 퇴화를 유도하려고 **일부러 고른 케이스**인데도 일치했습니다.

**일치하지 않았다면 어떻게 했을까.** 비교 보고서에 원칙이 명시돼 있습니다 — 계수가 달라도 **목적값(** `zmax` **)**과 **재구성 DN(** `Dict·coeff` **)**이 일치하면 "다르지만 동등한 해"로 판정하고 **허용한 차이로 기록**합니다.

현재는 미발생이라 **허용한 차이 0건**입니다. 다만 다른 데이터에서는 발생할 수 있으므로, **계수 하나하나가 항상 같으리라 기대해서는 안 됩니다.**

## 6. 변환 및 검증 범위

### 결과

| 세트 | 규모 | 통과/실패 | 판정 |
|---|---|---|---|
| set1 | 1픽셀 (단봉) | 15 / 0 | 모두 통과 |
| set2 | 2픽셀 (단봉 + 광폭 이중봉) | 15 / 0 | 모두 통과 |

엄격한 tolerance (float 1e-5) 적용

`basis_funcs` · `Dict` 구성: max_rel ~**1e-8** (순수 산술이라 매우 정확)

LP `coeffs` / `oem`: max_rel ~**3e-7**

회귀 테스트 **4/4** (IDL 기준값 고정 set1·set2 + 기저함수 delta/gaussian 가드)

**허용한 차이 0건**

### 외부에서 입력받은 값

`tresp`, `lgtaxis` — SSW `aia_get_response` 유래 ( `01_init` )

`image` — Gaussian EM forward model 산출 ( `00_input` )

나머지( `basis_funcs`, `Dict`, LP 전체)는 **포팅해서 검증**했습니다.

### 변환하거나 검증하지 않은 부분

| 항목 | 사유 |
|---|---|
| XRT/EIS/display/moments/plot 루틴 | DEM 계산과 직접 관련이 없어 제외. 사용한 SSW 설치본에서는 `make_xrt_*`도 구문 오류로 컴파일할 수 없었음 |
| `adaptive_tolfac` 재귀 | `adaptive_tolfac=0`인 단일 simplex 실행만 검증. 실패한 픽셀의 `tolfac`을 1.5배로 높여 다시 푸는 경로는 별도로 검증하지 않음 |
| bases_sigmas  기본값 | 튜토리얼 값 [0, 0.1, 0.2] 로 검증. 기본값 [0, 0.1, 0.2, 0.6] 은 미검 증이나 구조가 동일 |
| 맵 단위 처리 | 단일/2픽셀 규모. 실관측 맵은 후속 |

### 사용 시 주의사항

**해는 sparse합니다.** 매끄러운 DEM 곡선 대신 몇 개의 좁은 성분이 선택될 수 있습니다.

**Degenerate optimum에서는 coefficient가 유일하지 않습니다.** 구현을 비교할 때는 objective와 reconstructed DN도 함께 확인해야 합니다.

`tolfac`은 허용 오차 범위를 조절합니다. 값을 키우면 band가 넓어져 보통 더 sparse한 해가 나오고, 줄이면 관측값을 더 엄격히 따릅니다.

Basis dictionary에 없는 온도 구조는 표현할 수 없습니다. 검증에서는 `bases_sigmas`가 정한 세 종류의 폭만 사용했습니다.

LP 실패(`status = 1`) 또는 음수 해(`status = 10`)로 표시된 픽셀은 coefficient가 모두 0이므로 mask에서 제외해야 합니다.

**기록의 근거.** 코드 설명은 `idl/aia_sparse_em_init.pro`와 `python/aia_sparse_em.py`를
바탕으로 작성했다. 검증 수치는 변환 작업 당시의 비교 보고서를 따른다.

</details>

## 참고 논문

- Cheung et al. (2015), [Thermal Diagnostics with the Atmospheric Imaging Assembly onboard the Solar Dynamics Observatory: A Validated Method for Differential Emission Measure Inversions](https://doi.org/10.1088/0004-637X/807/2/143)
