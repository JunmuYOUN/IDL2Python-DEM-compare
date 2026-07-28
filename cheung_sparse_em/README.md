# cheung sparse_em 코드 설명

Cheung et al. 2015 — Sparse-EM / basis-pursuit DEM

SpaceAI-DEM / `work/cheung_sparse_em_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — 최소제곱이 아니라 선형계획법
2. 왜 L1인가 — sparsity의 의미
3. `aia_sparse_em_init`  — 사전(dictionary) 만들기
4. `aia_sparse_em_solve`  — 픽셀별 LP
5. 변환의 핵심 — IDL simplex를 scipy로 바꾸다
6. 검증 결과와 미변환 경계

## 1. 한눈에 — 최소제곱이 아니라 선형계획법

지금까지의 solver는 모두 **제곱오차를 줄이는** 문제였습니다. Cheung의 방법은 문제의 종류 자체가 다릅니다.

**DEM 복원을 선형계획법(Linear Programming)으로 푼다.**

**최소화:** 기저 계수의 합 (= L1 norm, 계수가 모두 ≥0이므로)

**제약:** 재구성한 DN이 관측 DN의 허용 밴드 `[y−tol, y+tol]`  안에 있을 것

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

"계수의 합을 최소화한다"는 게 왜 의미 있는 정규화일까요?

제약조건(데이터 밴드)을 만족하는 해는 무수히 많습니다. 그중에서 **계수 합이 가장 작은 것**을 고르면, 수학적으로 **대부분의 계수가 정확히 0인 해**가 선택됩니다. 이것이 L1 최소화의 잘 알려진 성질입니다.

**물리적 해석:** "이 픽셀의 플라즈마를 설명하는 데 필요한 **최소한의 온도 성분**만 쓴다."

Tikhonov 계열이 **"매끄러운 DEM"**을 선호한다면, sparse-EM은 **"성분 개수가 적은 DEM"**을 선호합니다. 코로나 플라즈마가 몇 개의 구별되는 온도 구조로 이뤄져 있다는 관점입니다.

실제 검증에서 이 성질이 확인됩니다 — 기저가 **3×nt개**나 되는데 실제로 0이 아닌 계수는 **set1에서 3개, set2에서 10개**뿐이었습니다.

## 3. aia_sparse_em_init — 사전(dictionary) 만들기

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

### ② 가우시안 기저 (σ > 0)

```
line = exp( −((logT − logT_m) / σ)² )
```

> **⚠ 지수의 형태에 주의.** 통상의 정규분포는 `exp(−x²/2σ²)` 인데, 여기는 `exp(−(x/σ)²)` 입니다. 즉 이 `σ` 는 표준편차가 아니라 **폭 파라미터**이고, 표준편차로 환산하면 `σ/√2` 입니다. 값을 물리적으로 해석할 때 이 차이를 감안해야 합니다.

### ③ 꼬리 절단

```
line[line < 0.04] = 0.0
```

진폭이 4% 미만인 꼬리를 잘라 **기저를 국소화**합니다. 기저가 온도축 전체에 퍼져 있으면 LP가 무의미하게 커지고, sparsity의 의미도 흐려집니다.

### ④ 사전 Dict

```
Dict = (tresp @ basis_funcs.T) * dictfac        # (nchan, nbasis), dictfac = 1e26
```

`Dict[c, k]`  = "기저 `k` 가 계수 1일 때 채널 `c` 가 받는 DN". 순방향 모델을 미리 계산해 표로 만든 것으로, compressed sensing 용어로 **dictionary**입니다.

`dictfac = 1e26` 은 **단위 환산 상수**입니다. EM 값이 10²⁶~10²⁸ 규모라 그대로 두면 LP 솔버의 수치 조건이 나빠지므로, 계수가 O(1) 근처에서 놀도록 스케일을 맞춥니다.

기저 개수는 `nsigmas × ntemp`  = **3 × nt**입니다. 같은 온도에 대해 폭이 다른 3가지 버전을 모두 후보로 두고, LP가 알아서 고르게 합니다.

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

`√y` 는 Poisson 잡음의 크기입니다. 즉 **"관측값에서 잡음의 1.4배 이내면 맞는 것으로 친다"**는 뜻이고, 다른 solver의 χ²≈1 조건과 **역할이 같습니다.** 다만 여기서는 목적함수가 아니라 **제약**으로 들어갑니다.

**③ 부등식 두 줄로 밴드 만들기**

```
Dict·x ≤ y + tol                    ← 상한
```

`Dict·x ≥ y − tol` ⟺ `−Dict·x ≤ −(y − tol)   ← 하한 (부호 뒤집어 상한 형태로)`

`linprog` 는 `≤` 형태만 받으므로, 하한을 표현하려고 `−Dict` 를 아래에 쌓았습니다. 그래서 `A_ub` 가 `[Dict;` `−Dict]` 입니다.

하한에는 `np.maximum(..., 0.0)` 이 한 번 더 걸려 있습니다 — `y − tol` 이 음수가 되면 **하한을 0으로 완화**합니다. 어두운 픽셀에서 제약이 모순되는 것을 막는 장치입니다.

### 결과 처리 58–69

```
if res.success:  x = res.x;  s = 0
else:            x = zeros;  s = 1                    # LP 실패
if x.min() < 0.0:                                     # 음수 해 → 통째로 0
    x = np.zeros(nbasis);  s = 10
coeffs[i,j,:] = x
zmax[i,j]     = -np.sum(x)                            # IDL r[0] 규약
status[i,j]   = s
```

> **⚠** `status = 10` **분기.** LP의 변수 하한이 이미 `x ≥ 0` 이라 이론상 음수 해는 나올 수 없습니다. 하지만 IDL `simplex` 는 수치오차로 아주 작은 음수를 낼 수 있고, 원본은 그럴 때 **해 전체를 버리고 status 10을** **기록**합니다. 변환본은 **실제로 발생하지 않더라도 이 규약을 그대로 재현**했습니다 — IDL과 같은 입력에서 같은 status 값이 나와야 parity가 성립하기 때문입니다.

`zmax = −Σx` 의 음수 부호도 IDL 관용구입니다. IDL `simplex` 는 **최대화** 문제로 정의되므로, 최소화하려면 목적함수의 부호를 뒤집어 넣고 결과도 음수로 나옵니다.

### EM 곡선 복원 73–76

```
oem[i,j,:] = coeffs[i,j,:] @ basis_funcs
```

계수에 기저함수를 곱해 **온도축 위의 실제 EM 곡선**을 만듭니다. `coeffs` 는 (3·nt)개, `oem` 은 nt개입니다.

## 5. 변환의 핵심 — IDL simplex를 scipy로 바꾸다

이 패키지에서 가장 위험했던 부분입니다.

```
IDL     내장 simplex  (Numerical Recipes의 simplx 구현)
        ↓
Python  scipy.optimize.linprog(method='highs')
```

수학적으로 **같은 문제**를 풀지만, **알고리즘 구현이 완전히 다릅니다.**

### 왜 위험한가 — 퇴화(degeneracy)

sparse L1 LP는 **최적점이 여러 개인 경우가 흔합니다.** 목적값(Σx)은 똑같은데 그 값을 달성하는 계수 조합이 여러 가지인 상황이죠. 이를 **degenerate optimum**이라 합니다.

simplex 계열 알고리즘은 그중 **어느 꼭짓점(vertex)에 착지하느냐**가 구현에 따라 달라집니다. IDL과 scipy가 **서로 다르지만 똑같이 옳은 답**을 낼 수 있다는 뜻입니다.

### 실제 결과 — 다행히 일치

| 세트 | 성격 | 결과 |
|---|---|---|
| set1 | 단봉 가우시안, nonzero 3개 | 동일 vertex, coeffs max_rel ≈ 3e-7 |
| set2 | 광폭 이중봉 — 퇴화가 잘 일어나는 케이 스, nonzero 10개 | 동일 vertex, coeffs max_rel ≈ 3e-7 |

두 케이스 모두 scipy HiGHS와 IDL simplex가 **같은 sparse vertex**를 float32 정밀도로 반환했습니다. 특히 set2는 퇴화를 유도하려고 **일부러 고른 케이스**인데도 일치했습니다.

**일치하지 않았다면 어떻게 했을까.** 인증서에 원칙이 명시돼 있습니다 — 계수가 달라도 **목적값(** `zmax` **)**과 **재구성 DN(** `Dict·coeff` **)**이 일치하면 "다르지만 동등한 해"로 판정하고 **waiver로 기록**합니다.

현재는 미발생이라 **waiver 0건**입니다. 다만 다른 데이터에서는 발생할 수 있으므로, **계수 하나하나가 항상 같으리라 기대해서는 안 됩니다.**

## 6. 검증 결과와 미변환 경계

### 결과

| 세트 | 규모 | pass/fail | 판정 |
|---|---|---|---|
| set1 | 1픽셀 (단봉) | 15 / 0 | ALL PASS |
| set2 | 2픽셀 (단봉 + 광폭 이중봉) | 15 / 0 | ALL PASS |

strict tolerance (float 1e-5) 적용

`basis_funcs` · `Dict`  구성: max_rel ~**1e-8** (순수 산술이라 매우 정확)

LP `coeffs` / `oem` : max_rel ~**3e-7**

회귀 테스트 **4/4** (오라클 고정 set1·set2 + 기저함수 delta/gaussian 가드)

**waiver 0건**

### 주입한 것

`tresp` , `lgtaxis`  — SSW `aia_get_response`  유래 ( `01_init` )

`image`  — 가우시안 EM 순방향모델 산출 ( `00_input` )

나머지( `basis_funcs` , `Dict` , LP 전체)는 **포팅해서 검증**했습니다.

### 미변환 경계 (인증서 명시)

| 항목 | 사유 |
|---|---|
| XRT/EIS/display/moments/plot 루틴 | DEM solve와 무관. 더불어 이 SSW 설치본에서 make_xrt_* 가 구문오류로 컴 파일 불가 → 코어만 트림해 변환 |
| adaptive_tolfac  재귀 | 단일 simplex pass( adaptive_tolfac=0 )로 검증. 실패 픽셀을 tolfac  1.5배 로 재귀 재풀이하는 경로는 같은 solve의 반복이라 별도 parity 미실시 |
| bases_sigmas  기본값 | 튜토리얼 값 [0, 0.1, 0.2] 로 검증. 기본값 [0, 0.1, 0.2, 0.6] 은 미검 증이나 구조가 동일 |
| 맵 단위 처리 | 단일/2픽셀 규모. 실관측 맵은 후속 |

### 알고 써야 할 것

**해가 희소합니다.** 매끄러운 DEM 곡선을 기대하면 안 됩니다 — 몇 개의 뾰족한 성분이 나옵니다

**퇴화 최적점에서는 계수가 유일하지 않습니다.** 재현성이 필요하면 목적값과 재구성 DN으로 비교해야 합니다

`tolfac = 1.4` 가 사실상의 정규화 손잡이입니다. 키우면 밴드가 넓어져 더 희소해지고, 줄이면 데이터를 더 엄격히 따릅니다

기저 사전에 없는 온도 구조는 표현할 수 없습니다 ( `bases_sigmas` 가 정하는 폭 3종만 후보)

LP 실패( `status = 1` )와 음수해( `status = 10` ) 픽셀은 **계수가 전부 0**입니다. 반드시 마스킹하세요

**출처 표시.** 3~5장의 코드 설명은 `inbox/aia_sparse_em_init.pro` 와 `converted/aia_sparse_em.py` 에서, 6장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1·2장의 방법 분류와 L1 sparsity 설명은 `work/DEM_SOLVER_METHODS.md` 와 코드에 근거한 해석이며, Cheung et al. 2015 원논문 대조는 수행하지 않았습니다.
