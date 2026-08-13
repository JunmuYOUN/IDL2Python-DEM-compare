# chianti_dem — CHIANTI 선 세기로 DEM 맞추기

`dem_fit`은 관측한 스펙트럼 선 세기를 이용해 DEM을 맞춘다. 온도축에 몇 개의 기준점을 고정하고,
각 기준점의 DEM 값만 바꾸면서 관측값과 계산값의 차이를 줄인다. 기준점 사이의 DEM은 스플라인
곡선으로 잇고, 값 조정에는 LM을 쓴다.

## 먼저 읽기

- 입력: 관측 선 세기, 오차, 각 선의 온도 기여함수, DEM 기준점 위치와 초깃값
- 출력: 기준점의 DEM 값과 χ²
- 핵심 가정: 기준점 위치는 고정되며 DEM 모양은 그 점들을 잇는 스플라인으로 표현됨

이 방법은 초깃값과 기준점 위치에 영향을 받는다. 시작값이 좋지 않으면 다른 해에 멈출 수 있고,
관측 선의 수는 맞출 DEM 값의 수보다 많아야 한다. 이번 변환은 IDL의 장력 스플라인과 유한차분
계산까지 옮겼다.

합성 자료 1세트에서 최종 출력 15개는 허용 오차 안에서 맞았다. 첫 반복의 중간 행렬 2개는
아주 작은 차이를 허용했지만 최종 DEM과 χ²는 통과했다. 실제 CHIANTI 선 목록을 이용한 관측
자료는 이번에 확인하지 않았다.

- IDL 원본: `idl/`
- Python 변환본: `python/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 한눈에 — 노드 5개짜리 스플라인을 LM으로 맞춘다
2. `idl_spline` — IDL SPLINE 재현
3. 세 개의 부품 — `functn` / `deriv` / `chisqr`
4. Levenberg-Marquardt 루프
5. LM 반복의 검증 — `xrt_iterative`와의 차이
6. 검증 결과와 허용한 차이

## 1. 한눈에 — 노드 5개짜리 스플라인을 LM으로 맞춘다

CHIANTI의 고전적 DEM 적합 루틴입니다. 관측은 **스펙트럼 선의 강도**(이미저 채널이 아님)이고, DEM은 **몇 개의 노드를 지나는 스플라인 곡선**으로 표현합니다.

```
미지수:  log_dem_mesh   = 노드에서의 log DEM 값 (검증에서 5개)
고정:    log_t_mesh     = 노드의 온도 위치
DEM(T) = spline(log_t_mesh, log_dem_mesh) 를 세밀한 온도격자에 펼친 것
```

구조는 `xrt_iterative`와 마찬가지로 spline node를 최적화하는 방식입니다. 다만 두 패키지는
Python으로 옮긴 범위가 다르며, 그 차이는 5장에서 설명합니다.

| 단계 | 내용 |
|---|---|
| DEM 표현 | 노드 5개를 지나는 spline-under-tension |
| 순방향 | ∫ G(T)·DEM(T) dT  → 선 강도 예측 |
| 미분 | 유한차분(node를 하나씩 1% perturb) |
| 최적화 | Levenberg-Marquardt —.pro  파일 안에 직접 구현됨 |
| 종료 | χ² 개선폭 < 1e-5 또는 20회 |

## 2. `idl_spline` — IDL SPLINE 재현

이 패키지에서 가장 이례적인 부분입니다.

이 저장소의 다른 패키지들은 IDL 라이브러리 함수를 **외부 입력**으로 처리했습니다 ( `aia_get_response`, `int_tabulated`, `gaussian_function`, `eit_line_map` …).

그런데 `dem_fit`은 **매 반복마다 SPLINE을 호출**합니다. 외부 입력만으로는 해결이 안 되므로, **IDL SPLINE 자체를** **Python으로 완전히 재구현**했습니다 — 45줄짜리 `idl_spline` 함수가 그것입니다.

### 보통의 스플라인이 아닙니다

IDL의 `SPLINE`은 3차 스플라인이 아니라 **spline-under-tension**(장력 스플라인, 기본 `σ = 1.0` )입니다.

```
dels   = sigmap * delx1
exps   = np.exp(dels)
sinhs  = 0.5 * (exps - 1./exps)          ← sinh
diag1  = sinhin * (dels*0.5*(exps+1./exps) - sinhs)   ← cosh 포함
```

코드 전반에 `sinh` **·** `cosh`가 등장하는 이유입니다. 3차 다항식 대신 **쌍곡선 함수**를 씁니다. 장력 파라미터 `σ`가 커질수록 곡선이 직선에 가까워지고, 0에 가까우면 일반 3차 스플라인이 됩니다.

구조는 표준 스플라인과 같습니다:

1. **양 끝점 기울기 추정** (14–20줄) — 3점 유한차분 공식

2. **삼중대각 시스템 풀기** (33–39줄) — 전진 소거 후 후진 대입

3. **보간값 계산** (40–45줄) — `searchsorted`로 구간을 찾아 쌍곡선 보간

> **⚠ 전체가** `float32` **입니다.** 파일 맨 위의 `f32 = np.float32`가 코드 곳곳에 박혀 있는 이유입니다. IDL `SPLINE`이 float32로 계산하므로, 중간값 하나라도 float64로 새면 결과가 어긋납니다. 제어실험에서 **max_rel 1.6e-7**로 일치를 확인했는데, 이는 **float32의 ULP(최소 단위) 수준**입니다 — 즉 **가능한 한계까지** **일치**시킨 것이고, 이보다 더 좁힐 수는 없습니다.

## 3. 세 개의 부품

### functn — 순방향 모델 61–68

```
def functn(s):
    dem  = idl_spline(s.log_t_mesh, s.log_dem_mesh, s.log_dem_temp)   # ① 노드 → 곡선
    dem  = f32(10.0) ** dem                                            # ② log 해제
    dlnt = np.log(f32(10.0) ** s.d_dem_temp)                           # ③ 적분 폭
    for i in range(s.n_obs):
        out[i] = np.sum(s.ch_tot_contr[:, i] * dem * s.dem_temp) * dlnt  # ④ 적분
    return out
```

④의 적분은 xrt와 같은 변수변환입니다:

```
I = ∫ G(T)·DEM(T) dT ,   T = 10^t 이므로 dT = ln(10)·T·dt
  = Σ  G · DEM · T · (ln10 · dt)
              └ dem_temp ┘  └── dlnt ──┘
```

`dlnt = ln(10^d_dem_temp) = d_dem_temp · ln(10)` 이므로 정확히 일치합니다.

### deriv — 유한차분 Jacobian 70–80

```
da = np.maximum(f32(0.01) * s.log_dem_mesh, f32(1e-6))   # 노드값의 1%, 최소 1e-6
yfit = s.functn()
for j in range(nt):
    s.log_dem_mesh[j] = sav[j] + da[j]       # 노드 하나만 흔들고
    d[:, j] = (s.functn() - yfit) / da[j]    # 차이를 나눠 미분 근사
    s.log_dem_mesh[j] = sav[j]               # 원복
```

> **⚠ 해석적 미분이 아니라 유한차분입니다.** 노드 하나당 `functn()`을 한 번씩 더 부르므로, 반복마다 **스플라인 계산이 (노드 수 + 1)회** 일어납니다.

그리고 이것이 6장 허용한 차이의 원인입니다 — `idl_spline`의 **float32 ULP 수준 차이**가 **작은 수로 나누는** **유한차분**에서 증폭됩니다.

### chisqr  82–84

```
dy = s.obs_int - s.functn()
return np.float32(np.sum(dy ** 2 / s.obs_sig ** 2))
```

표준 χ²입니다. **reduced가 아닌 총합**이라는 점에 주의하세요.

## 4. Levenberg-Marquardt 루프 86–148

Bevington의 고전적 `CURFIT` 구조를 따릅니다.

### 곡률행렬과 기울기 쌓기 100–109

```
weight = 1.0 / np.maximum(s.obs_sig, f32(1.0))      # ← 93줄
for i in range(s.n_obs):
    for j in range(nterms):
        beeta[j] += weight[i] * dy[i] * deriv1[i, j]              # 기울기
        for k in range(j + 1):                                     # ← 하삼각만
            alpha[j, k] += weight[i] * deriv1[i, j] * deriv1[i, k] # 곡률 (JᵀJ)
alpha = alpha.T.copy()
```

`alpha`는 simple_reg_dem에서 본 `JᵀJ`와 같은 역할이고, `beeta`는 `Jᵀ·잔차`입니다.

> **⚠ 원본의 특징 두 가지** (변환본은 그대로 재현했습니다):

① **가중치가** `1/σ` **입니다** — `chisqr`가 쓰는 `1/σ²`와 다릅니다. 게다가 `max(σ, 1.0)` 바닥값이 걸려 있어, **σ < 1인 선은 가중치가 1로 고정**됩니다. LM은 하강 방향만 맞으면 수렴하므로 동작에는 문제가 없지만, 교과서적 Gauss-Newton은 아닙니다.

② `alpha`를 **하삼각만 채운 뒤 transpose**합니다. 대칭화가 아니라 전치이므로 결과는 **상삼각 행렬**이 되고 하삼각은 0으로 남습니다.

### flambda 재시도 루프 111–139

```
while True:
    # ① alpha를 대각 1로 정규화하고 flambda를 더함
    array[j,k] = alpha[j,k] / sqrt(alpha[j,j]*alpha[k,k])
    array[j,j] = 1.0 + flambda
    array = np.linalg.inv(array)
    # ② 스텝 적용
    log_dem_mesh[j] = b[j] + Σ_k beeta[k]*array[j,k]/sqrt(alpha[j,j]*alpha[k,k]) * scale
    chisqr = s.chisqr()
    if chisqr < chisq1:                    # ③ 개선됨 → 신뢰도 높이고 탈출
        flambda = flambda / 10.0;  break
    else:                                   # ④ 악화됨 → 더 보수적으로 재시도
        flambda = flambda * 10.0
        if scale > 0.1: scale = scale / 2.0
        continue
```

`flambda` **가 LM의 핵심입니다.** 대각선에 `1 + flambda`를 넣으면:

• `flambda`가 작으면 → Gauss-Newton에 가까움 (과감한 큰 단계)

• `flambda`가 크면 → 경사하강에 가까움 (안전한 작은 단계)

개선되면 과감하게( `/10` ), 실패하면 보수적으로( `×10` ) 전환하며 **스스로 보폭을 조절**합니다. simple_reg_dem의 trust region + line search와 목적이 같습니다.

`NaN` 처리(123–128줄)와 `flambda > 1e10` 탈출(135–136줄)이 무한루프 방지 장치입니다.

### 종료 조건 99, 143

```
while (it <= niter) and (dchisq >= dchisq_min) and not failed:
    ...
    dchisq = chisq1 - chisqr      # 이번 반복의 χ² 개선폭
```

χ² 개선폭이 `1e-5` 미만이 되거나 20회를 넘으면 종료합니다. `nfree = n_obs − nterms ≤ 0` 이면(관측보다 노드가 많으면) **시작도 안 하고 실패 처리**합니다(96–97줄).

## 5. LM 반복의 검증 — `xrt_iterative`와의 차이

`xrt_iterative`에서는 MPFIT을 이번 변환에서 제외했지만, `chianti_dem`에서는 LM 반복까지
옮겨 최종 결과를 검증했습니다. 두 경우의 차이는 다음과 같습니다.

|  | xrt_iterative | chianti_dem |
|---|---|---|
| LM의 소재 | MPFIT 라이브러리 호출 | dem_fit.pro  안에 직접 작성 |
| 변환 방법 | scipy로 대체해야 함 | 한 줄씩 그대로 포팅 |
| iteration path | MPFIT 구현에 따라 달라질 수 있음 | 같은 연산 순서를 재현할 수 있음 |
| 이번 변환 | LM 반복 제외 | LM 반복 포함 |

Path dependence 자체보다 optimizer 구현을 바꾸는 것이 결과 재현에 더 큰 영향을 줍니다.

`dem_fit`의 LM loop는 source에 명시되어 있어 같은 순서와 자료형으로 옮길 수 있었습니다.
반면 `xrt_iterative`는 외부 MPFIT 구현에 의존하므로 같은 iteration path를 보장하기 어렵습니다.

실제로 **반복 횟수 21회까지** IDL과 일치했습니다.

## 6. 검증 결과와 허용한 차이

### 결과 (set1)

| 항목 | 값 |  |  |
|---|---|---|---|
| 세트 | 합성 기여함수 31T × 10라인 (Gaussian), 스플라인 DEM 노드 [22, 23, 22.5, 21.5, 20.5] |  |  |
|  |  | [22, 23, 22.5, |  |
| 초기값 | 평탄 22  (5노드 전부) |  |  |
| 판정 | 모두 통과 — 15 통과 / 2 차이 허용 / 0 실패 (exit 0) |  |  |
| 회복 결과 | mesh ≈ [22, 23, 22.5, 21.5, 20.5]  — 참 DEM 회복, chisqr 3.9e-6, 21 회 |  |  |
| IDL SPLINE 재구현 | 제어실험 max_rel 1.6e-7 (float32 ULP 수준) |  |  |
| 01_forward, 99_final | 엄격한 1e-5 통과 |  |  |
| 회귀 테스트 | 3/3 |  |  |

|  | 21.5, 20.5] |
|---|---|

### 허용한 차이 2건 — 02_iter1/alpha, 02_iter1/beeta

| 항목 | 내용 |
|---|---|
| 관측 차이 | max_rel 6.3e-5 (alpha) / 4.1e-5 (beeta) |
| 원인 | idl_spline 이 IDL SPLINE과 float32 ULP(~1e-7)까지만 일치 → 유한차분 미분에서 증폭 → 곱셈·누적합에서 재증폭 |
| 영향 | 최종 출력에는 없음. `01_forward`와 `99_final`은 기준을 통과했고, 중간 곡률행렬의 작은 차이는 최종 해에 남지 않았음 |
| 완화기준 | rtol 1e-4 / atol 1e-3 |

**이 허용한 차이는 피할 수 없는 종류입니다.** `idl_spline`은 이미 float32가 표현할 수 있는 한계까지 일치시켰습니다. 그 아래로 내려갈 방법이 없고, 유한차분은 정의상 **작은 차이를 작은 수로 나눠** 상대오차를 키웁니다.

중요한 점은 **최종 DEM·강도·χ²가 정해 둔 기준을 통과했다**는 것입니다. 이 시험에서는 중간값의 작은 차이가 다음 LM 반복에서 줄어 최종 해에 남지 않았습니다.

### 변환 범위

| 항목 | 처리 |
|---|---|
| dem_fit  + dem_functn  + dem_deriv  + dem_chisqr | ✅ 변환·검증 (IDL SPLINE 재구현 포함) |
| run_data2dem_reg | ❌ `demreg` 패키지에서 별도로 다룸 |
| run_mcmc_dem | ❌ `pintofale_mcmc` 패키지에서 별도로 다룸 |
| `mpfit_dem` | ❌ 대안 MPFIT 변형 — 변환하지 않음 |
| get_contributions, calc_dmm_*, change_abund, plot_* | ❌ CHIANTI 데이터 준비 / GUI 입력 준비 코드 |

### 사용 시 주의사항

**노드 개수가 관측 선 개수보다 적어야 합니다.** `nfree ≤ 0` 이면 즉시 실패합니다

노드 **위치**( `log_t_mesh` )는 고정입니다 — 값만 적합됩니다. 위치 선택이 결과를 좌우합니다

LM 가중치가 `1/max(σ, 1)` 이라 **σ가 1보다 작은 선들은 동등하게 취급**됩니다

Spline tension은 `σ = 1.0`인 경우만 검증했습니다.

검증 자료는 **합성 자료 1세트**입니다. 실제 관측 선 목록과 CHIANTI contribution function은 확인하지 않았습니다.

LM은 초기값에 따라 **local minimum**에 수렴할 수 있습니다.

**기록의 근거.** 코드 설명은 `idl/dem_fit.pro`와 `python/chianti_dem_fit.py`를 바탕으로
작성했다. 가중치와 `alpha` 전치에 관한 설명은 코드를 읽고 정리한 것이며, 원저자의 의도를
확인한 내용은 아니다. 검증 수치는 변환 작업 당시의 비교 보고서를 따른다.

</details>

## 참고 자료

- Dere et al. (1997), [CHIANTI — an atomic database for emission lines. I. Wavelengths greater than 50 Å](https://doi.org/10.1051/aas:1997368) — CHIANTI 데이터베이스의 첫 논문
- [CHIANTI 논문 목록](https://www.chiantidatabase.org/chianti_papers.html) — 사용하는 CHIANTI 버전에 맞는 인용 논문 확인용
- `dem_fit.pro` 머리말은 LM 계산 설명으로 Bevington의 *Data Reduction and Error Analysis for the Physical Sciences*, pp. 237–239를 가리킨다. 논문이 아니라 책이다.
