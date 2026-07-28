# chianti_dem (dem_fit) 코드 설명

Del Zanna / Dere — CHIANTI 스플라인 DEM, Levenberg-Marquardt 적합

SpaceAI-DEM / `work/chianti_dem_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — 노드 5개짜리 스플라인을 LM으로 맞춘다
2. `idl_spline`  — IDL SPLINE을 통째로 재구현하다
3. 세 개의 부품 — `functn`  / `deriv`  / `chisqr`
4. Levenberg-Marquardt 루프
5. 왜 이 LM은 parity가 됐나 — xrt와의 결정적 차이
6. 검증 결과와 waiver

## 1. 한눈에 — 노드 5개짜리 스플라인을 LM으로 맞춘다

CHIANTI의 고전적 DEM 적합 루틴입니다. 관측은 **스펙트럼 선의 강도**(이미저 채널이 아님)이고, DEM은 **몇 개의 노드를 지나는 스플라인 곡선**으로 표현합니다.

```
미지수:  log_dem_mesh   = 노드에서의 log DEM 값 (검증에서 5개)
고정:    log_t_mesh     = 노드의 온도 위치
DEM(T) = spline(log_t_mesh, log_dem_mesh) 를 세밀한 온도격자에 펼친 것
```

구조는 xrt_iterative(spline knot + 최적화)와 같은 계열입니다. 하지만 **변환 결과가 정반대로 갈렸습니다** — 그 이유가 이 문서의 5장입니다.

| 단계 | 내용 |
|---|---|
| DEM 표현 | 노드 5개를 지나는 spline-under-tension |
| 순방향 | ∫ G(T)·DEM(T) dT  → 선 강도 예측 |
| 미분 | 유한차분 (노드를 하나씩 1%씩 흔들어봄) |
| 최적화 | Levenberg-Marquardt — .pro  파일 안에 직접 구현됨 |
| 종료 | χ² 개선폭 < 1e-5 또는 20회 |

## 2. idl_spline — IDL SPLINE을 통째로 재구현하다

이 패키지에서 가장 이례적인 부분입니다.

이 저장소의 다른 패키지들은 IDL 라이브러리 함수를 **주입**으로 처리했습니다 ( `aia_get_response` , `int_tabulated` , `gaussian_function` , `eit_line_map`  …).

그런데 `dem_fit` 은 **매 반복마다 SPLINE을 호출**합니다. 주입으로는 해결이 안 되므로, **IDL SPLINE 자체를** **Python으로 완전히 재구현**했습니다 — 45줄짜리 `idl_spline`  함수가 그것입니다.

### 보통의 스플라인이 아닙니다

IDL의 `SPLINE` 은 3차 스플라인이 아니라 **spline-under-tension**(장력 스플라인, 기본 `σ = 1.0` )입니다.

```
dels   = sigmap * delx1
exps   = np.exp(dels)
sinhs  = 0.5 * (exps - 1./exps)          ← sinh
diag1  = sinhin * (dels*0.5*(exps+1./exps) - sinhs)   ← cosh 포함
```

코드 전반에 `sinh` **·** `cosh` 가 등장하는 이유입니다. 3차 다항식 대신 **쌍곡선 함수**를 씁니다. 장력 파라미터 `σ` 가 커질수록 곡선이 직선에 가까워지고, 0에 가까우면 일반 3차 스플라인이 됩니다.

구조는 표준 스플라인과 같습니다:

1. **양 끝점 기울기 추정** (14–20줄) — 3점 유한차분 공식

2. **삼중대각 시스템 풀기** (33–39줄) — 전진 소거 후 후진 대입

3. **보간값 계산** (40–45줄) — `searchsorted` 로 구간을 찾아 쌍곡선 보간

> **⚠ 전체가** `float32` **입니다.** 파일 맨 위의 `f32 = np.float32` 가 코드 곳곳에 박혀 있는 이유입니다. IDL `SPLINE` 이 float32로 계산하므로, 중간값 하나라도 float64로 새면 결과가 어긋납니다. 제어실험에서 **max_rel 1.6e-7**로 일치를 확인했는데, 이는 **float32의 ULP(최소 단위) 수준**입니다 — 즉 **가능한 한계까지** **일치**시킨 것이고, 이보다 더 좁힐 수는 없습니다.

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

> **⚠ 해석적 미분이 아니라 유한차분입니다.** 노드 하나당 `functn()` 을 한 번씩 더 부르므로, 반복마다 **스플라인 계산이 (노드 수 + 1)회** 일어납니다.

그리고 이것이 6장 waiver의 원인입니다 — `idl_spline` 의 **float32 ULP 수준 차이**가 **작은 수로 나누는** **유한차분**에서 증폭됩니다.

### chisqr  82–84

```
dy = s.obs_int - s.functn()
return np.float32(np.sum(dy ** 2 / s.obs_sig ** 2))
```

표준 χ²입니다. **reduced가 아닌 총합**이라는 점에 주의하세요.

## 4. Levenberg-Marquardt 루프 86–148

Bevington의 고전적 `CURFIT`  구조를 따릅니다.

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

`alpha` 는 simple_reg_dem에서 본 `JᵀJ` 와 같은 역할이고, `beeta` 는 `Jᵀ·잔차`입니다.

> **⚠ 원본의 특징 두 가지** (변환본은 그대로 재현했습니다):

① **가중치가** `1/σ` **입니다** — `chisqr` 가 쓰는 `1/σ²` 와 다릅니다. 게다가 `max(σ, 1.0)`  바닥값이 걸려 있어, **σ < 1인 선은 가중치가 1로 고정**됩니다. LM은 하강 방향만 맞으면 수렴하므로 동작에는 문제가 없지만, 교과서적 Gauss-Newton은 아닙니다.

② `alpha` 를 **하삼각만 채운 뒤 transpose**합니다. 대칭화가 아니라 전치이므로 결과는 **상삼각 행렬**이 되고 하삼각은 0으로 남습니다.

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

`flambda` **가 LM의 핵심입니다.** 대각선에 `1 + flambda` 를 넣으면:

• `flambda` 가 작으면 → Gauss-Newton에 가까움 (과감한 큰 스텝)

• `flambda` 가 크면 → 경사하강에 가까움 (안전한 작은 스텝)

개선되면 과감하게( `/10` ), 실패하면 보수적으로( `×10` ) 전환하며 **스스로 보폭을 조절**합니다. simple_reg_dem의 trust region + line search와 목적이 같습니다.

`NaN`  처리(123–128줄)와 `flambda > 1e10`  탈출(135–136줄)이 무한루프 방지 장치입니다.

### 종료 조건 99, 143

```
while (it <= niter) and (dchisq >= dchisq_min) and not failed:
    ...
    dchisq = chisq1 - chisqr      # 이번 반복의 χ² 개선폭
```

χ² 개선폭이 `1e-5`  미만이 되거나 20회를 넘으면 종료합니다. `nfree = n_obs − nterms ≤ 0` 이면(관측보다 노드가 많으면) **시작도 안 하고 실패 처리**합니다(96–97줄).

## 5. 왜 이 LM은 parity가 됐나 — xrt와의 결정적 차이

xrt_iterative에서는 MPFIT(Levenberg-Marquardt)을 **변환 경계**로 두고 제외했습니다. 그런데 여기서는 같은 LM인데 **ALL PASS**입니다. 모순처럼 보이지만 이유가 명확합니다.

|  | xrt_iterative | chianti_dem |
|---|---|---|
| LM의 소재 | MPFIT 라이브러리 호출 | dem_fit.pro  안에 직접 작성 |
| 변환 방법 | scipy로 대체해야 함 | 한 줄씩 그대로 포팅 |
| 경로 | 구현이 다르니 다른 길로 감 | 같은 연산 순서 → 같은 길 |
| 결과 | parity 불가 → 경계 | parity 성립 |

**경로 의존성 자체가 문제인 게 아닙니다.** 문제는 **구현을 갈아끼울 때** 생깁니다.

알고리즘이 소스에 그대로 적혀 있으면, 아무리 경로 의존적이어도 **같은 순서로 같은 연산을 재현**해 동일한 경로를 밟을 수 있습니다. 라이브러리 블랙박스를 다른 라이브러리로 바꾸는 순간 그게 불가능해집니다.

실제로 **반복 횟수 21회까지** IDL과 일치했습니다.

## 6. 검증 결과와 waiver

### 결과 (set1)

| 항목 | 값 |  |  |
|---|---|---|---|
| 세트 | 합성 기여함수 31T × 10라인 (Gaussian), 스플라인 DEM 노드 [22, 23, 22.5, 21.5, 20.5] |  |  |
|  |  | [22, 23, 22.5, |  |
| 초기값 | 평탄 22  (5노드 전부) |  |  |
| 판정 | ALL PASS — 15 pass / 2 waived / 0 fail (exit 0) |  |  |
| 회복 결과 | mesh ≈ [22, 23, 22.5, 21.5, 20.5]  — 참 DEM 회복, chisqr 3.9e-6 , 21 회 |  |  |
| IDL SPLINE 재구현 | 제어실험 max_rel 1.6e-7 (float32 ULP 수준) |  |  |
| 01_forward , 99_final | strict 1e-5 통과 |  |  |
| 회귀 테스트 | 3/3 |  |  |

|  | 21.5, 20.5] |
|---|---|

### waiver 2건 — 02_iter1/alpha , 02_iter1/beeta

| 항목 | 내용 |
|---|---|
| 관측 차이 | max_rel 6.3e-5 (alpha) / 4.1e-5 (beeta) |
| 원인 | idl_spline 이 IDL SPLINE과 float32 ULP(~1e-7)까지만 일치 → 유한차분 미분에서 증폭 → 곱셈·누적합에서 재증폭 |
| 영향 | 없음. 순방향( 01_forward )과 최종 적합( 99_final )은 strict 통과. LM이 스스로 보정하므 로 중간 곡률행렬의 미세차가 최종해에 남지 않음 |
| 완화기준 | rtol 1e-4 / atol 1e-3 |

**이 waiver는 피할 수 없는 종류입니다.** `idl_spline` 은 이미 float32가 표현할 수 있는 한계까지 일치시켰습니다. 그 아래로 내려갈 방법이 없고, 유한차분은 정의상 **작은 차이를 작은 수로 나눠** 상대오차를 키웁니다.

중요한 것은 **인증 대상 출력(최종 DEM·강도·χ²)이 strict 기준을 통과**했다는 사실입니다. LM처럼 **자기보정하는 반복**은 중간값의 미세한 흔들림을 흡수합니다.

### 스코프 — 무엇을 변환하고 무엇을 뺐나

| 항목 | 처리 |
|---|---|
| dem_fit  + dem_functn  + dem_deriv  + dem_chisqr | ✅ 변환·검증 (IDL SPLINE 재구현 포함) |
| run_data2dem_reg | ❌ 중복 — demreg와 동일. work/demreg_parity/ 에서 이미 인증 |
| run_mcmc_dem | ❌ MCMC — 확률적이라 bit-parity 불가 ( work/SPECIAL_CASES.md ) |
| mpfit_dem | ❌ 대안 MPFIT 변형 — 미변환 |
| get_contributions , calc_dmm_* , change_abund , plot_* | ❌ CHIANTI 데이터 준비 / GUI 스캐폴딩 |

### 알고 써야 할 것

**노드 개수가 관측 선 개수보다 적어야 합니다.** `nfree ≤ 0` 이면 즉시 실패합니다

노드 **위치**( `log_t_mesh` )는 고정입니다 — 값만 적합됩니다. 위치 선택이 결과를 좌우합니다

LM 가중치가 `1/max(σ, 1)` 이라 **σ가 1보다 작은 선들은 동등하게 취급**됩니다

스플라인 장력 `σ = 1.0` 으로만 검증했습니다

검증은 **합성 데이터 1세트**입니다. 실관측 선목록과 실제 CHIANTI 기여함수는 후속 과제

초기값에 따라 **지역 최소점**에 빠질 수 있습니다 (LM의 일반적 성질)

**출처 표시.** 2~4장의 코드 설명은 `converted/chianti_dem_fit.py` 와 `inbox/dem_fit.pro` 에서, 6장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 4장의 가중치· `alpha`  전치에 대한 지적은 변환된 코드를 읽고 기술한 것으로, 원저자의 의도나 이론적 근거를 확인한 것은 아닙니다. 5장의 xrt 대비와 1장의 방법 분류는 두 인증서와 코드에 근거한 해석입니다.
