# firdem 코드 설명

Plowman et al. 2012 — FIR DEM (Fast, Iterative, Regularized)

SpaceAI-DEM / `work/firdems_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — simple_reg_dem의 조상
2. 삼각 기저 — DEM을 텐트로 표현하기
3. `firdem_regularize_data`  — 해가 아니라 **데이터**를 정규화한다
4. `first_pass_dem`  — 한 줄로 끝나는 역산
5. 변환 범위 — 무엇을 옮기고 무엇을 남겼나
6. IDL ↔ Python 대조
7. 검증 결과

## 1. 한눈에 — simple_reg_dem의 조상

이름의 **FIR**은 **F**ast, **I**terative, **R**egularized입니다. 같은 저자(Plowman)가 8년 뒤 내놓은 것이 simple_reg_dem이고, `work/DEM_SOLVER_METHODS.md` 는 후자를 **"이 방법의 후대 단일파일 버전"**으로 기록하고 있습니다.

두 코드를 나란히 놓으면 무엇이 계승되고 무엇이 버려졌는지 보입니다:

|  | firdem (2012) | simple_reg_dem (2020) |
|---|---|---|
| DEM 표현 | 삼각 기저함수 전개 | 온도격자 free-form (내부적으로 같은 삼각 형) |
| 정규화 대상 | 데이터를 흔든다 | 해에 페널티를 건다 |
| 역산 | SVD 유사역 한 번 | Cholesky 반복 |
| 양수성 | 2000회 반복으로 음수 제거 | log 치환으로 공짜 |
| 파일 규모 | 826줄 + 라이브러리 의존 다수 | 86줄, 의존성 0 |

**simple_reg_dem의 "simple"이 무엇에 대한 대답인지** 이 표가 보여줍니다. firdem에서 가장 골치 아팠던 두 가지 — 2000회 양수화 반복과 라이브러리 의존성 — 을 **log 공간 치환 하나로 동시에 없앤 것**이 후속작입니다.

## 2. 삼각 기저 — DEM을 텐트로 표현하기

```
def firdem_triangle(x, x0, b, h):
    sep = np.abs(x - x0)
    y = h * (0.5 * b - sep) / (0.5 * b)
    y = np.where(y < 0, 0.0, y)
    return y
```

`x0` 를 꼭짓점으로, 밑변 폭 `b` , 높이 `h` 인 **삼각형**입니다. 바깥은 0으로 자릅니다.

```
def firdem_triangle_basis(nb, logt):
    logdt = (logthi - logtlo) / (nb + 1)              # 기저 간격
    t0a   = logdt + logtlo + logdt * np.arange(nb)    # 꼭짓점 위치들
    for j in range(nb):
        basis[:, j] = firdem_triangle(logt, t0a[j], 2.0 * logdt, 1.0)
    return basis, t0a
```

밑변이 `2·logdt` 이고 꼭짓점 간격이 `logdt` 이므로, **이웃한 삼각형이 정확히 절반씩 겹칩니다.**

```
      /\    /\    /\    /\
     /  \  /  \  /  \  /  \        각 지점에서 두 삼각형의 합 = 1
    /    \/    \/    \/    \
   ─────┴──────┴──────┴─────  logT
      t0a[0] t0a[1] t0a[2]
```

이것이 **선형 유한요소의 hat function**입니다. simple_reg_dem에서 `Bij` (FEM mass matrix)가 암묵적으로 쓰던 바로 그 표현을, firdem은 **명시적인 기저 행렬**로 들고 다닙니다.

표현 방식이 같으므로 두 코드가 그리는 DEM 곡선의 형태(격자점 사이를 직선으로 연결)도 같습니다.

`logdt = (logthi − logtlo)/(nb+1)` 에서 `nb+1` **로 나누는 것**에 주목하세요. 꼭짓점이 양 끝에서 `logdt` 만큼 안쪽으로 들어가므로, **경계에서 삼각형이 잘리지 않습니다.**

## 3. firdem_regularize_data — 해가 아니라 데이터를 정규화한다

이 함수가 firdem의 가장 독특한 아이디어입니다.

보통의 정규화는 **해**에 "매끄러워라"는 페널티를 겁니다 (simple_reg_dem의 `regmat` ).

firdem은 반대로 **데이터를 잡음 허용 범위 안에서 살짝 흔들어**, 그 흔든 데이터를 그냥 유사역으로 풀면 얌전한 답이 나오도록 만듭니다.

### 조기 반환 두 가지 51–54

```
if np.sum((datavec0 / sigmas) ** 2) < chi2_end:
    return np.zeros(nchan)          # ① 데이터가 잡음 수준 → DEM 없음
if chi2_end / nchan < 1.0e-4:
    return datavec0                 # ② 정규화 예산이 0 → 원본 그대로
```

`chi2_end` 는 **"데이터를 얼마나 흔들어도 되는가"의 예산**입니다 (오라클에서 `chisqr_cvf` 로 계산해 주입). ①은 데이터 자체가 그 예산보다 작아 통째로 잡음이라는 뜻이고, ②는 예산이 없어 손댈 수 없다는 뜻입니다.

### 반복 루프 64–87

```
for k in range(niter_max):                                   # 최대 50회
    datavec  = datavec0 / tr_norm                            # ① 매번 원본에서 다시 시작
    data2vec = (a_inv @ datavec) / tr_norm
    if k == 0:
        alpha = np.sqrt(9.0*chi2_end / np.sum(sigmas**2 * data2vec**2))   # ② 초기 추정
    b_inv   = np.linalg.inv(alpha * a_inv_scaled + sigs2_diag)            # ③ 행렬 역산
    datavec = datavec - (alpha * (b_inv @ data2vec)) / tr_norm            # ④ 데이터 수정
    chi2    = np.sum(((datavec0 - datavec*tr_norm) / sigmas) ** 2)        # ⑤ 얼마나 흔들었나
```

**① 누적이 아니라 매번 새로 계산**

루프 첫 줄에서 `datavec` 을 **원본으로 되돌립니다.** 즉 이 반복은 답을 조금씩 개선하는 게 아니라, `alpha` **값 하나를 찾는 탐색**입니다. 매 회차가 독립적인 시도예요.

**③ Tikhonov 형태의 축소(shrinkage)**

```
b_inv = inv( α · A_scaled + diag(1/σ²) )
```

`α` 가 크면 정규화가 세지고, 0이면 데이터가 그대로 남습니다. simple_reg_dem에서 `regmat` 이 `JᵀJ` 에 더해지던 것과 **수식 모양이 같습니다** — 다만 작용 대상이 해가 아니라 데이터입니다.

**⑤** `chi2` **의 의미가 다릅니다**

```
chi2 = Σ ( (원본데이터 − 수정된데이터) / σ )²
```

> **⚠ 여기서** `chi2` **는 "모델이 데이터에 얼마나 맞나"가 아닙니다.** **"데이터를 원본에서 얼마나 멀리 옮겼** **나"**입니다. 다른 solver의 χ²와 의미가 정반대이므로 읽을 때 주의해야 합니다.

### alpha 탐색 — 5배 증가 후 이분법 73–87

```
if bisect_start == 1:                       # 이분법 단계
    if chi2 < chi2_end:  alpha_low  = alpha
    if chi2 > chi2_end:  alpha_high = alpha
    alpha = alpha_low + 0.5*(alpha_high - alpha_low)
if chi2 < chi2_end and bisect_start == 0:   # 아직 여유 있음 → 과감히 키움
    alpha = 5 * alpha
if chi2 > chi2_end and bisect_start == 0:   # 처음 넘어섬 → 이분법 시작
    bisect_start = 1
    alpha_high = alpha
    alpha = alpha_low + 0.5*(alpha_high - alpha_low)
if abs(chi2 - chi2_end)/chi2_end < chi2_tol:   # 5% 이내면 종료
    break
```

2단계 탐색입니다:

1. **확장 단계** — 예산을 아직 안 썼으면 `alpha` 를 **5배씩** 키워 빠르게 범위를 찾습니다

2. **이분 단계** — 예산을 넘는 순간 상한이 확정되고, 그때부터 **구간을 반씩 좁혀** 정확히 `chi2_end` 에 착지합니다

**목표는 χ²를 최소화하는 게 아니라 정확히** `chi2_end` **에 맞추는 것**입니다. "잡음이 허용하는 만큼만 데이터를 흔들고, 그 이상은 흔들지 않는다" — 개념적으로 **discrepancy principle**이고, simple_reg_dem의 χ²=1 targeting과 같은 철학입니다.

## 4. first_pass_dem — 한 줄로 끝나는 역산

```
def first_pass_dem(data_out, a_inv, normfac, basis22):
    return basis22 @ (a_inv @ (data_out / normfac))
```

정규화된 데이터가 준비되면, DEM 복원 자체는 **행렬 곱 두 번**입니다.

| 단계 | 하는 일 |
|---|---|
| / normfac | 채널별 정규화 — 스케일 맞추기 |
| a_inv @ | SVD 유사역. 데이터 → 저해상도 기저 계수 |
| basis22 @ | spline 보간 행렬. 저해상도 계수 → 고해상도 기저 계수 |

> **⚠** `a_inv` **는 이 코드가 만들지 않습니다.** 응답함수의 spline 보간, `int_tabulated`  기저 정규화, `svdc`  특이값 분해를 거쳐 만들어지는데, 전부 IDL 라이브러리 함수라 **오라클에서 주입**했습니다. 검증 대상은 **이** **곱셈 산술 자체**입니다.

여기까지가 "first pass"이고, 이 결과 `dem_initial` 에는 **음수 구간이 섞여 있을 수 있습니다.** 원본은 이후 `firdem_iterate` 로 그걸 제거하는데 — 그게 다음 장의 주제입니다.

## 5. 변환 범위 — 무엇을 옮기고 무엇을 남겼나

인증서의 판정은 **"PASS (scoped)"**입니다. 전체가 아니라 **결정론적으로 검증 가능한 핵심**만 변환했습니다.

### ✅ 변환·검증한 것

`firdem_regularize_data`  — 반복 Tikhonov 데이터 정규화 (double 정밀도 행렬 역산)

`first_pass_dem`  — `basis22 # (a_inv # (data/normfac))`

`firdem_triangle_basis`  — 고해상도 삼각 기저

### 🔌 주입한 라이브러리 경계 (재유도 안 함)

`a_struc`  — 응답 spline 보간( `spl_interp` ), 기저· `a_array`  정규화( `int_tabulated` ), `svdc`  → `a_inv` `basis22`  ( `spl_interp` ), `a2_array`  ( `int_tabulated` ), `chisqr_cvf`

### ❌ 변환하지 않은 것 — firdem_iterate

음수 EM을 선형 외삽으로 제거하며 χ² 개선 조건에서 **최대 2000회** 반복하는 루틴입니다.

**왜 제외했나.** 이 루프는 **경로 의존적**입니다. 매 반복의 `break`  조건과 외삽 방향이 직전 값에 의존하므로, 부동소수점 차이가 **반복 횟수와 종료 지점 자체**를 바꿉니다. 게다가 앞단의 라이브러리 함수

( `spl_interp` , `svdc` )가 bit-exact로 재현되지 않는 이상 입력부터 미세하게 달라집니다. 따라서 최종 `dem_out.coffs` 의 strict parity는 **원리적으로 보장 불가**이며, 별도 후속 대상으로 문서화했습니다.

**이 판단은 이 저장소 전반의 원칙입니다.** `work/CONVERSION_SUMMARY.md` 에 **"경로의존 반복/최적화는 bit-parity** **불가"**로 정리돼 있고, 같은 이유로 **xrt_iterative(MPFIT)**, **vdem(GCV/Brent)**, **pintofale(MCMC)**도 최적화 루프를 경계로 두었습니다. 대응 방식은 공통 — **solver가 최적화하는 물리(순방향모델 + χ²)를 검증하고, 최적화 루프는** **문서화된 경계로 남긴다.**

## 6. IDL ↔ Python 대조

| IDL | Python |
|---|---|
| M # v  (2D # 1D) | M @ v |
| a # b  (1D # 1D) | np.outer(a, b) |
| invert(A, /double) | np.linalg.inv(A) |
| diag_matrix(1/sigmas)^2 | np.diag((1/sigmas)**2) |
| 원소별 > | np.maximum |
| findgen(nb) | np.arange(nb, dtype=np.float64) |
| y < 0  인덱싱 후 0 대입 | np.where(y < 0, 0.0, y) |

## 7. 검증 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일 픽셀 set1 — 합성 Gaussian DEM, aia_get_response  응답 |
| 판정 | ALL PASS — 21 pass / 0 waived / 0 fail ( compare_probes  exit 0) |
| 03_regdata | 정규화된 데이터 — rtol 1e-8 통과 |
| 04_firstpass | dem_initial  (~1e28) — rtol 1e-8 통과 |
| 02_setup | basis2 / t0a  — atol 2e-5 통과 (아래 설명) |
| 회귀 테스트 | 5/5 PASS |

### float32 그리드 문제 — waiver가 아닌 이유

`basis2` 에서 약 `8.5e-6` 의 차이가 관측되어 `atol 2e-5`  기준을 적용했습니다.

원인은 `a_struc.t` **(온도축)가 IDL float32**라는 데 있습니다. 삼각형의 **가장자리**에서 `0.5b − |x − x0|` 가 0에 가까워지는데, 이때 float32 온도값의 반올림 오차가 그대로 드러납니다.

즉 **계산이 틀린 게 아니라 입력 격자의 정밀도 한계**입니다. 그래서 인증서는 이를 **waiver가 아니라 "IDL** **데이터형에 맞춘 정책"**으로 분류했습니다 — firdem의 최종 waiver 건수는 **0**입니다.

### 한계

검증 입력은 **합성 Gaussian 단일 픽셀** 1세트입니다

`firdem_iterate` 를 거치지 않았으므로, 이 코드의 출력 `dem_initial` 에는 **음수가 남아 있을 수 있습니다.** 최종 DEM이 아닙니다

맵 단위 처리, 실관측 다세트는 후속 과제

`a_inv`  등 주입값이 바뀌면(다른 응답함수·기저 개수) 재검증이 필요합니다

**출처 표시.** 2~4·6장의 코드 설명은 `converted/firdem.py` 와 그 docstring에서, 5·7장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1장의 simple_reg_dem과의 비교는 `work/DEM_SOLVER_METHODS.md` · `work/ROADMAP.md` 의 기록과 두 코드의 대조에 근거한 해석이며, Plowman et al.2012 원논문 대조는 수행하지 않았습니다.
