# firdems — FIR DEM의 첫 계산 단계

FIR은 Fast, Iterative, Regularized의 머리글자다. DEM을 여러 삼각형 모양의 합으로 나타내고,
관측 오차 범위 안에서 자료를 안정화한 뒤 첫 DEM을 계산한다.

## 먼저 읽기

이번 Python 변환에는 FIR 전체 과정이 들어 있지 않다. 포함된 부분은 삼각형 후보 모양 만들기,
관측자료 정규화, 그리고 첫 번째 DEM 계산이다. 원본의 최종 반복 함수 `firdem_iterate`는 옮기지
않았다.

따라서 Python의 `dem_initial`은 **최종 FIR DEM이 아니다.** 음수 값이 남을 수 있으며, 원래
알고리즘은 뒤의 반복 단계에서 이를 더 다듬는다. 이 차이를 모르고 결과를 최종 DEM으로 사용하면
안 된다.

합성 Gaussian DEM 한 세트에서 21개 비교 항목과 회귀 테스트 5개를 통과했다. 응답행렬과
역행렬 등은 외부에서 계산한 값을 입력해 시험했다. 다른 기기나 다른 온도 후보 수를 쓰면 그
입력값을 다시 만들고 결과를 확인해야 한다.

- IDL 전체 배포본과 원문 안내: `idl/firdems/`
- Python으로 옮긴 일부 기능: `python/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. `simple_reg_dem`과의 관계
2. 삼각 기저 — DEM을 텐트로 표현하기
3. `firdem_regularize_data` — 해가 아니라 **데이터**를 정규화한다
4. `first_pass_dem` — 한 줄로 끝나는 역산
5. 변환 범위 — 무엇을 옮기고 무엇을 남겼나
6. IDL ↔ Python 대조
7. 검증 결과

## 1. `simple_reg_dem`과의 관계

이름의 **FIR**은 **F**ast, **I**terative, **R**egularized입니다. 같은 저자 Plowman은 약 8년 뒤
`simple_reg_dem`을 발표했습니다. 두 코드는 정규화에 같은 생각을 일부 공유하지만 구현과 출력은
서로 다릅니다.

두 코드를 나란히 놓으면 무엇이 계승되고 무엇이 버려졌는지 보입니다:

|  | firdem (2012) | simple_reg_dem (2020) |
|---|---|---|
| DEM 표현 | 삼각 기저함수 전개 | 온도격자 free-form (내부적으로 같은 삼각 형) |
| regularization 대상 | data를 noise 범위에서 조정 | 해에 penalty를 적용 |
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

`x0`를 꼭짓점으로, 밑변 폭 `b`, 높이 `h` 인 **삼각형**입니다. 바깥은 0으로 자릅니다.

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

`logdt = (logthi − logtlo)/(nb+1)`에서 `nb+1` **로 나누는 것**에 주목하세요. 꼭짓점이 양 끝에서 `logdt` 만큼 안쪽으로 들어가므로, **경계에서 삼각형이 잘리지 않습니다.**

## 3. firdem_regularize_data — 해가 아니라 데이터를 정규화한다

이 함수가 firdem의 가장 독특한 아이디어입니다.

일반적인 regularization은 해의 roughness에 penalty를 줍니다(`simple_reg_dem`의 `regmat`).

FIR DEM은 data를 noise 허용 범위 안에서 조정한 뒤 pseudoinverse를 적용했을 때 안정적인 해가
나오도록 만듭니다.

### 조기 반환 두 가지 51–54

```
if np.sum((datavec0 / sigmas) ** 2) < chi2_end:
    return np.zeros(nchan)          # ① 데이터가 잡음 수준 → DEM 없음
if chi2_end / nchan < 1.0e-4:
    return datavec0                 # ② 정규화 예산이 0 → 원본 그대로
```

`chi2_end`는 data를 원래 값에서 얼마나 조정할 수 있는지를 정하는 한계입니다. IDL에서는
`chisqr_cvf`로 계산하며 이번 검증에서는 그 값을 외부 입력으로 사용했습니다. 첫 번째 조기
종료는 signal이 noise와 구분되지 않는 경우이고, 두 번째는 허용된 조정 폭이 0인 경우입니다.

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

각 반복의 첫 줄에서 `datavec`을 원본으로 되돌립니다. 따라서 이 loop는 이전 해를 갱신하는
과정이 아니라 `alpha`를 찾는 독립적인 trial의 반복입니다.

**③ Tikhonov 형태의 축소(shrinkage)**

```
b_inv = inv( α · A_scaled + diag(1/σ²) )
```

`α`가 크면 정규화가 세지고, 0이면 데이터가 그대로 남습니다. simple_reg_dem에서 `regmat`이 `JᵀJ`에 더해지던 것과 **수식 모양이 같습니다** — 다만 작용 대상이 해가 아니라 데이터입니다.

**⑤** `chi2` **의 의미가 다릅니다**

```
chi2 = Σ ( (원본데이터 − 수정된데이터) / σ )²
```

> 여기서 `chi2`는 model fit의 χ²가 아니라 **조정한 data가 원래 data에서 얼마나 멀어졌는지**를
> 나타냅니다. 다른 solver의 χ²와 정의가 다릅니다.

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

1. **확장 단계** — 예산을 아직 안 썼으면 `alpha`를 **5배씩** 키워 빠르게 범위를 찾습니다

2. **이분 단계** — 허용 한계를 넘으면 상한을 정하고 구간을 절반씩 줄여 `chi2_end`에 수렴합니다.

목표는 χ²를 최소화하는 것이 아니라 `chi2_end`에 맞추는 것입니다. 즉 noise가 허용하는
범위까지만 data를 조정하는 discrepancy principle이며, `simple_reg_dem`의 χ²≈1 기준과
비슷한 역할을 합니다.

## 4. `first_pass_dem` — 한 줄로 끝나는 역산

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

> **⚠** `a_inv` **는 이 코드가 만들지 않습니다.** 응답함수의 spline 보간, `int_tabulated` 기저 정규화, `svdc` 특이값 분해를 거쳐 만들어지는데, 전부 IDL 라이브러리 함수라 **IDL에서 계산한 값을 입력**했습니다. 검증 대상은 **이** **곱셈 산술 자체**입니다.

여기까지가 first pass이고, 이 결과인 `dem_initial`에는 **음수 구간이 남을 수 있습니다.** 원본은 이후 `firdem_iterate`에서 이를 줄여 나갑니다. 이 부분은 다음 장에서 설명합니다.

## 5. 변환 범위 — 무엇을 옮기고 무엇을 남겼나

비교 보고서의 판정은 **"PASS (범위 제한)"**입니다. 전체가 아니라 **결정론적으로 검증 가능한 핵심**만 변환했습니다.

### ✅ 변환·검증한 것

`firdem_regularize_data` — 반복 Tikhonov 데이터 정규화 (double 정밀도 행렬 역산)

`first_pass_dem` — `basis22 # (a_inv # (data/normfac))`

`firdem_triangle_basis` — 고해상도 삼각 기저

### 🔌 외부 입력으로 남긴 라이브러리 부분 (재유도 안 함)

`a_struc` — 응답 spline 보간( `spl_interp` ), 기저· `a_array` 정규화( `int_tabulated` ), `svdc` → `a_inv` `basis22` ( `spl_interp` ), `a2_array` ( `int_tabulated` ), `chisqr_cvf`

### ❌ 변환하지 않은 것 — firdem_iterate

음수 EM을 선형 외삽으로 제거하며 χ² 개선 조건에서 **최대 2000회** 반복하는 루틴입니다.

이 loop는 path-dependent합니다. 각 반복의 `break` 조건과 extrapolation 방향이 직전 값에
의존하므로 작은 부동소수점 차이도 반복 횟수와 종료점을 바꿀 수 있습니다. 또한 앞단의 라이브러리 함수

`spl_interp`와 `svdc`의 결과가 비트 단위까지 같지 않으면 반복 계산의 입력부터 조금씩
달라집니다. 따라서 최종 `dem_out.coffs`의 비트 단위 일치는 보장하지 않으며, 이번 검증
범위에서도 제외했습니다.

반복 최적화는 아주 작은 반올림 차이만 있어도 다음 계산 순서가 달라질 수 있습니다. 그래서 이
저장소에서는 실제로 옮긴 forward model과 χ² 부분을 확인하고, 옮기지 않은 반복 계산은 각
README에 분명히 적었습니다. `xrt_iterative`의 MPFIT과 `vdem`의 GCV/Brent도 같은 기준을
적용했습니다.

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
| 판정 | 모두 통과 — 21개 일치 / 0개 차이 허용 / 0개 실패 |
| 03_regdata | 정규화된 데이터 — rtol 1e-8 통과 |
| 04_firstpass | `dem_initial`(~1e28) — rtol 1e-8 통과 |
| 02_setup | basis2 / t0a  — atol 2e-5 통과 (아래 설명) |
| 회귀 테스트 | 5/5 PASS |

### float32 그리드 문제 — 허용한 차이가 아닌 이유

`basis2`에서 약 `8.5e-6`의 차이가 관측되어 `atol 2e-5` 기준을 적용했습니다.

원인은 `a_struc.t` **(온도축)가 IDL float32**라는 데 있습니다. 삼각형의 **가장자리**에서 `0.5b − |x − x0|`가 0에 가까워지는데, 이때 float32 온도값의 반올림 오차가 그대로 드러납니다.

즉 **계산이 틀린 게 아니라 입력 격자의 정밀도 한계**입니다. 그래서 비교 보고서는 이를 **허용한 차이가 아니라 "IDL** **데이터형에 맞춘 정책"**으로 분류했습니다 — firdem의 최종 허용한 차이 건수는 **0**입니다.

### 한계

검증 입력은 **합성 Gaussian 단일 픽셀** 1세트입니다

`firdem_iterate`를 거치지 않았으므로 `dem_initial`에는 음수가 남을 수 있습니다. 이 값은 최종 FIR DEM이 아닙니다.

맵 단위 처리와 여러 실관측 자료는 검증하지 않았습니다.

응답함수나 basis 수가 달라져 `a_inv` 등의 외부 입력값이 바뀌면 다시 검증해야 합니다.

**기록의 근거.** 코드 설명은 `idl/firdems/`의 원본과 `python/firdem.py`를 바탕으로 작성했다.
`simple_reg_dem`과의 비교는 두 구현을 직접 대조해 정리했다. 검증 수치는 변환 작업 당시의
비교 보고서를 따른다.

</details>

## 참고 논문

- Plowman, Kankelborg & Martens (2013), [Fast Differential Emission Measure Inversion of Solar Coronal Data](https://doi.org/10.1088/0004-637X/771/1/2) ([arXiv](https://arxiv.org/abs/1204.6306))

원본 `idl/firdems/readme.txt`는 2012년의 예비 배포본이지만, 논문은 2013년에 출판됐다.
