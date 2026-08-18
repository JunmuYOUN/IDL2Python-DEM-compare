# trace_dem — TRACE 3채널용 단순 선형 DEM

TRACE에는 이 계산에서 쓰는 채널이 세 개뿐이다. 이 코드는 DEM의 모양도 미리 정한 세 성분으로
줄여서, 채널 밝기 세 개로 계수 세 개를 한 번에 계산한다.

## 먼저 읽기

계산은 빠르고 단순하다. 노출시간으로 밝기를 보정한 뒤 미리 만든 3×3 역행렬을 곱하고, 세
계수를 이용해 DEM과 평균 온도를 계산한다. 반복 계산이나 정규화는 없다.

단순한 만큼 제한도 분명하다.

- 미리 정한 세 모양 밖의 DEM은 표현할 수 없다.
- 계수가 음수가 되지 않도록 막는 조건이 없다.
- χ²나 잔차를 내놓지 않아 결과가 관측에 잘 맞는지 코드 안에서 판단할 수 없다.
- 원본 주석에 적힌 pedestal은 입력 전에 빼야 한다.

13개 비교 항목과 회귀 테스트 4개를 통과했다. TRACE 응답함수와 적분표를 만드는 과정은 옮기지
않았고, 미리 계산한 `K_inv`와 `tau`를 입력해 시험했다. 문헌 제목으로 확인된 방법 논문은 아직
찾지 못했으며, `Kankelborg 1998`은 원본 코드의 작성자와 날짜 정보다.

- IDL 원본: `idl/`
- Python 변환본: `python/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 방법 요약 — 세 개의 basis coefficient 계산
2. 기저함수 3개 — 이 방법의 전부
3. 코드 흐름
4. IDL ↔ Python 대조
5. 변환 및 검증 범위

## 1. 방법 요약 — 세 개의 basis coefficient 계산

TRACE는 EUV 채널이 **3개뿐**입니다 (171, 195, 284 Å). AIA의 6채널보다도 적습니다.

채널이 세 개뿐이므로 자유로운 DEM을 복원하면 매우 ill-posed합니다. 이 구현은 DEM 자체의
자유도를 세 개로 제한합니다.

**미지수를 아예 3개로 줄인다.**

DEM을 온도격자 위 자유 함수로 두지 않고, **미리 정한 기저함수 3개의 선형결합**으로 표현합니다. 그러면 관측 3개 · 미지수 3개 → **정정(well-determined)**. 답이 딱 하나로 결정됩니다.

그 결과 이 코드에는 다른 DEM solver에 반드시 있는 것들이 **전부 없습니다**:

| 보통 DEM solver | trace_dem |
|---|---|
| 정규화 항 (regularization) | 없음 — 기저를 3개로 줄인 것 자체가 정규화 |
| 반복 (iteration) | 없음 — 행렬 한 번 곱하면 끝 |
| 양수성 강제 | 없음 — 계수가 음수로 나올 수 있음 |
| χ² 판정 | 없음 — 방정식 수 = 미지수 수라 잔차가 0 |

DEM 복원의 핵심 계산이 **3×3 행렬 곱셈 한 줄**입니다:

```
emc = K_inv @ p
```

대신 대가가 있습니다. **DEM의 모양을 미리 가정**했으므로, 실제 플라즈마가 그 3개 기저로 표현되지 않는 모양이면 복원할 수 없습니다.

## 2. 기저함수 3개 — 이 방법의 전부

DEM 곡선을 이렇게 표현합니다:

```
DEM(T) = emc₀·A₀·sbasis(0,T) + emc₁·A₁·sbasis(1,T) + emc₂·A₂·sbasis(2,T)
         └── 픽셀마다 구하는 값 ──┘  └─ 고정된 모양 ─┘
```

우리가 구하는 건 계수 `emc` 3개뿐이고, `sbasis`의 모양은 고정입니다.

### sigmoid — 종 모양 만들기

```
def sigmoid(theta):
    result = np.sin(theta) ** 2
    result[(theta > PI) | (theta < 0.0)] = 0.0
    return result
```

`sin²θ`는 `θ = 0 ~ π` 구간에서 0에서 올라가 1을 찍고 다시 0으로 내려오는 **매끄러운 종 모양**입니다. 그 바깥은 0으로 잘라냅니다. 이름은 sigmoid지만 실제 모양은 S자가 아니라 **봉우리**입니다.

### trace_sbasis — 봉우리 3개 배치

```
logTa, logTb, logTc = 5.95, 6.31, 5.0
logT0 = {0: logTa, 1: (logTa+logTb)/2.0, 2: logTb}[index]   # 봉우리 위치
omega = 0.5 * PI / (logTb - logTa)                          # 폭 조절
theta = omega * (logT - logT0) + 0.5 * PI
result = sigmoid(theta)
```

| 기저 | 봉우리 위치 logT | 대응 온도 |
|---|---|---|
| index 0 | 5.95 | 약 0.9 MK — 저온 |
| index 1 | 6.13 | 약 1.3 MK — 중온 |
| index 2 | 6.31 | 약 2.0 MK — 고온 |

`theta = omega·(logT − logT0) + π/2`를 풀어보면, `logT`가 `logT0`에서 **±0.36** 벗어날 때 `theta`가 0과 π에 닿습니다. 즉 각 봉우리의 폭이 `logT` 기준 ±0.36으로 고정됩니다.

### 0번 기저만 특별 처리

```
if index == 0:
    result[logT < logTa] = 1.0    # 봉우리 왼쪽을 1로 채움 (평평하게 연장)
    result[logT < logTc] = 0.0    # 단, logT < 5.0 은 0으로 잘라냄
```

0번은 봉우리가 아니라 **왼쪽이 평평한 계단 모양**이 됩니다. TRACE 채널들이 잘 못 보는 저온 영역(logT 5.0~5.95)을 이 기저가 통째로 담당하도록 만든 설계입니다.

## 3. 코드 흐름

### 입출력 38–39

```
def trace_dem(data, exptimes, K_inv, tau):
    """Return (emc, Tavg, emtot). data (nx,ny,3) DN; K_inv (3,3); tau (3,)."""
```

| 이름 | 크기 | 내용 |
|---|---|---|
| data | (nx, ny, 3) | TRACE 3채널 이미지 (DN). pedestal은 미리 빼둔 상태여야 함 |
| exptimes | (3) | 채널별 노출시간 |
| K_inv | (3, 3) | 외부 입력값. 기저↔관측 변환행렬의 역행렬 |
| tau | (3) | 외부 입력값. 각 기저의 대표온도 |
| emc | (nx, ny, 3) | 출력. 기저별 emission measure 계수 |
| Tavg | (nx, ny) | 출력. EM 가중 평균온도 맵 |
| emtot | (nx, ny) | 출력. 총 emission measure 맵 |

> **`K_inv`와 `tau`는 이 함수가 계산하지 않습니다.** IDL 원본은 내부에서
> `trace_dem_setup`을 호출하지만, 이 과정은 SSW 응답함수와 IDL 적분 루틴에 의존합니다.
> 이번 검증에서는 두 값을 IDL에서 계산해 입력했습니다. 자세한 내용은 5장을 참고하십시오.

### DN → DN/s 47

```
p = data / exptimes[None, None, :]              # DN/s
```

노출시간으로 나눠 **초당 카운트**로 바꿉니다. 노출시간이 다른 이미지들을 같은 기준에 놓기 위한 정규화입니다.

`[None, None,:]`는 (3,)짜리 `exptimes`를 (1,1,3)으로 만들어 브로드캐스팅하는 numpy 표현입니다. IDL은 채널마다 루프를 돌았습니다.

### 핵심 — 픽셀별 역산 50–53

```
for x in range(nx):
    for y in range(ny):
        emc[x, y, :] = K_inv @ p[x, y, :]                        # ①
        Tavg[x, y] = np.sum(emc[x,y,:] * tau) / np.sum(emc[x,y,:])  # ②
```

**①** `emc = K_inv @ p` **— DEM 복원 전체**

`K` 행렬은 "기저 `i`가 채널 `c`에 얼마나 밝게 보이나"를 담은 3×3 표입니다. 순방향은

```
관측(3개) = K × 계수(3개)
```

3×3 정방행렬이므로 **그냥 역행렬을 곱하면 끝**입니다. 반복도, 정규화도, 근사도 없습니다.

simple_reg_dem이 픽셀당 41×41 Cholesky를 수십 번 도는 자리에서, trace_dem은 **3×3 곱셈 한 번**으로 끝냅니다. 미지수를 3개로 고정한 대가이자 보상입니다.

**②** `Tavg` **— EM 가중 평균온도**

```
Tavg = Σ(emc_i × tau_i) / Σ(emc_i)
```

각 기저의 대표온도 `tau`를, 그 기저가 차지하는 양 `emc`로 가중평균한 값입니다. "이 픽셀의 플라즈마는 대략 몇 도인가"를 하나의 숫자로 요약합니다.

> **⚠ 분모가 0이 될 수 있습니다.** `emc` 합이 0에 가까우면 `Tavg`가 발산합니다. 이 코드는 양수성을 강제하지 않으므로 `emc`에 음수가 섞일 수 있고, 그러면 합이 0 근처가 되는 픽셀이 생깁니다. IDL 원본도 동일하게 처리 없이 나눕니다 — 변환본은 이 동작을 그대로 재현했습니다. **사용 시 마스킹 필요.**

### 총 emission measure 54

```
emtot = np.sum(emc, axis=2)                      # IDL total(emc,3)
```

세 basis coefficient를 더해 시선방향의 전체 EM을 계산합니다. IDL의 `total(emc,3)`에서 3은
세 번째 차원을 뜻하며, NumPy는 0부터 세므로 `axis=2`에 해당합니다.

## 4. IDL ↔ Python 대조

| IDL | Python |
|---|---|
| K_inv # reform(p(x,y,*),3) | K_inv @ p[x,y,:] |
| total(emc(x,y,*)*tau) | np.sum(emc[x,y,:]*tau) |
| total(emc, 3)  (3번째 차원) | np.sum(emc, axis=2) |
| where((theta gt pi) or (theta lt 0.0)) | (theta > PI) \| (theta < 0.0)  불린 인덱싱 |
| data(*,*,i)/index(i).sht_mdur  루프 | data / exptimes[None,None,:]  브로드캐스팅 |
| fltarr  (float32) | 내부 float64로 계산 |

### π 상수의 정밀도

```
IDL  (trace_dem_setup.pro 92, 103줄)   pi = 3.14159265359
Python (trace_dem.py 13줄)             PI = 3.14159265359
```

IDL의 내장 `!pi`는 **float32**(3.14159274…)라, 이걸 쓴 코드를 numpy로 옮기면 미세한 차이가 납니다 (실제로 **vdem 변환에서 이 문제가 발생**했습니다).

Kankelborg는 `!pi` 대신 **12자리 상수를 직접 사용했습니다.** Python에서도 같은 상수를 쓰면
이 부분의 반올림 결과를 그대로 재현할 수 있습니다.

## 5. 변환 및 검증 범위

### 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일 5×5 격자 (합성 emc  → forward model → TRACE data) |
| 판정 | 모두 통과 — 13 통과 / 0 실패 |
| 중간값 저장 전후 결과 | 3개 출력 전부 max\|emc − base\| = 0 |
| 회귀 테스트 | 4/4 PASS |

검증 지점은 `01_basis` (기저함수 `sb0/sb1/sb2/T` )와 `99_final` ( `emc/Tavg/emtot` ) 두 곳이고, 둘 다 통과했습니다. 즉 **기저함수 생성과 선형 역산 전체**가 IDL과 일치합니다.

### 옮기지 않은 부분

`K_inv`, `A`, `tau`, `K`는 `trace_dem_setup`이 만드는데, 이 계산이 두 가지에 의존합니다:

`trace_t_resp` — SSW의 TRACE 온도응답 함수

`int_tabulated` — IDL 고유의 5점 Newton-Cotes 적분. 비균일 격자를 spline으로 재샘플한 뒤 적분하는데, **정확한 복제가 어렵습니다**

그래서 이 산출물들은 재구현하지 않고 **IDL에서 계산한 값을 입력**했습니다 (외부 입력을 다시 읽는 방식으로 `02_setup` 검증).

`int_tabulated`와 TRACE 응답함수는 응답행렬을 준비하는 부분이므로 이번 변환에서 제외했습니다.
대신 IDL에서 계산한 `K_inv`와 `tau`를 입력해, 기저함수 생성과 선형 inversion 자체를
검증했습니다. 다른 패키지의 `aia_get_response`, `eit_line_map`, `gaussian_function`도 같은
방식으로 다뤘습니다.

### 사용 시 주의사항

DEM은 세 개의 parameter로 제한되므로 좁은 온도 peak나 다봉 분포를 표현할 수 없습니다.

Positivity를 강제하지 않으므로 `emc`가 음수가 될 수 있고, 이 경우 `Tavg`의 분모가 0에 가까워질 수 있습니다.

IDL 원본 주석에 따라 pedestal을 미리 빼야 합니다. 기준값은 원래 해상도에서 약 83 DN, 2×2 합산 자료에서 약 95 DN입니다.

Basis가 약 `logT=5.0~6.67`만 포함하므로 이 범위 밖의 plasma는 결과에 반영되지 않습니다.

Coefficient 세 개를 채널 세 개에 맞추므로 residual과 χ²가 없습니다. 결과의 물리적 타당성은 별도로 확인해야 합니다.

**기록의 근거.** 코드 설명은 `idl/trace_dem.pro`, `idl/trace_dem_setup.pro`,
`python/trace_dem.py`를 바탕으로 작성했다. 검증 수치는 변환 작업 당시의 비교 보고서를 따른다.

</details>

## 참고 문헌 상태

IDL 원본에는 Kankelborg와 1998년 9월 25일이라는 코드 정보가 있지만, 논문 제목이나 학술지
정보는 적혀 있지 않다. 이번 정리에서는 이 코드의 방법을 직접 설명한다고 확인된 논문을 찾지
못했으므로 `Kankelborg 1998`을 논문처럼 표기하지 않았다.
