# trace_dem 코드 설명

Kankelborg 1998 — TRACE 3채널 선형 DEM

SpaceAI-DEM / `work/trace_dem_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — 다른 DEM과 무엇이 다른가
2. 기저함수 3개 — 이 방법의 전부
3. 코드 한 줄 한 줄
4. IDL ↔ Python 대조
5. 검증 결과와 변환 경계

## 1. 한눈에 — 다른 DEM과 무엇이 다른가

TRACE는 EUV 채널이 **3개뿐**입니다 (171, 195, 284 Å). AIA의 6채널보다도 적습니다.

보통 이 정도면 더 심한 ill-posed 문제가 되어야 하는데, Kankelborg의 선택은 정반대였습니다.

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

우리가 구하는 건 계수 `emc`  3개뿐이고, `sbasis` 의 모양은 고정입니다.

### sigmoid — 종 모양 만들기

```
def sigmoid(theta):
    result = np.sin(theta) ** 2
    result[(theta > PI) | (theta < 0.0)] = 0.0
    return result
```

`sin²θ` 는 `θ = 0 ~ π`  구간에서 0에서 올라가 1을 찍고 다시 0으로 내려오는 **매끄러운 종 모양**입니다. 그 바깥은 0으로 잘라냅니다. 이름은 sigmoid지만 실제 모양은 S자가 아니라 **봉우리**입니다.

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

`theta = omega·(logT − logT0) + π/2` 를 풀어보면, `logT` 가 `logT0` 에서 **±0.36** 벗어날 때 `theta` 가 0과 π에 닿습니다. 즉 각 봉우리의 폭이 `logT`  기준 ±0.36으로 고정됩니다.

### 0번 기저만 특별 처리

```
if index == 0:
    result[logT < logTa] = 1.0    # 봉우리 왼쪽을 1로 채움 (평평하게 연장)
    result[logT < logTc] = 0.0    # 단, logT < 5.0 은 0으로 잘라냄
```

0번은 봉우리가 아니라 **왼쪽이 평평한 계단 모양**이 됩니다. TRACE 채널들이 잘 못 보는 저온 영역(logT 5.0~5.95)을 이 기저가 통째로 담당하도록 만든 설계입니다.

## 3. 코드 한 줄 한 줄

### 입출력 38–39

```
def trace_dem(data, exptimes, K_inv, tau):
    """Return (emc, Tavg, emtot). data (nx,ny,3) DN; K_inv (3,3); tau (3,)."""
```

| 이름 | 크기 | 내용 |
|---|---|---|
| data | (nx, ny, 3) | TRACE 3채널 이미지 (DN). pedestal은 미리 빼둔 상태여야 함 |
| exptimes | (3) | 채널별 노출시간 |
| K_inv | (3, 3) | 주입값. 기저↔관측 변환행렬의 역행렬 |
| tau | (3) | 주입값. 각 기저의 대표온도 |
| emc | (nx, ny, 3) | 출력. 기저별 emission measure 계수 |
| Tavg | (nx, ny) | 출력. EM 가중 평균온도 맵 |
| emtot | (nx, ny) | 출력. 총 emission measure 맵 |

> **⚠ K_inv와 tau는 이 함수가 계산하지 않습니다.** IDL 원본은 내부에서 `trace_dem_setup` 을 호출해 만들지만, 그 계산이 SSW 응답함수와 IDL 고유 적분루틴에 의존해 **변환 경계**로 두고 오라클에서 주입했습니다. 자세한 내용은 5장.

### DN → DN/s 47

```
p = data / exptimes[None, None, :]              # DN/s
```

노출시간으로 나눠 **초당 카운트**로 바꿉니다. 노출시간이 다른 이미지들을 같은 기준에 놓기 위한 정규화입니다.

`[None, None, :]` 는 (3,)짜리 `exptimes` 를 (1,1,3)으로 만들어 브로드캐스팅하는 numpy 표현입니다. IDL은 채널마다 루프를 돌았습니다.

### 핵심 — 픽셀별 역산 50–53

```
for x in range(nx):
    for y in range(ny):
        emc[x, y, :] = K_inv @ p[x, y, :]                        # ①
        Tavg[x, y] = np.sum(emc[x,y,:] * tau) / np.sum(emc[x,y,:])  # ②
```

**①** `emc = K_inv @ p` **— DEM 복원 전체**

`K`  행렬은 "기저 `i` 가 채널 `c` 에 얼마나 밝게 보이나"를 담은 3×3 표입니다. 순방향은

```
관측(3개) = K × 계수(3개)
```

3×3 정방행렬이므로 **그냥 역행렬을 곱하면 끝**입니다. 반복도, 정규화도, 근사도 없습니다.

simple_reg_dem이 픽셀당 41×41 Cholesky를 수십 번 도는 자리에서, trace_dem은 **3×3 곱셈 한 번**으로 끝냅니다. 미지수를 3개로 고정한 대가이자 보상입니다.

**②** `Tavg` **— EM 가중 평균온도**

```
Tavg = Σ(emc_i × tau_i) / Σ(emc_i)
```

각 기저의 대표온도 `tau` 를, 그 기저가 차지하는 양 `emc` 로 가중평균한 값입니다. "이 픽셀의 플라즈마는 대략 몇 도인가"를 하나의 숫자로 요약합니다.

> **⚠ 분모가 0이 될 수 있습니다.** `emc`  합이 0에 가까우면 `Tavg` 가 발산합니다. 이 코드는 양수성을 강제하지 않으므로 `emc` 에 음수가 섞일 수 있고, 그러면 합이 0 근처가 되는 픽셀이 생깁니다. IDL 원본도 동일하게 처리 없이 나눕니다 — 변환본은 이 동작을 그대로 재현했습니다. **사용 시 마스킹 필요.**

### 총 emission measure 54

```
emtot = np.sum(emc, axis=2)                      # IDL total(emc,3)
```

기저 3개의 계수를 그냥 더합니다. 시선방향 전체 플라즈마 양이죠. IDL의 `total(emc,3)` 에서 **3은 "3번째 차원"**이라는 뜻이고, numpy는 0부터 세므로 `axis=2` 가 됩니다.

## 4. IDL ↔ Python 대조

| IDL | Python |
|---|---|
| K_inv # reform(p(x,y,*),3) | K_inv @ p[x,y,:] |
| total(emc(x,y,*)*tau) | np.sum(emc[x,y,:]*tau) |
| total(emc, 3)  (3번째 차원) | np.sum(emc, axis=2) |
| where((theta gt pi) or (theta lt 0.0)) | (theta > PI) \| (theta < 0.0)  불린 인덱싱 |
| data(*,*,i)/index(i).sht_mdur  루프 | data / exptimes[None,None,:]  브로드캐스팅 |
| fltarr  (float32) | 내부 float64로 계산 |

### π 상수 — 운 좋게 피해간 함정

```
IDL  (trace_dem_setup.pro 92, 103줄)   pi = 3.14159265359
Python (trace_dem.py 13줄)             PI = 3.14159265359
```

IDL의 내장 `!pi` 는 **float32**(3.14159274…)라, 이걸 쓴 코드를 numpy로 옮기면 미세한 차이가 납니다 (실제로 **vdem 변환에서 이 문제가 발생**했습니다).

그런데 Kankelborg는 `!pi`  대신 **12자리 상수를 직접 써넣었습니다.** 덕분에 이 코드는 그 함정을 피해갔고, Python도 같은 상수를 그대로 쓰면 정확히 일치합니다.

## 5. 검증 결과와 변환 경계

### 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일 5×5 격자 (합성 emc  → 순방향모델 → TRACE data) |
| 판정 | ALL PASS — 13 pass / 0 fail |
| 계측 무결성 | 3개 출력 전부 max\|emc − base\| = 0 |
| 회귀 테스트 | 4/4 PASS |

검증 지점은 `01_basis` (기저함수 `sb0/sb1/sb2/T` )와 `99_final` ( `emc/Tavg/emtot` ) 두 곳이고, 둘 다 통과했습니다. 즉 **기저함수 생성과 선형 역산 전체**가 IDL과 일치합니다.

### 변환 경계 — 무엇을 옮기지 않았나

`K_inv` , `A` , `tau` , `K` 는 `trace_dem_setup` 이 만드는데, 이 계산이 두 가지에 의존합니다:

`trace_t_resp`  — SSW의 TRACE 온도응답 함수

`int_tabulated`  — IDL 고유의 5점 Newton-Cotes 적분. 비균일 격자를 spline으로 재샘플한 뒤 적분하는데, **정확한 복제가 어렵습니다**

그래서 이 산출물들은 재구현하지 않고 **오라클에서 주입**했습니다 (injection round-trip으로 `02_setup`  검증).

**이 경계가 타당한 이유.** 변환의 목적은 **DEM 알고리즘**을 옮기는 것이지 IDL 적분 루틴을 재현하는 게 아닙니다. `int_tabulated` 는 응답행렬을 만드는 전처리일 뿐이고, 알고리즘의 핵심(기저함수 + 선형 역산)은 전부 변환·검증되었습니다. 이 저장소 전반에서 반복 사용된 **"오라클 입력 주입" 패턴**입니다 — `aia_get_response` , `eit_line_map` , `gaussian_function`  등도 같은 방식으로 처리했습니다.

### 이 solver를 쓸 때 알아야 할 것

**DEM 모양이 3-파라미터로 고정**됩니다. 좁은 온도 피크나 다봉 분포는 표현 불가

**양수성 보장 없음** — `emc` 에 음수가 나올 수 있고, 그러면 `Tavg`  분모가 0에 가까워질 수 있음

**pedestal을 미리 빼야 함** (full res 약 83 DN, 2×2 summed 약 95 DN). IDL 원본 주석의 명시 사항

기저가 `logT`  5.0~6.67 정도만 덮습니다. 그 바깥 온도의 플라즈마는 보이지 않음

잔차가 0이므로 **적합도를 판단할 수단이 없습니다** (χ² 없음). 결과가 물리적으로 말이 되는지는 사용자가 따로 확인해야 함

**출처 표시.** 3·4장의 코드 설명은 `inbox/trace_dem.pro` , `inbox/trace_dem_setup.pro` , `converted/trace_dem.py` 에서, 5장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1·2장의 설명 구성과 다른 solver와의 대비는 코드에서 읽어낸 해석이며, Kankelborg 1998 원논문 대조는 수행하지 않았습니다.
