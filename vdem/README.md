# vdem 코드 설명

Newton, Emslie & Mariska 1995 — Velocity DEM (속도 DEM), GCV 정규화

SpaceAI-DEM / `work/vdem_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 먼저 — 이것은 온도 DEM이 아닙니다
2. 속도 격자 만들기 — Doppler 변환
3. 커널 — 속도에서 파장으로
4. 3차 차분 정규화
5. 선형 역산과 오차막대
6. `!pi` 가 float32였다 — 이 패키지의 핵심 발견
7. 검증 결과와 경계

## 1. 먼저 — 이것은 온도 DEM이 아닙니다

> **⚠ 물리량이 다릅니다.** 이 저장소의 나머지 10개는 모두 **온도** 분포(DEM(T))를 구합니다. vdem이 구하는 것은 **시선방향 속도 분포**입니다 — 단위는 `photons/s/(km/s)` .

같은 "DEM"이라는 이름을 쓰지만 **결과를 서로 비교하거나 섞어 쓸 수 없습니다.**

### 무엇을 푸는가

관측은 **스펙트럼 선 하나의 프로파일**입니다. 파장별 밝기 곡선이죠.

```
관측된 선 프로파일 = ∫ (속도 v인 플라즈마의 양) × (그 플라즈마가 만드는 스펙트럼 모양) dv
```

플라즈마가 시선방향으로 움직이면 Doppler 효과로 선이 이동하고, 여러 속도의 플라즈마가 섞여 있으면 선이 넓어지거나 비대칭이 됩니다. **그 프로파일을 풀어 속도별 분포를 복원**하는 것이 목표입니다.

수학적 구조는 온도 DEM과 **완전히 같습니다** — 커널로 적분된 분포를 역산하는 ill-posed 문제. 다만 적분 변수가 온도가 아니라 **속도**이고, 커널이 응답함수가 아니라 **Gaussian 선 프로파일**일 뿐입니다.

|  | 온도 DEM (예: demreg) | vdem |
|---|---|---|
| 관측 | 채널별 세기 (6개) | 파장별 세기 (스펙트럼) |
| 미지수 | 온도별 DEM | 속도별 분포 |
| 커널 | 기기 온도응답 | Gaussian 선 프로파일 |
| 정규화 | 0차 / 1차 도함수 | 3차 차분 |
| λ 선택 | discrepancy principle | GCV (이번 검증에서는 우회) |

## 2. 속도 격자 만들기 — Doppler 변환 61–69

```
diff = wave[1] - wave[0]
dvel = C_LIGHT * diff / L0                            # 파장 간격 → 속도 간격
vmax = -1.0 * C_LIGHT * (np.min(wave) - L0) / L0
vmin = -1.0 * C_LIGHT * (np.max(wave) - L0) / L0
nvel = int(np.trunc((vmax - vmin) / dvel))
velocity[0] = vmin + dvel / 2.0                       # 첫 칸의 중심
for i in range(1, nvel):
    velocity[i] = velocity[i-1] + dvel
```

Doppler 공식 `v = −c·(λ − λ₀)/λ₀` 를 그대로 적용합니다. `L0` 는 정지 파장입니다.

**부호가 뒤집힙니다** — 긴 파장(적색편이)이 **음의** 속도(멀어짐)에 대응하도록 `−1` 이 붙습니다. 그래서 `vmax` 가 `min(wave)` 에서, `vmin` 이 `max(wave)` 에서 나옵니다

`velocity[i]` 를 **누적 덧셈**으로 만듭니다. `vmin + i*dvel` 로 계산하는 것과 부동소수점 결과가 미세하게 다르므로, parity를 위해 원본 방식을 유지했습니다

> **⚠** `fix` **→** `np.trunc` **.** IDL `fix` 는 **0 방향 절삭**입니다. Python `int()` 도 같지만, `np.floor` 로 옮기면 음수에서 달라집니다.

## 3. 커널 — 속도에서 파장으로 72–80

```
v_instr = width_sumer * C_LIGHT / L0                  # 기기 폭 → 속도 단위
center  = L0 * (1.0 - velocity / C_LIGHT)             # 속도 i의 선 중심 파장
sig     = L0 * np.sqrt(K_BOLTZ*tbar/m + xi**2 + v_instr**2) / C_LIGHT
for j in range(nflux):        # 파장
    for i in range(nvel):     # 속도
        val = dvel / (sqrt(2π) * sig) * exp( -(wave[j] - center[i])**2 / (2*sig**2) )
        basekernel[i, j] = 0.0 if val < 1.0e-15 else val
```

### 선폭 sig 의 세 성분

```
sig ∝ √( k·T̄/m  +  ξ²  +  v_instr² )
         └ 열운동 ┘  └난류┘  └ 기기 ┘
```

| 성분 | 내용 |
|---|---|
| k·T̄/m | 열적 Doppler 폭. 이온 질량 m 이 클수록 좁아짐 |
| ξ² | 비열적 난류 속도 (입력 파라미터) |
| v_instr² | 기기 분해능 ( width_sumer ) |

세 폭을 제곱합으로 더합니다 — 독립적인 Gaussian들이 합성될 때의 표준 규칙입니다.

`basekernel[i,j]`  = "속도 `i` 인 플라즈마가 파장 `j` 에 만드는 밝기".

**온도 DEM의 응답행렬** `Rij` **와 정확히 같은 역할**입니다.

`val < 1e-15 → 0`  절단은 커널 꼬리를 잘라 **수치적으로 무의미한 값이 행렬에 남지 않게** 합니다.

## 4. 3차 차분 정규화 85–91

```
ke3 = np.zeros((nvel, nvel-3))
for n in range(nvel-3):
    ke3[n,   n] = -1.0
    ke3[n+1, n] =  3.0
    ke3[n+2, n] = -3.0
    ke3[n+3, n] =  1.0
baseh = ke3 @ ke3.T
h = baseh * 1e16 * (1/6) / 1e4
```

`(−1, 3, −3, 1)` 은 **3차 유한차분**의 계수입니다. 즉 `ke3ᵀ·x` 는 해의 3차 도함수를 근사하고, `xᵀ·h·x` 는 **3차 도함수의 제곱합**이 됩니다.

| 정규화 차수 | 선호하는 해 | 사용 패키지 |
|---|---|---|
| 0차 | 값이 작은 해 | demreg |
| 1차 | 기울기가 완만한 해 | simple_reg_dem |
| 3차 | 국소적으로 2차식에 가까운 해 | vdem |

| 3차 차분은 매우 관대한 제약입니다. 3차 도함수가 0인 함수 = 2차 이하 다항식이므로, 직선·포물선은 전 혀 벌점을 받지 않습니다. 속도 분포가 한두 개의 매끄러운 봉우리일 것이라는 사전지식과 잘 맞습니다. 그 대신 행렬의 밴드 폭이 넓어집니다 — simple_reg_dem의 regmat 이 3줄(tridiagonal)이었다면, 여기 h = ke3·ke3ᵀ는 7줄짜리입니다. |  |  |
|---|---|---|
|  | h |  |

## 5. 선형 역산과 오차막대 97–111

```
for time in range(specnum):                      # 스펙트럼마다
    kernel = basekernel / eflux[:, time][None, :]   # ① 오차 가중
    A      = kernel @ kernel.T                       # ② KᵀK
    inten  = flux[:, time] / eflux[:, time]          # ③ 가중된 데이터
    data   = kernel @ inten                          # ④ Kᵀg
    initial = np.trace(A) / TrH                      # ⑤ 스케일 맞추기
    q       = weight * initial                       # ⑥ 정규화 세기
    matrix0 = A + q * h                              # ⑦ 정규화된 정규방정식
    x       = np.linalg.solve(matrix0, data)         # ⑧ 풀기
    inverse    = np.linalg.inv(matrix0)
    covariance = inverse @ A @ inverse               # ⑨ 공분산
    evdem[:, time] = SCALE * np.sqrt(np.diag(covariance))   # ⑩ 오차막대
    retflux[:, time] = (x @ basekernel) / SCALE            # ⑪ 순방향 재계산
```

### ①~⑧ — 익숙한 구조

`matrix0 = A + q·h` 는 simple_reg_dem의 `amat = JᵀJ + regmat` 과 **정확히 같은 형태**입니다. 차이는 이 문제가 **선형**이라 반복이 필요 없다는 점입니다 — 한 번 풀면 끝입니다.

**⑤** `initial = Tr(A) / Tr(H)` **— 정규화 세기의 기준**

두 행렬의 **대각합 비율**입니다. 이렇게 잡으면 `q·h` 가 `A` 와 **비슷한 크기 규모**가 되어, 데이터 항과 정규화 항이 균형을 이룹니다. `weight = 1.0` 이 기본값이므로 `q = initial` 이 됩니다.

원래 이 `q` 는 **GCV로 자동 선택**되는 값입니다. 이번 검증에서는 **WEIGHT 모드**로 고정해 결정론적으로 만들었습니다 (7장 참조).

### ⑨~⑩ — 오차막대

```
covariance = (A + qh)⁻¹ · A · (A + qh)⁻¹
evdem      = scale × √( diag(covariance) )
```

정규화된 역산의 **선형 오차 전파**입니다. 데이터 오차가 해로 어떻게 퍼지는지 계산합니다.

**오차 추정이 있는 solver입니다.** simple_reg_dem·trace_dem·aia_teem에는 없던 출력이고, demreg의

`edem` ·dem_sites의 `demerr` 와 같은 계열입니다. 다만 **정규화 편향(bias)은 포함하지 않은 분산만**의 추정입니다.

### ⑪ 순방향 재계산

`retflux = x @ basekernel` 로 복원된 속도분포에서 스펙트럼을 다시 만듭니다. 원래 관측과 비교하면 적합도를 눈으로 확인할 수 있습니다.

## 6. !pi 가 float32였다 — 이 패키지의 핵심 발견

이 변환에서 가장 값진 발견이고, 다른 IDL 코드를 옮길 때도 반드시 알아야 할 내용입니다.

> **⚠ IDL의** `!pi` **는 float32입니다** — 값이 `3.14159274` 입니다.

double이 필요하면 `!dpi` 를 써야 합니다.

> **⚠ IDL의** `e` **리터럴도 전부 float32입니다.** `1.38054e-16` , `2.997925e10` , `1.495979e13` , `1.6713e-24` , `0.99` , `1e16` , `1e4`  — 모두 float32로 저장됩니다. double을 원하면 `1.38054d-16` 처럼 `d` 를 써야 합니다.

Python은 이 모든 상수가 **기본적으로 float64**입니다. 그냥 옮기면 **모든 물리상수가 미세하게 다른 값**이 되고, 결과에 **약 3e-8의 계통적 차이**가 남습니다.

### 어떻게 해결했나

```
def _f32(x):
    """IDL `e` literals are FLOAT32; reproduce their exact value in float64 arithmetic."""
    return float(np.float32(x))
K_BOLTZ = _f32(1.38054e-16)      # float32로 반올림한 값을 float64에 담음
C_LIGHT = _f32(2.997925e10)
_PI32   = float(np.float32(np.pi))
```

**요령:** 값을 `np.float32()` 로 **한 번 통과시켜 IDL과 같은 비트로 만든 뒤**, `float()` 로 float64에 담아 이후 산술은 double로 합니다.

IDL도 **상수는 float32지만 계산은 double로 승격**하는 경우가 많기 때문에, "상수만 float32, 연산은 float64"가 정확한 재현입니다.

### 더 미묘한 부분 — 연산 순서까지

```
_SSD32   = np.float32(0.99) * _AU32          # float32 × float32 → float32
_SSD2    = float(_SSD32 * _SSD32)            # float32 제곱 → float32
_SCALE_F32 = float((np.float32(4.0) * np.float32(np.pi)) * np.float32(_SSD2))
```

IDL에서 `sun_sat_dist = 0.99*au` 와 그 제곱은 **float32 연산으로 수행**됩니다. 반면 `flux`  스케일링은 double 변수에 누적되므로 double입니다. 같은 `4·π·ssd²` 라도 **쓰이는 위치에 따라 정밀도가 다릅니다.**

주석 37–39줄이 이 구분을 설명하고 있습니다 — `evdem` / `retflux` 에서는 **standalone float32 곱**이고, `flux` 스케일링에서는 double 누적입니다.

### 다른 패키지에서도 반복된 교훈

| 패키지 | float32 이슈 |
|---|---|
| vdem | !pi 와 e  리터럴 |
| aschwanden | fit  전체를 float32로 — float64면 χ가 1e-3 어긋남 |
| dem_sites | gaussian_function  커널 |
| chianti_dem | IDL SPLINE 전체 |
| trace_dem | 원저자가 !pi  대신 상수를 직접 써서 우연히 회피 |

## 7. 검증 결과와 경계

### 결과

| 항목 | 값 |
|---|---|
| 세트 | set1 — 합성 시프트/광폭 방출선 프로파일, 단일 스펙트럼, nvel = 50 |
| 판정 | ALL PASS — 23 pass / 0 fail |
| 계측 무결성 | max\|vdem − base\| = 0 |
| 검증 지점 | 00_input , 01_kernel , 02_h , 03_setup , 99_final  — 전부 strict 통과 |
| 회귀 테스트 | 4/4 |

float32 수정 **전**에는 약 `3e-8` 의 계통차가 있었고, 상수를 IDL과 맞춘 **후** strict double parity를 달성했습니다.

### 경계 — 검증하지 않은 것

| 항목 | 사유 |
|---|---|
| GCV / Brent 최적화 ( gcvcompute  + gcvfunction ) | WEIGHT 모드로 우회. Brent 반복은 경로 의존적이라 bit-parity가 어렵습니다 ( firdem_iterate 와 같은 성격). 코어 선형 역산은 완전 검증됨 |
| CAXIX 분기 | d13 위성선, excitcoef.dat  필요 |
| VERBOSE 플롯 | 시각화 |
| prep_vdem_bcs  / prep_vdem_sumer | 데이터 준비 스캐폴딩 |
| 실관측 다스펙트럼 | 단일 합성 스펙트럼으로만 검증. SUMER/BCS 실데이터는 후속 |

**대화형 프롬프트 우회.** 원본 `vdem.pro` 는 `read` 로 사용자 입력을 받는 부분이 있어 헤드리스 실행이 불가능했습니다. `ION_IN`  키워드로 우회해 배치 실행했습니다.

### 알고 써야 할 것

**결과는 속도 분포입니다.** 온도 DEM과 혼동하면 안 됩니다

`weight` 를 지정하지 않고 GCV를 쓰면 **이 검증 범위를 벗어납니다**

`tbar` (가정 온도), `xi` (난류), `width_sumer` (기기폭)가 커널 폭을 결정합니다 — **이 값이 틀리면 결과 전체가 틀립니다**

속도 격자 개수 `nvel` 은 관측 파장 범위와 간격에서 **자동 결정**됩니다. 직접 지정할 수 없습니다

`eflux` (오차)가 0인 파장이 있으면 **0으로 나누기**가 발생합니다

오차막대 `evdem` 은 **분산만**의 추정이며 정규화로 인한 편향은 포함하지 않습니다

**출처 표시.** 2~6장의 코드 설명은 `inbox/vdem.pro` 와 `converted/vdem.py` (주석 포함)에서, 7장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1장의 온도 DEM과의 구조 비교와 4장의 정규화 차수 비교는 `work/DEM_SOLVER_METHODS.md` 와 코드에 근거한 해석이며, Newton, Emslie & Mariska 1995 원논문 대조는 수행하지 않았습니다.
