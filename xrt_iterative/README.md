# xrt_iterative 코드 설명

Weber et al. — Hinode/XRT 반복 DEM (spline knot + MPFIT)

SpaceAI-DEM / `work/xrt_iterative_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — 1391줄 중 35줄만 변환한 이유
2. 전체 방법 — spline knot을 MPFIT으로 맞춘다
3. 변환한 것 — 순방향 모델과 χ²
4. 코드 한 줄 한 줄
5. 변환 경계 — MPFIT을 왜 옮기지 않았나
6. 검증 결과

## 1. 한눈에 — 1391줄 중 35줄만 변환한 이유

`xrt_dem_iterative2.pro` 는 이 저장소에서 **가장 큰 IDL 소스(1391줄)**입니다. 변환된 Python은 **35줄**입니다.

이것은 게으름이 아니라 **의도된 스코핑**이고, 인증서의 판정도 **"PASS (scoped core)"**입니다.

**이 저장소의 원칙:**

경로 의존적인 최적화 루프는 bit-parity가 **원리적으로 불가능**하다.

따라서 **solver가 최적화하는 대상 — 순방향 모델과 χ² — 을 검증하고, 최적화 루프 자체는 문서화된 경계로 남긴다.**

같은 판단이 이 저장소의 세 패키지에 적용됐습니다:

| 패키지 | 경계로 남긴 최적화 | 검증한 핵심 |
|---|---|---|
| xrt_iterative | MPFIT (Levenberg-Marquardt) | 순방향 모델 + χ² |
| vdem | GCV / Brent 최소화 | 코어 선형 역산 (q 고정) |
| firdems | firdem_iterate  (2000회) | 정규화 + first-pass |

## 2. 전체 방법 — spline knot을 MPFIT으로 맞춘다

Weber의 방법은 **파라메트릭 순방향 적합**(aschwanden과 같은 계열)이지만, DEM의 모양 가정이 훨씬 유연합니다.

1. DEM 곡선을 log 공간에서 몇 개의 knot(제어점)을 지나는 spline으로 표현
2. knot 값들을 조금씩 바꿔가며
3. 순방향 모델로 관측 세기를 예측하고
4. χ²를 계산해서
5. MPFIT(Levenberg-Marquardt)이 χ²가 최소가 되는 knot 값을 찾는다

|  | aschwanden aia_teem | xrt_iterative |
|---|---|---|
| DEM 모양 | Gaussian 하나 (파라미터 3개) | spline — knot 개수만큼 자유도 |
| 탐색 방법 | 격자 전수 탐색 | MPFIT (경사 기반 최적화) |
| 경로 의존성 | 없음 → 완전 parity 가능 | 있음 → parity 불가 |

두 코드가 **같은 "파라메트릭 적합" 계열인데 parity 가능성이 갈린 이유**가 여기 있습니다. 전수 탐색은 순서가 고정이라 결정론적이지만, MPFIT은 매 스텝의 방향과 크기가 **직전 계산의 부동소수점 값**에 의존합니다. 1e-16 차이가 스텝 수와 종료점을 바꿉니다.

## 3. 변환한 것 — 순방향 모델과 χ²

위 흐름에서 **3번과 4번**이 변환 대상입니다. MPFIT이 **반복 호출하는 함수**이며, 물리가 들어있는 곳이죠.

```
i_mod = 10^dem ## (emis × 10^t × ln(10^dt))
di    = i_mod − i_obs
chisq = Σ ( di² × weights / i_err² )
```

MPFIT이 이 함수를 수백 번 호출하며 `dem` 을 조정합니다. **이 함수가 IDL과 동일하면, 최적화가 향하는 목표 지형도 동일합니다.** 다만 그 지형을 걸어가는 경로는 구현마다 달라질 수 있습니다.

## 4. 코드 한 줄 한 줄

### 입출력 12–19

```
def xrt_iter_demstat(line_emis, line_t, line_nt, i_obs, i_err, dem, t,
                     weights=None, abunds=None):
```

| 이름 | 크기 | 내용 |
|---|---|---|
| line_emis | (n_line, ·) | 스펙트럼 선별 emissivity. 선마다 고유 온도격자 위에 정의됨 |
| line_t | (n_line, ·) | 각 선의 온도격자 |
| line_nt | (n_line) | 각 선의 유효 격자점 개수 (선마다 다름) |
| i_obs , i_err | (n_line) | 관측 세기와 오차 |
| dem | (nt) | log10(DEM)  — 로그값이 입력 |
| t | (nt) | 공통 온도격자 ( log10 T , 균일 간격) |
| weights | (n_line) | 선별 가중치. 0이면 그 선을 χ²에서 제외 |
| abunds | (n_line) | 원소 존재비 배율 |

### 공통 격자로 재샘플 23–25

```
for i in range(n_line):
    lnt = int(line_nt[i])
    emis[i, :] = np.maximum(np.interp(t, line_t[i, :lnt], line_emis[i, :lnt]), 0.0)
```

선마다 emissivity 테이블의 온도격자가 다르므로, **공통 격자** `t` **로 선형보간**합니다. 그리고 `np.maximum(...,` `0.0)` 으로 **음수를 잘라냅니다** — 보간이 음수를 만들 수 있고 emissivity는 음수일 수 없기 때문입니다.

> **⚠ IDL 함정.** 원본의 `> 0.0` 은 비교가 아니라 **원소별 max**입니다.

### 적분 가중치 만들기 26

```
p[i, :] = emis[i, :] * 10.0 ** t * np.log(10.0 ** dt)
```

이 한 줄이 **적분의 변수 변환**입니다. 유도해보면:

```
구하려는 것:   I = ∫ emis(T) · DEM(T) dT
t = log₁₀T 이므로   T = 10^t,   dT = ln(10) · 10^t · dt
I = Σ  emis · DEM · 10^t · ln(10) · dt
              └──────── p ────────┘
```

코드의 `np.log(10.0 ** dt)` 는 `ln(10^dt) = dt · ln(10)` 이므로, `p = emis × 10^t × dt × ln(10)` 가 되어 위 식과 정확히 일치합니다.

`dt · ln(10)` 을 직접 쓰지 않고 `ln(10^dt)` 로 우회한 것은 IDL 원본의 표기를 그대로 옮긴 것입니다. 수학적으로 동일하지만 **부동소수점 연산 순서가 달라지면 결과가 미세하게 달라질 수 있어**, parity를 위해 원본 형태를 유지했습니다.

### 순방향 모델과 χ² 27–32

```
ldem  = 10.0 ** dem                 # log10(DEM) → DEM
i_mod = (p @ ldem) * abunds         # 예측 세기 (n_line)
di    = i_mod - i_obs               # 잔차
chisq_arr = di ** 2 * weights / i_err ** 2
chisq_arr[weights == 0.0] = 0.0     # ← 가중치 0인 선을 명시적으로 배제
chisq = np.sum(chisq_arr)
```

> **⚠ 31줄이 왜 필요한가.** `weights = 0` 이면 `di² × 0` 이라 이미 0일 것 같지만, `i_err` 도 0인 선이 있으면 `0 × 0 / 0 = NaN` 이 됩니다. 그 NaN 하나가 `np.sum`  전체를 NaN으로 오염시키죠. 이 줄은 **비활성 선을** **확실히 배제하는 안전장치**이고, IDL 원본의 동작을 그대로 옮긴 것입니다.

### float32 유지 21–22

```
emis = np.zeros((n_line, nt), dtype=np.float32)     # IDL FLTARR
p    = np.zeros((n_line, nt), dtype=np.float32)
```

IDL이 `FLTARR` (float32)로 잡으므로 Python도 맞췄습니다. 이 저장소에서 반복 확인된 원칙입니다 — **"IDL의** **dtype을 그대로 따라가야 parity가 성립한다."**

## 5. 변환 경계 — MPFIT을 왜 옮기지 않았나

### ❌ 제외: MPFIT / Levenberg-Marquardt 반복

`xrt_dem_iter_solver`  → MPFIT으로 spline knot을 조정하는 루프입니다.

**bit-parity가 불가능한 이유.** MPFIT은 매 스텝에서 수치 미분으로 Jacobian을 만들고, damping 파라미터를 조정하며 방향과 보폭을 정합니다. 이 모든 결정이 **직전 계산 결과에 의존**하므로, 1e-16 수준의 부동소수점 차이가 **스텝 개수·종료 시점·최종 knot 값**을 바꿉니다. IDL MPFIT과 scipy의 구현이 다르므로, 같은 지형에서도 **다른 경로로 다른 지점에 착지**합니다.

### ❌ 제외: CHIANTI/XRT 응답 구성 ( make_xrt_temp_resp )

스캐폴딩입니다. DEM 알고리즘이 아니라 입력을 만드는 전처리이므로, 오라클에서 주입했습니다.

### ❌ 제외: Monte Carlo 오차막대

원본에는 `RANDOMU` 로 데이터를 교란해 오차막대를 만드는 기능이 있으나, **이 코어에는 난수가 없습니다.** 확률적이라 parity 대상이 아니며 비활성 상태로 검증했습니다.

## 6. 검증 결과

| 항목 | 값 |
|---|---|
| 세트 | set1 — 합성 5개 선, 26점 logT  격자 위 Gaussian emissivity |
| 판정 | ALL PASS — 12 pass / 0 waived / 0 fail ( compare_probes  exit 0) |
| 계측 무결성 | max\|i_mod − base\| = 0 , max\|chisq − base\| = 0 |
| 정확도 | i_mod / chisq  < 1e-6 일치. chisq = 18.2523041  정확 일치 |
| 회귀 테스트 | 3/3 |

**검증 설계의 요점:** `i_obs` 를 **다른 DEM**에서 순방향 모델로 만들었습니다. 그래야 `di ≠ 0` 이 되어 χ² 계산 경로가 실제로 검증됩니다. 같은 DEM을 썼다면 잔차가 0이라 대부분의 산술이 무의미해졌을 겁니다.

### 검증의 한계 — 정직하게

`interpol` **이 실질적으로 검증되지 않았습니다.** 테스트에서 각 선의 고유 격자가 공통 격자와 동일했기 때문에, 보간이 **항등 연산(identity)**으로만 실행됐습니다. 격자가 실제로 다른 경우의 보간 정확도는 미검증입니다.

인증서는 "범위 내에서 linear `interpol`  == `np.interp` "라고 기록하고 있으나, **범위 밖 외삽 동작**은 두 구현이 다를 수 있습니다 (IDL `interpol` 은 선형 외삽, `np.interp` 는 **끝값으로 고정**). 실제 CHIANTI 데이터로 쓸 때 확인이 필요한 지점입니다.

### 알고 써야 할 것

이 함수는 **DEM을 구해주지 않습니다.** 주어진 DEM에 대해 세기와 χ²를 계산할 뿐입니다. 실제 역산에는 최적화 루프가 별도로 필요합니다

`dem`  입력은 `log10(DEM)` 입니다. DEM 자체가 아닙니다

`t` 는 **균일 간격**이어야 합니다 — `dt = t[1] − t[0]` 를 한 번만 계산해 전 구간에 씁니다

`weights = 0` 인 선은 χ²에서 빠지지만 `i_mod` 는 계산됩니다

**출처 표시.** 4장의 코드 설명은 `converted/xrt_demstat.py` 와 그 docstring에서, 5·6장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1~3장의 방법 설명과 다른 패키지와의 대비는 `work/DEM_SOLVER_METHODS.md` · `work/CONVERSION_SUMMARY.md` 와 코드에 근거한 해석이며, Weber et al. 원논문 대조는 수행하지 않았습니다. `np.interp` 와 IDL `interpol` 의 외삽 동작 차이는 두 함수의 일반적 사양에 근거한 지적으로, 이 데이터셋에서 실측 확인한 것은 아닙니다.
