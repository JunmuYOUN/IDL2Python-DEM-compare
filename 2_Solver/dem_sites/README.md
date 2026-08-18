# dem_sites — 오차를 되돌려 보정하는 반복 DEM

SITES는 큰 행렬의 역행렬을 직접 구하지 않는다. 임시 DEM으로 관측값을 다시 계산한 다음,
실제 관측과의 차이를 온도 구간에 나누어 더하는 작업을 반복한다.

## 먼저 읽기

한 번의 반복은 다음 순서로 진행된다.

1. 현재 DEM으로 각 채널의 밝기를 계산한다.
2. 실제 밝기와 계산한 밝기의 차이를 구한다.
3. 그 차이를 응답함수에 따라 온도 구간에 나누어 되돌린다.
4. 결과를 부드럽게 만들고 음수는 0으로 바꾼다.
5. 차이가 충분히 작아지면 멈춘다.

관측 오차뿐 아니라 응답함수의 오차도 입력할 수 있다. 핵심 함수 `dem_sites`는 한 번에 한
픽셀을 계산한다. 여러 픽셀을 빠르게 처리하는 Grid-SITES는 IDL 원본에는 있지만 이번 Python
변환에는 포함하지 않았다.

합성 Gaussian DEM 한 픽셀에서 19개 비교 항목과 회귀 테스트 4개를 통과했다. 관측값이 0인
채널은 나눗셈 오류를 만들 수 있으므로 먼저 확인해야 한다. 여러 봉우리나 급격한 변화가 있는
DEM과 실제 관측 대형 맵은 이번에 확인하지 않았다.

- IDL 원본과 원문 안내: `idl/`
- Python 변환본: `python/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 방법 요약 — iterative inversion
2. 준비 단계 — 가중치와 배분표
3. 반복 계산의 핵심
4. 오차 추정과 마무리
5. IDL ↔ Python 대조
6. 검증 결과와 한계

## 1. 방법 요약 — iterative inversion

SITES는 `simple_reg_dem`의 Cholesky 분해나 `trace_dem`의 3×3 inverse처럼 선형계를 한 번에
풀지 않습니다.

**SITES는 행렬을 아예 뒤집지 않습니다.**

대신 현재 DEM으로 관측값을 예측하고, residual을 온도 공간으로 back-project해 DEM을 보정하는
과정을 수렴할 때까지 반복합니다. Landweber/Richardson–Lucy 계열과 비슷한 iterative
deconvolution 방식입니다.

구조상 이런 성질이 따라옵니다:

| 항목 | SITES의 처리 |
|---|---|
| 역행렬 | 없음. 행렬 곱셈과 덧셈만 사용 → 조건수 문제에서 자유로움 |
| 정규화 | 별도 penalty 항 대신 매 반복의 Gaussian smoothing과 early stopping 사용 |
| 양수성 | 매 반복 음수를 0으로 잘라냄 ( >0 ) |
| 오차 추정 | 있음. demerr 를 함께 반환 (simple_reg_dem에는 없던 것) |
| 적용 단위 | 단일 픽셀. 맵 처리는 상위 래퍼(Grid-SITES)의 역할 |

## 2. 준비 단계 — 가중치와 배분표

### 입출력 28–30

```
def dem_sites(obs_in, err_in, response, response_err, delta_temp,
              convergence=1e-2, ker=None):
    """Return (demmain, demerr, obsmod, irep). response is [nt,nwl]."""
```

| 이름 | 크기 | 내용 |
|---|---|---|
| obs_in | (nwl) | 이 픽셀의 채널별 관측값 |
| err_in | (nwl) | 관측 오차 |
| response | (nt, nwl) | 온도응답 행렬 |
| response_err | (nwl) | 응답 자체의 상대 불확실도 — 기기 보정의 신뢰도 |
| delta_temp | (nt) | 온도 구간 폭 |
| convergence | 스칼라 | 종료 기준. 기본 1e-2  = 잔차 1% |
| demmain | (nt) | 출력. 복원된 DEM |
| demerr | (nt) | 출력. DEM 오차 |
| obsmod | (nwl) | 출력. 최종 DEM으로 다시 계산한 관측 예측값 |
| irep | 정수 | 출력. 몇 번 반복했나 |

### 채널 가중치 40

```
wt = 1.0 / np.sqrt((err_in / obs_in) ** 2 + response_err ** 2)
```

두 종류의 불확실도를 합쳐 채널별 신뢰도를 만듭니다:

`err_in/obs_in` — **상대 측정오차**. 신호 대비 잡음이 얼마나 큰가

`response_err` — **응답함수 불확실도**. 이 채널의 기기 보정을 얼마나 믿을 수 있나

둘을 제곱합해 루트를 씌운 뒤 역수를 취하므로, **잡음이 적고 응답이 잘 알려진 채널일수록 큰 가중치**를 받습니다.

`response_err`를 별도 입력으로 받는 것이 이 solver의 특징입니다. 대부분의 DEM 코드는 측정오차만 보지만, SITES는 **"응답함수 자체도 틀릴 수 있다"**를 명시적으로 모델에 넣습니다.

### gres — 배분표 41–42

```
num  = response * wt[None, :]
gres = num / np.sum(num, axis=1)[:, None]
```

`np.sum(num, axis=1)`은 **채널 방향** 합이므로, `gres`는 각 온도 행이 **합쳐서 1**이 되도록 정규화됩니다.

```
gres[t, w] = "온도 t의 신호 중 채널 w가 담당하는 비율" (신뢰도 가중)
           채널1   채널2   채널3  ...    합
온도 t  [  0.42    0.31    0.27  ... ]  = 1.0
```

나중에 반복 루프에서 **잔차를 온도공간으로 되돌릴 때 이 비율대로 배분**합니다. 알고리즘의 핵심 장치입니다.

### 스무딩 커널 45–48

```
nker = max(0.08 * nt, 0.5)
if ker is None:
    ker = gaussian_function(nker)
```

커널 폭을 **온도격자 개수의 8%**로 잡습니다 (최소 0.5). 격자가 촘촘하면 커널도 비례해 넓어지므로, 물리적 스무딩 폭이 격자 해상도와 무관하게 유지됩니다.

```
def gaussian_function(sigma):
    n = int(2 * np.ceil(3 * sigma) + 1)      # ±3σ 범위, 홀수 길이
    x = np.arange(n) - (n - 1) / 2.0         # 중앙이 0
    g = np.exp(-x ** 2 / (2.0 * sigma ** 2))
    return (g / np.sum(g)).astype(np.float32)   # 합=1로 정규화
```

> **⚠** `float32` **가 의도적입니다.** IDL의 `gaussian_function(σ, /norm)`이 float32를 반환하므로 Python도 맞췄습니다. 이 커널이 코드 전체에서 **유일한 float32 요소**이고, 여기가 어긋나면 반복이 진행되며 차이가 누적됩니다.

### 미리 계산해두는 행렬 52–54

```
res2    = response * delta_temp[:, None]     # [nt,nwl] 응답 × 온도폭
totres2 = np.sum(res2 ** 2, axis=0)          # [nwl]    채널별 제곱합
res3    = res2 * gres                        # [nt,nwl] 배분표 반영
```

|  | 뜻 |
|---|---|
| res2 | res2[t,w] × DEM[t]  = 온도구간 t가 채널 w에 기여하는 양. 순방향 모델의 핵심 |
| totres2 | 채널별 정규화 상수. 반복에서 단계 크기를 정함 |
| res3 | forward response와 배분표의 곱. residual back-projection에 사용 |

## 3. 반복 계산의 핵심

```
while True:
    dem     = (obs / totres2)[None, :] * res3                    # ① back-projection
    conv    = np.convolve(np.sum(dem, axis=1), ker, mode="same") # ② smoothing
    demmain = np.maximum(demmain + conv, 0.0)                    # ③ 누적 + 양수화
    obs     = obs_in - np.sum(demmain[:, None] * res2, axis=0)   # ④ 새 잔차
    convobs = np.sum(np.abs(obs * obswt)) / totwt                # ⑤ 수렴 측정
    irep   += 1
    if irep > 300 or convobs < convergence:
        break
```

루프 진입 직전 `obs = obs_in` (전체 관측), `demmain = 0` (빈 DEM)에서 시작합니다.

### ① Back-projection 63

```
dem = (obs / totres2)[None, :] * res3
```

`obs`는 "아직 설명하지 못한 관측량"입니다. 이를 `totres2`로 나눠 단계 크기를 정하고, `res3`를 곱해 **온도공간에 펼칩니다.**

결과 `dem`은 `[nt, nwl]` — "채널 w의 잔차가 온도 t에 얼마를 기여해야 하는가"입니다. 아직 채널별로 나뉜 상태입니다.

### ② Smoothing 64

```
conv = np.convolve(np.sum(dem, axis=1), ker, mode="same")
```

`np.sum(dem, axis=1)`로 채널 방향을 합쳐 온도축 벡터 `[nt]`를 만든 뒤, **Gaussian 커널로 흐립니다.**

**여기가 SITES의 정규화입니다.** simple_reg_dem이 `regmat`을 방정식에 더해 진동을 억제한다면, SITES는 **매 반복마다 증분을 물리적으로 뭉개서** 같은 효과를 냅니다. 페널티 함수도, λ도, 행렬도 없습니다.

### ③ 누적하고 음수 잘라내기 65

```
demmain = np.maximum(demmain + conv, 0.0)
```

흐린 증분을 기존 DEM에 **더합니다.** 그리고 **음수가 된 부분은 0으로 잘라냅니다.**

> **IDL 문법 주의.** 원본의 `...>0`은 비교가 아니라 0과 비교해 큰 값을 고르는 연산이며,
> NumPy의 `np.maximum`에 해당합니다. 이 연산이 매 반복에서 positivity constraint를 적용합니다.

### ④ 새 잔차 계산 66

```
obs = obs_in - np.sum(demmain[:, None] * res2, axis=0)
```

갱신된 DEM으로 **순방향 모델**을 돌려 관측을 예측하고, 실제 관측에서 뺍니다. 남은 것이 다음 반복이 설명할 몫입니다.

`axis=0`은 온도 방향 합이므로 결과는 `[nwl]` — 채널별 잔차입니다.

### ⑤ 수렴 판정 67, 72

```
obswt   = wt / obs_in                                # 루프 밖에서 미리 계산 (60줄)
totwt   = np.sum(wt)
convobs = np.sum(np.abs(obs * obswt)) / totwt
if irep > 300 or convobs < convergence:
    break
```

`obs × obswt = obs/obs_in × wt` 이므로, `convobs`는 **가중 평균 상대잔차**입니다. 기본값 `1e-2`는 **"평균적으로** **관측값의 1% 이내로 맞으면 그만"**이라는 뜻입니다.

**이 종료 조건이 두 번째 정규화입니다.** 잔차를 0까지 밀어붙이면 잡음까지 맞추게 되므로, 1%에서 멈추는 것 자체가 overfit 방지 장치입니다. 개념적으로 **early stopping**이고, simple_reg_dem의 χ²=1 targeting과 목적이 같습니다.

안전장치로 **300회 상한**도 걸려 있습니다. 실제 검증에서는 **8회** 만에 수렴했습니다.

## 4. 오차 추정과 마무리

### DEM 오차 57–58, 76

```
e      = (err_in / obs_in) ** 2 + response_err ** 2     # 채널별 상대분산
demerr = np.sqrt(np.sum(e[None, :] * gres, axis=1))     # 배분표로 온도공간에 전파
...
demerr = demerr * demmain                                # 상대 → 절대
```

채널별 상대 불확실도를 `gres` (온도별 배분 비율)로 가중합해 **온도별 상대오차**를 만들고, 마지막에 DEM을 곱해 **절대오차**로 바꿉니다.

**simple_reg_dem에는 없던 출력입니다.** 오차막대가 필요한 연구에서는 이 점이 SITES의 실질적 장점입니다. 다만 이것은 **선형 전파 근사**이지, 반복 과정의 비선형성(음수 잘라내기, 조기 종료)까지 반영한 엄밀한 불확실도는 아닙니다.

### 순방향 재계산 75

```
obsmod = np.sum((demmain * delta_temp)[:, None] * response, axis=0)
```

최종 DEM으로 관측을 다시 예측합니다. `obs_in`과 비교하면 **적합도를 사용자가 직접 확인**할 수 있습니다. trace_dem에는 없던 진단 수단입니다.

## 5. IDL ↔ Python 대조

이 코드는 IDL의 `rebin` 브로드캐스팅 관용구가 많아 변환 시 **축 방향을 헷갈리기 쉽습니다.**

| IDL | Python |
|---|---|
| rebin(reform(v,1,nwl),nt,nwl)  — 행 복제 | v[None,:] |
| rebin(v_nt,nt,nwl)  — 열 복제 | v[:, None] |
| total(x, 2)  — 2번째 차원(nwl) 합 | np.sum(x, axis=1) |
| total(x, 1)  — 1번째 차원(nt) 합 | np.sum(x, axis=0) |
| convol(a, ker, /edge_zero) | np.convolve(a, ker, 'same') |
|...>0  (원소별 max) | np.maximum(..., 0.0) |
| repeat... endrep until cond | while True:... if cond: break |
| temporary(demmain)  (메모리 최적화) | 불필요 — 대응물 없음 |
| arg_present(demerr)  (호출자가 요청했나) | 불필요 — 항상 계산해 반환 |
| gaussian_function(σ, /norm) | 재구현 ( 2*ceil(3σ)+1, float32) |

> **⚠** `convol` **의 경계 처리.** IDL은 `edge_zero=1, edge_trunc=0`로 **가장자리를 0으로 패딩**합니다. `np.convolve(..., 'same')`이 같은 동작이라 일치했지만, IDL `convol`은 기본값이 다르므로 다른 코드를 옮길 때는 매번 키워드를 확인해야 합니다.

## 6. 검증 결과와 한계

### 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일 픽셀 set1 — 합성 Gaussian DEM 순방향 → 관측값 |
| 판정 | 모두 통과 — 19 통과 / 0 실패 |
| 중간값 저장 전후 결과 | max\|dem − base\| = 0 |
| 일치 확인 | 반복 횟수 irep = 8  일치, DEM 피크 9.567e21  일치 |
| 회귀 테스트 | 4/4 (IDL 기준값 고정 + gaussian_function  float32 일치 + max 연산자) |

`irep = 8` **이 정확히 일치한 것이 중요합니다.** 반복 횟수는 종료 조건에 걸리는 시점이 결정하므로, 중간 계산이 조금이라도 어긋나면 반복 수가 달라집니다. 8회가 맞았다는 건 **매 반복의 수치가 IDL과 같은 궤적을 그렸다**는 뜻입니다.

### 외부에서 입력받은 값

`obs_in`, `err_in`, `response`, `response_err`, `delta_temp`는 IDL의 `aia_get_response`와
`dem_aiainterpolresponse_sites`로 계산한 값을 사용했습니다.

Kernel `ker`는 IDL의 `gaussian_function`으로 계산해 외부 입력으로 사용했습니다. 알고리즘
본체는 double precision으로 비교했고, `gaussian_function` 자체는 별도로 float32 결과를
확인했습니다.

### 한계 (비교 보고서 명시)

검증 입력은 **합성 Gaussian DEM 한 픽셀**입니다. 여러 실관측 자료는 확인하지 않았습니다.

맵 단위 처리 함수인 **Grid-SITES**(`dem_gridsites.pro`)는 이번 Python 변환에 포함하지 않았습니다.

다봉 또는 급경사 DEM에서는 검증하지 않았습니다.

`obs_in`이 0이면 `err_in/obs_in`과 `obswt`에서 0으로 나누게 됩니다. IDL 원본도 같으므로 입력 전에 해당 채널을 처리해야 합니다.

**기록의 근거.** 코드 설명은 `idl/dem_sites.pro`와 `python/dem_sites.py`를 바탕으로 작성했다.
검증 수치는 변환 작업 당시의 비교 보고서를 따른다.

</details>

## 참고 논문

- Morgan & Pickering (2019), [SITES: Solar Iterative Temperature Emission Solver for Differential Emission Measure Inversion of EUV Observations](https://doi.org/10.1007/s11207-019-1525-4)
- Pickering & Morgan (2019), [GRID-SITES: Gridded Solar Iterative Temperature Emission Solver for Fast Differential Emission Measure Inversion](https://doi.org/10.1007/s11207-019-1526-3)

원본 `idl/readme.txt`도 이 방법을 사용할 때 위 두 논문을 인용해 달라고 요청한다.

---

## 2026-08-18 재검증·수정 기록

실관측 프레임(2011-02-15T01:49:50, disk center 128×128 = 16,384픽셀) 대조에서
`gaussian_function`의 기본 커널 폭이 IDL과 달랐다(21 vs 23). 원인은 IDL 8.6
`lib/gaussian_function.pro`의 폭 규칙 `width = 2*((CEIL(3σ)) OR 1)+1` 중
홀수 강제(`OR 1`) 단계 누락 — `ceil(3σ)`가 짝수인 σ 구간(예: nt=41 → σ=3.28)에서만
발현되어 기존 검증(커널 IDL 주입 방식)에서는 드러나지 않았다.

수정: `gaussian_function`을 IDL 라이브러리 소스 충실 구현(float32 체인 + glibc
`expf`)으로 교체, `nker` 계산의 float32 리터럴 체인 반영. 수정 후 커널은 IDL과
비트 단위로 일치하며, 16,384픽셀 solver-core·end-to-end 대조에서 dem·demerr·
obsmod·irep 전 항목 통과.
