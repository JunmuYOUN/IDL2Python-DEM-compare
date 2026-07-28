# dem_sites 코드 설명

Morgan & Pickering 2019 — SITES, 반복 deconvolution DEM

SpaceAI-DEM / `work/dem_sites_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — 행렬을 뒤집지 않는 DEM
2. 준비 단계 — 가중치와 배분표
3. 반복 루프 — 알고리즘의 심장
4. 오차 추정과 마무리
5. IDL ↔ Python 대조
6. 검증 결과와 한계

## 1. 한눈에 — 행렬을 뒤집지 않는 DEM

지금까지 본 solver들은 모두 **연립방정식을 풀었습니다.** simple_reg_dem은 41×41 Cholesky를, trace_dem은 3×3 역행렬을 썼죠.

**SITES는 행렬을 아예 뒤집지 않습니다.**

대신 "지금 DEM으로 관측을 예측해보고 → 안 맞는 만큼(residual)을 온도공간에 되뿌려(back-project) 더한다"를 **수렴할 때까지 반복**합니다. Landweber / Richardson-Lucy 계열의 **반복 deconvolution**입니다.

구조상 이런 성질이 따라옵니다:

| 항목 | SITES의 처리 |
|---|---|
| 역행렬 | 없음. 행렬 곱셈과 덧셈만 사용 → 조건수 문제에서 자유로움 |
| 정규화 | 페널티 항이 아니라 매 반복마다 가우시안으로 뭉개기 + 조기 종료 |
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

`err_in/obs_in`  — **상대 측정오차**. 신호 대비 잡음이 얼마나 큰가

`response_err`  — **응답함수 불확실도**. 이 채널의 기기 보정을 얼마나 믿을 수 있나

둘을 제곱합해 루트를 씌운 뒤 역수를 취하므로, **잡음이 적고 응답이 잘 알려진 채널일수록 큰 가중치**를 받습니다.

`response_err` 를 별도 입력으로 받는 것이 이 solver의 특징입니다. 대부분의 DEM 코드는 측정오차만 보지만, SITES는 **"응답함수 자체도 틀릴 수 있다"**를 명시적으로 모델에 넣습니다.

### gres — 배분표 41–42

```
num  = response * wt[None, :]
gres = num / np.sum(num, axis=1)[:, None]
```

`np.sum(num, axis=1)` 은 **채널 방향** 합이므로, `gres` 는 각 온도 행이 **합쳐서 1**이 되도록 정규화됩니다.

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

> **⚠** `float32` **가 의도적입니다.** IDL의 `gaussian_function(σ, /norm)` 이 float32를 반환하므로 Python도 맞췄습니다. 이 커널이 코드 전체에서 **유일한 float32 요소**이고, 여기가 어긋나면 반복이 진행되며 차이가 누적됩니다.

### 미리 계산해두는 행렬 52–54

```
res2    = response * delta_temp[:, None]     # [nt,nwl] 응답 × 온도폭
totres2 = np.sum(res2 ** 2, axis=0)          # [nwl]    채널별 제곱합
res3    = res2 * gres                        # [nt,nwl] 배분표 반영
```

|  | 뜻 |
|---|---|
| res2 | res2[t,w] × DEM[t]  = 온도구간 t가 채널 w에 기여하는 양. 순방향 모델의 핵심 |
| totres2 | 채널별 정규화 상수. 반복에서 스텝 크기를 정함 |
| res3 | 순방향 응답에 배분표를 곱해둔 것. 역방향(back-projection) 전용 |

## 3. 반복 루프 — 알고리즘의 심장

```
while True:
    dem     = (obs / totres2)[None, :] * res3                    # ① 되뿌리기
    conv    = np.convolve(np.sum(dem, axis=1), ker, mode="same") # ② 뭉개기
    demmain = np.maximum(demmain + conv, 0.0)                    # ③ 누적 + 양수화
    obs     = obs_in - np.sum(demmain[:, None] * res2, axis=0)   # ④ 새 잔차
    convobs = np.sum(np.abs(obs * obswt)) / totwt                # ⑤ 수렴 측정
    irep   += 1
    if irep > 300 or convobs < convergence:
        break
```

루프 진입 직전 `obs = obs_in` (전체 관측), `demmain = 0` (빈 DEM)에서 시작합니다.

### ① 되뿌리기 (back-projection) 63

```
dem = (obs / totres2)[None, :] * res3
```

`obs` 는 "아직 설명하지 못한 관측량"입니다. 이를 `totres2` 로 나눠 스텝 크기를 정하고, `res3` 를 곱해 **온도공간에 펼칩니다.**

결과 `dem` 은 `[nt, nwl]`  — "채널 w의 잔차가 온도 t에 얼마를 기여해야 하는가"입니다. 아직 채널별로 나뉜 상태입니다.

### ② 뭉개기 (스무딩) 64

```
conv = np.convolve(np.sum(dem, axis=1), ker, mode="same")
```

`np.sum(dem, axis=1)` 로 채널 방향을 합쳐 온도축 벡터 `[nt]` 를 만든 뒤, **가우시안 커널로 흐립니다.**

**여기가 SITES의 정규화입니다.** simple_reg_dem이 `regmat` 을 방정식에 더해 진동을 억제한다면, SITES는 **매 반복마다 증분을 물리적으로 뭉개서** 같은 효과를 냅니다. 페널티 함수도, λ도, 행렬도 없습니다.

### ③ 누적하고 음수 잘라내기 65

```
demmain = np.maximum(demmain + conv, 0.0)
```

흐린 증분을 기존 DEM에 **더합니다.** 그리고 **음수가 된 부분은 0으로 잘라냅니다.**

> **⚠ IDL 함정.** 원본은 `...>0` 인데, 이는 비교가 아니라 **0과 비교해 큰 값**입니다 ( `np.maximum` ). 이 한 글자가 **양수성 강제 장치 전체**입니다 — simple_reg_dem이 log 치환으로 해결한 문제를, SITES는 매 반복 잘라내기로 해결합니다.

### ④ 새 잔차 계산 66

```
obs = obs_in - np.sum(demmain[:, None] * res2, axis=0)
```

갱신된 DEM으로 **순방향 모델**을 돌려 관측을 예측하고, 실제 관측에서 뺍니다. 남은 것이 다음 반복이 설명할 몫입니다.

`axis=0` 은 온도 방향 합이므로 결과는 `[nwl]`  — 채널별 잔차입니다.

### ⑤ 수렴 판정 67, 72

```
obswt   = wt / obs_in                                # 루프 밖에서 미리 계산 (60줄)
totwt   = np.sum(wt)
convobs = np.sum(np.abs(obs * obswt)) / totwt
if irep > 300 or convobs < convergence:
    break
```

`obs × obswt = obs/obs_in × wt` 이므로, `convobs` 는 **가중 평균 상대잔차**입니다. 기본값 `1e-2` 는 **"평균적으로** **관측값의 1% 이내로 맞으면 그만"**이라는 뜻입니다.

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

최종 DEM으로 관측을 다시 예측합니다. `obs_in` 과 비교하면 **적합도를 사용자가 직접 확인**할 수 있습니다. trace_dem에는 없던 진단 수단입니다.

## 5. IDL ↔ Python 대조

이 코드는 IDL의 `rebin`  브로드캐스팅 관용구가 많아 변환 시 **축 방향을 헷갈리기 쉽습니다.**

| IDL | Python |
|---|---|
| rebin(reform(v,1,nwl),nt,nwl)  — 행 복제 | v[None, :] |
| rebin(v_nt,nt,nwl)  — 열 복제 | v[:, None] |
| total(x, 2)  — 2번째 차원(nwl) 합 | np.sum(x, axis=1) |
| total(x, 1)  — 1번째 차원(nt) 합 | np.sum(x, axis=0) |
| convol(a, ker, /edge_zero) | np.convolve(a, ker, 'same') |
| ...>0  (원소별 max) | np.maximum(..., 0.0) |
| repeat ... endrep until cond | while True: ... if cond: break |
| temporary(demmain)  (메모리 최적화) | 불필요 — 대응물 없음 |
| arg_present(demerr)  (호출자가 요청했나) | 불필요 — 항상 계산해 반환 |
| gaussian_function(σ, /norm) | 재구현 ( 2*ceil(3σ)+1 , float32) |

> **⚠** `convol` **의 경계 처리.** IDL은 `edge_zero=1, edge_trunc=0` 로 **가장자리를 0으로 패딩**합니다. `np.convolve(..., 'same')` 이 같은 동작이라 일치했지만, IDL `convol` 은 기본값이 다르므로 다른 코드를 옮길 때는 매번 키워드를 확인해야 합니다.

## 6. 검증 결과와 한계

### 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일 픽셀 set1 — 합성 가우시안 DEM 순방향 → 관측값 |
| 판정 | ALL PASS — 19 pass / 0 fail |
| 계측 무결성 | max\|dem − base\| = 0 |
| 일치 확인 | 반복 횟수 irep = 8  일치, DEM 피크 9.567e21  일치 |
| 회귀 테스트 | 4/4 (오라클 고정 + gaussian_function  float32 일치 + max 연산자) |

`irep = 8` **이 정확히 일치한 것이 중요합니다.** 반복 횟수는 종료 조건에 걸리는 시점이 결정하므로, 중간 계산이 조금이라도 어긋나면 반복 수가 달라집니다. 8회가 맞았다는 건 **매 반복의 수치가 IDL과 같은 궤적을 그렸다**는 뜻입니다.

### 주입한 것 (변환 경계)

`obs_in` , `err_in` , `response` , `response_err` , `delta_temp`  — 오라클이 `aia_get_response` (SSW)와 `dem_aiainterpolresponse_sites` 로 만든 값

**커널** `ker`  — IDL 라이브러리 유틸 `gaussian_function` 의 출력. 주입해서 알고리즘 본체를 double 정밀도로 검증하고, `gaussian_function`  자체는 **별도로 float32 일치를 확인**했습니다

### 한계 (인증서 명시)

검증 입력은 **합성 가우시안 DEM 단일 픽셀**입니다. 실관측 다세트는 후속 과제

**Grid-SITES**(맵 단위 처리, `dem_gridsites.pro`  343줄)는 이번 범위 밖

다른 모양의 DEM(다봉·급경사)에 대한 거동은 미검증

`obs_in` 이 0이면 `err_in/obs_in` 과 `obswt` 에서 **0으로 나누기**가 발생합니다. IDL 원본도 동일 — 사용 전 필터링 필요

**출처 표시.** 2~5장의 코드 설명은 `inbox/dem_sites.pro` 와 `converted/dem_sites.py` 에서, 6장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1장의 방법 분류와 다른 solver와의 대비는 `work/DEM_SOLVER_METHODS.md` 와 코드에 근거한 해석이며, Morgan & Pickering 2019 원논문 대조는 수행하지 않았습니다.
