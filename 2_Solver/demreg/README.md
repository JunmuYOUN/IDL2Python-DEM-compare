# demreg — DEM과 오차를 함께 구하는 정규화 방법

`demreg`는 관측값을 잘 재현하면서도 매끄러운 DEM을 구한다. GSVD로 regularization parameter를
고르고, 필요하면 positivity를 만족하도록 허용 χ²를 조정한다.

## 먼저 읽기

이 패키지는 새로 처음 옮긴 코드가 아니다. 기존 Python 변환본을 IDL 결과와 다시 비교했고,
그 과정에서 다음 세 가지 오류를 고쳤다.

- `logT` 구간의 가운데 값을 잘못 계산하던 문제
- 응답함수 범위 밖에서 끝값으로 고정하던 문제
- 온도 해상도 `elogt`를 잘못 계산하던 문제

주요 출력은 DEM, DEM 오차 `edem`, 온도 해상도 `elogt`, χ²다. `elogt`가 온도 격자 간격보다
훨씬 크면 그 온도에서 얻은 정보가 서로 독립적이지 않다는 뜻이다. `chisq`도 반드시 확인해야
한다. positivity 조건을 맞추는 동안 허용 χ²가 커질 수 있기 때문이다.

현재는 `dem_norm0`을 직접 넣어 쓰는 것이 안전하다. 이를 생략할 때 들어가는 초기값 자동 계산
경로는 제대로 작동하지 않거나 확인되지 않았다. 단일 픽셀과 2×2 합성 맵, 회귀 테스트 5개는
통과했다. 256픽셀 이상에서 쓰는 병렬 실행 경로는 이번 시험 규모에 포함되지 않았다.

- 사용 권장 변환본: `python/`
- 비교를 위해 남긴 기존 변환본: `python_ref/`
- IDL 원본: `idl/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 한눈에 — 4개 모듈과 전체 흐름
2. GSVD를 쓰는 이유
3. `dem_reg_map` — 판별원리로 λ 찾기
4. `dem_pix` — 2-pass 재가중과 positivity 반복
5. 오차와 온도해상도 — 이 solver만의 출력
6. **재검증으로 찾은 버그 3건**
7. 검증 결과와 확인하지 않은 경로

## 1. 한눈에 — 4개 모듈과 전체 흐름

이 패키지는 다른 폴더와 작업 방식이 다릅니다. IDL을 처음부터 새로 옮긴 것이 아니라,
**기존에 공개되어 있던 Python 변환본을 IDL 기준값과 다시 비교**했습니다. 그 과정에서 실제
오류 3건을 찾아 고쳤습니다(6장).

| 모듈 | 줄수 | 역할 |
|---|---|---|
| dn2dem_pos.py | 198 | 진입점. 온도격자·응답행렬 준비 → demmap_pos  호출 |
| demmap_pos.py | 303 | 픽셀 분배 + dem_pix (실제 solver) |
| dem_inv_gsvd.py | 65 | GSVD (일반화 특이값 분해) |
| dem_reg_map.py | 52 | 정규화 파라미터 λ 결정 |

### 전체 흐름

```
[준비 — dn2dem_pos]
   logt      = 온도 bin 중점        ← 버그 ①이 있던 곳
   tr        = 응답을 logt로 보간   ← 버그 ②가 있던 곳
   rmatrix   = 응답 × dlogT
[픽셀마다 — dem_pix]
   1st pass:  L = diag(1/√dlogT)  →  GSVD  →  λ  →  가이드 DEM
              (검증 경로에서는 dem_norm0=ones로 대체)
   반복(최대 10회):
       L = diag(√dlogT / √|dem_reg|)      ← 재가중
       GSVD → λ(판별원리) → filt → kdag → DEM
       음수가 있으면 rgt *= 1.5 하고 다시   ← 음수 방지
   출력: dem, edem(오차), elogt(온도해상도)  ← 버그 ③이 있던 곳
        chisq, dn_reg(재구성 DN)
```

**demreg의 특징 세 가지**

① **GSVD**로 λ를 바꿔가며 해를 즉시 재계산 — 반복 최소화

② **재가중(2-pass)** — 0차 제약의 편향을 줄임

③ **오차와 온도해상도를 함께 출력** — 이 저장소에서 가장 풍부한 진단 정보

## 2. GSVD를 쓰는 이유

풀려는 문제는 익숙한 형태입니다:

```
minimize  ‖K·ξ − g‖²  +  λ·‖L·ξ‖²
        └ 데이터 적합 ┘   └ 정규화 ┘
```

정규방정식은 `(KᵀK + λLᵀL)·ξ = Kᵀg` 이고, simple_reg_dem처럼 λ를 고정하고 풀 수도 있습니다. **문제는 λ를** **모른다는 것**입니다.

λ 후보 수십~수백 개를 시험할 때마다 행렬을 다시 분해하면 픽셀당 계산량이 매우 커집니다.

### GSVD의 해결책

`K`와 `L`을 **동시에** 분해합니다:

```
K = U · diag(α) · W⁻¹
L = V · diag(β) · W⁻¹          (α² + β² = 1)
```

이렇게 하면 임의의 λ에 대한 해가 **닫힌 형태**로 나옵니다:

```
ξ_λ = Σ_i  [ α_i / (α_i² + λ·β_i²) ] · (g·u_i) · w_i
             └────── filt ──────┘
```

**분해는 한 번, λ 평가는 무한정.** 대괄호 안의 `filt`만 λ에 의존하므로, λ를 바꿔도 **스칼라 몇 개만 다시 계산**하면 됩니다. 이것이 demreg가 판별원리로 λ를 탐색할 수 있는 이유입니다.

코드에서는 `dem_inv_gsvd`가 이를 SVD로 우회 구현합니다:

```
AB1 = A @ inv(B)              # K·L⁻¹
u, s, v = svd(C)              # 그 SVD
beta  = 1./np.sqrt(1+s**2)    # α²+β²=1 을 만족하도록 환산
alpha = s * beta
```

`K·L⁻¹`의 특이값 `s`가 곧 `α/β` 이므로, `α² + β² = 1` 정규화를 걸면 위 식이 나옵니다.

## 3. dem_reg_map — 판별원리로 λ 찾기

```
sigs = sigmaa[:nf] / sigmab[:nf]              # α/β
maxx = max(sigs);  minx = min(sigs)**2 * 1e-2
step = (log(maxx) - log(minx)) / (nmu - 1)
mu   = exp(arange(nmu) * step) * minx          # ① 로그 등간격 후보
for kk in range(nf):
    coef = data @ U[kk,:] - sigmaa[kk]
    for ii in range(nmu):
        arg[kk,ii] = ( mu[ii]*sigmab[kk]**2 * coef
                       / (sigmaa[kk]**2 + mu[ii]*sigmab[kk]**2) )**2    # ② 잔차 예측
discr = np.sum(arg, axis=0) - np.sum(err**2) * reg_tweak                 # ③ 판별식
opt   = mu[np.argmin(np.abs(discr))]                                     # ④ 선택
```

**① 후보 생성**

특이값 비율 `α/β`의 범위를 기준으로 **로그 등간격** 후보를 만듭니다. λ는 수십 자릿수를 오갈 수 있으므로 로그 스케일이 필수입니다.

**②③ 판별식(discrepancy)**

```
discr(μ) = [μ를 썼을 때 예상되는 잔차 제곱합]  −  [잡음 수준] × reg_tweak
```

**Morozov 판별원리(discrepancy principle)**입니다. "잔차가 잡음 수준과 같아지는 λ를 고른다" — `reg_tweak = 1` 이면 **χ² = 1을 목표로 한다는 뜻**이고, simple_reg_dem이 반복의 방향·보폭으로 달성하던 것을 demreg는 **λ를 직접 탐색해** 달성합니다.

**④ 선택**

`|discr|`이 최소인 `μ`를 고릅니다. 이분법이나 근찾기가 아니라 **격자 위 argmin**이므로 **결정론적**입니다 — 결과 일치에 유리한 성질입니다.

기본 인자는 `nmu = 500` 이지만, `dem_pix`는 `nmu = 42`를 넘겨 호출합니다. 실제로 쓰이는 값은 42입니다.

## 4. `dem_pix` — 2-pass 재가중과 positivity 반복

### 왜 재가중이 필요한가

demreg의 제약행렬 `L`은 **0차**입니다 — 도함수가 아니라 DEM **값 자체**를 억제합니다.

0차 제약을 그대로 쓰면 **DEM이 큰 온도 구간이 집중적으로 벌점**을 받습니다. DEM은 온도에 따라 수십 자릿수를 오가므로, 이는 물리적으로 부당한 왜곡입니다. (simple_reg_dem이 1차 도함수 제약을 택해 이 문제를 아예 피한 것과 대조됩니다.)

### 해결 — 가이드 DEM으로 다시 가중

```
1st pass:  L = diag( 1/√dlogT )                 → 대략적인 가이드 DEM 획득
2nd pass:  L = diag( √dlogT / √|dem_reg| )      → DEM이 큰 곳은 제약을 풂
```

`dem_reg`가 분모에 있으므로, **가이드 DEM이 큰 온도에서는** `L` **이 작아져 제약이 완화**됩니다. 0차 제약의 편향을 상쇄하는 장치입니다.

**이번 검증에서는** 두 경우 모두 `dem_norm0 = ones`로 두었습니다. 즉 1st pass의 GSVD guide DEM 생성을 건너뛰고 균일한 초기 가중치로 시작했습니다. 초기 추정값을 GSVD로 만드는 경로는 **검증 범위에 포함하지 않았습니다**(7장).

### positivity 반복 258–272

```
while (ndem > 0) and (piter < max_iter):
    L = np.diag(np.sqrt(dlogt) / np.sqrt(abs(dem_reg)))     # ① 재가중
    sva, svb, U, V, W = dem_inv_gsvd(rmatrixin.T, L)        # ② GSVD
    lamb = dem_reg_map(sva, svb, U, W, dn, edn, rgt, nmu)   # ③ λ
    for kk in range(nf):
        filt[kk,kk] = sva[kk] / (sva[kk]**2 + svb[kk]**2*lamb)
    kdag = W @ (filt.T @ U[:nf,:nf])                         # ④ 유사역
    dem_reg_out = (kdag @ dn).squeeze()
    ndem = len(dem_reg_out[dem_reg_out < 0])                 # ⑤ 음수 개수
    rgt  = rgt_fact * rgt                                     # ⑥ 목표 χ²를 1.5배로
    piter += 1
```

**양수성 처리 방식이 독특합니다.** 음수를 잘라내맵(dem_sites), log로 치환하맵(simple_reg_dem) 않습니다.

대신 **목표 χ²를** `rgt_fact = 1.5` **배씩 키웁니다.** 목표 χ²가 커지면 정규화가 강해지고, 해가 더 매끄러워져 진동으로 인한 음수가 사라집니다. **"음수가 나온다는 건 정규화가 부족하다는 신호"**라는 해석입니다.

|  | ⚠ 대가: 양수해를 얻는 대신 데이터 적합도를 희생합니다. max_iter = 10 을 다 쓰면 목표 χ²가 1.5¹⁰ ≈ 58 배가 됩니다. 그런 픽셀의 DEM은 매우 매끄럽지만 데이터를 잘 설명하지 못합니다. 반환된 chisq 를 반드시 확인해야 합니다. |
|---|---|
|  | 58 |

## 5. 오차와 온도해상도 — 이 solver만의 출력

### DEM 오차 285–287

```
delxi2 = kdag @ kdag.T
edem   = np.sqrt(np.diag(delxi2))
```

`kdag`가 데이터→DEM 변환이고 데이터가 이미 오차로 나뉘어 있으므로( `dn = dnin/ednin` ), `kdag·kdagᵀ`의 대각이 곧 **DEM의 분산**입니다. 선형 오차 전파입니다.

### 온도해상도 elogt  290–303

```
kdagk = kdag @ rmatrixin.T          # 해상도 행렬 (resolution matrix)
for kk in range(nt):
    rr = np.interp(ltt, logt, kdagk[:,kk])
    hm_mask = (rr >= max(kdagk[:,kk]) / 2.)          # 반치폭 마스크
    elogt[kk] = dlogt[kk]
    if np.sum(hm_mask) > 1:
        elogt[kk] = max(ltt[hm_mask][-1] - ltt[hm_mask][0], dlogt[kk])
```

`kdagk = kdag·K` **가 완벽한 역산이라면 단위행렬**이어야 합니다. 실제로는 정규화 때문에 각 열이 **퍼진 봉우리**가 되는데, 그 **반치전폭(FWHM)**이 "이 온도의 DEM 값이 실제로는 얼마나 넓은 온도 범위를 평균한 것인가"를 말해줍니다.

즉 **온도 분해능**입니다. 이 값이 `dlogT`보다 훨씬 크면 해당 온도의 DEM 값을 독립적인
측정값으로 해석하기 어렵습니다.

이 저장소의 11개 중 **온도해상도를 계산해주는 것은 demreg뿐**입니다. simple_reg_dem 문서에서 지적한 "불확실도 출력이 없다"는 한계를 demreg는 갖지 않습니다.

## 6. 재검증으로 찾은 버그 3건

이 패키지의 가장 중요한 결과입니다. 커뮤니티 Python 포팅이 IDL과 다르게 동작하고 있었습니다.

### 버그 ① — 온도 bin 중점 공식 dn2dem_pos.py 77–81

```
포팅:  logt[i] = log10(temps[0]) + dlogt[i] * (i + 0.5)      # 균일 간격 가정
IDL :  get_edges(alog10(temps), /mean)
       = ( log10(temps[i]) + log10(temps[i+1]) ) / 2         # 국소 중점
```

두 식은 **간격이 정확히 균일할 때만** 같습니다. 그런데 `temps`는 IDL `findgen` (float32)에서 나오므로 **미세한 지터**가 있습니다.

**입력 차이 ~5e-6 → 출력 DEM 100% 차이.**

이것이 **ill-posed 문제의 본질**을 보여주는 가장 극적인 사례입니다. 조건수가 나쁜 역산에서는 입력의 6자리째 차이가 결과를 **두 배로** 바꿀 수 있습니다. "어차피 반올림 오차인데 무슨 상관이냐"는 직관이 여기서는 완전히 틀립니다.

```
# 수정
logt = 0.5 * (np.log10(temps[:-1]) + np.log10(temps[1:]))
```

### 버그 ② — 응답함수 보간의 외삽 dn2dem_pos.py 143–148

```
포팅:  np.interp(...)                 → 범위 밖은 끝값으로 고정(clamp)
IDL :  interpol(...)                  → 범위 밖은 선형 외삽(extrapolate)
```

`logt`의 **최상단 bin 중점이** `max(tresp_logt)` **를 아주 살짝 넘습니다.** 그 한 점에서 두 함수가 다른 값을 주고, 역산이 그 차이를 증폭시켰습니다.

```
# 수정
scipy.interpolate.interp1d(tresp_logt, truse[:,i], kind='linear',
                           fill_value='extrapolate')(logt)
```

**일반화할 교훈:** `np.interp`와 IDL `interpol`은 **범위 안에서는 같지만 범위 밖에서 다릅니다.** IDL 코드를 옮길 때 **경계에서 단 한 점이라도 벗어나는지** 반드시 확인해야 합니다. (같은 지적을 xrt_iterative 문서에서도 했습니다 — 그쪽은 실측 확인이 안 된 잠재 위험으로 남아 있습니다.)

### 버그 ③ — 온도해상도 elogt 계산 demmap_pos.py 186–189, 296–303

| 항목 | 포팅 | IDL |
|---|---|---|
| ltt  격자 | arange(51)/(52-1)  + 1e-8 | findgen(51)/(51-1) |
| FWHM 폭 | 범위를 /2 | 나누지 않음 |
| 하한 | 없음 | dlogt 로 floor |
| 판정 조건 | > 0 | > 1  (2점 이상일 때만) |

네 곳이 모두 달랐습니다. 결과적으로 온도해상도가 **절반으로 과소평가**되고, 격자도 IDL과 어긋나 있었습니다.

```
# 수정
ltt = min(logt) + (max(logt)-min(logt)) * np.arange(51)/(51-1.0)
...
if (np.sum(hm_mask) > 1):
    elogt[kk] = max(ltt[hm_mask][-1] - ltt[hm_mask][0], dlogt[kk])
```

### 부수 수정

astropy 버전 비호환 `timing print` → `try/except`로 감쌈

`chisq` 스칼라 packaging 정합

### 수정하지 않은 것 — kdag / kdagk 전치

```
IDL   : kdag ## dn      →  저장 규약 [nf, nt]
numpy : kdag @  dn      →  저장 규약 [nt, nf]
```

두 배열이 **서로 전치**되어 저장됩니다. 중간값 비교 세트에서 제외했고, `kdagk`의 `transpose_match = true`와 **전 출력 통과**로 값의 정확성이 증명되므로 **수정하지 않았습니다.** 저장 규약 차이일 뿐 계산은 동일합니다.

## 7. 검증 결과와 확인하지 않은 경로

### 결과

| 세트 | 규모 | 통과 / 실패 | 판정 |
|---|---|---|---|
| g3 | 단일 픽셀 | 18 / 0 | 모두 통과 |
| map | 2×2 (4 px) | 18 / 0 | 모두 통과 |

입력: 합성 Gaussian DEM ( `d1=4e22`, `s1=0.15`, 중심 단일 `6.5` / 맵 `[6.3, 6.5, 6.7, 6.9]` )

응답: `demreg/idl/aia_resp.dat` (저장된 `aia_get_response` ), 코로나 6채널

격자: `temps = 10^(5.7 + findgen(30)/20)` → **nt = 29, nf = 6**

tolerance: 출력 double **rtol 1e-8** / atol 1e-6~1e-8

중간값 저장 전후 결과: 3세트 `max|dem − base| = 0`

회귀 테스트 **5/5** (IDL 기준값 고정 + **3버그 가드**)

회귀 테스트에 **"3버그 가드"**가 포함된 것이 중요합니다. 같은 버그가 다시 들어오면 테스트가 잡아냅니다.

### 확인하지 않은 실행 경로

**초기 guess GSVD 경로와 gloci 경로는 검증되지 않았습니다.** 게다가 그 블록에는 **포팅 버그가 남아 있는** **것이 관측**되었습니다:

• `if (test_dem_reg).shape[0] == nt` — `test_dem_reg`는 길이 1이고 `nt`는 29이므로 **항상 False** → 블록 전체가 **현재 작동하지 않는 코드**

• 그 안에 `np.arrange` (오타, `arange` ), `fcofmx` (오타, `fcofmax` ) 등이 존재

• `kdag` 조립식도 본 루프와 다름 ( `U.T@filt` vs `filt.T@U` )

`dem_norm0 = ones` **로 쓰는 한 무해**하지만, 그 경로를 활성화하면 즉시 오류가 납니다.

### 사용 시 주의사항

`dem_norm0` **을 지정하고 쓰세요.** 미지정 시 현재 작동하지 않는 코드 경로로 들어갑니다

`chisq`를 반드시 확인해야 합니다. Positivity 반복 과정에서 목표 χ²가 최대 58배까지 커질 수 있습니다.

`elogt`가 `dlogT`보다 훨씬 크면 그 온도의 DEM 값을 독립적인 정보로 해석하기 어렵습니다.

검증은 **최대 2×2픽셀의 합성 자료**로 수행했습니다. 여러 실관측 자료는 확인하지 않았습니다.

256픽셀 이상에서는 `ProcessPoolExecutor`를 사용합니다. 이 병렬 경로는 이번 검증 규모(4픽셀 이하)에서 실행되지 않았습니다.

**기록의 근거.** 코드 설명은 `idl/`의 원본과 `python/`의 네 모듈을 바탕으로 작성했다.
오류 수정 내용은 Python 코드의 `PARITY FIX` 주석에서도 확인할 수 있다. 검증 수치는 변환 작업
당시의 비교 보고서를 따른다.

</details>

## 참고 논문

- Hannah & Kontar (2012), [Differential emission measures from the regularized inversion of Hinode and SDO data](https://doi.org/10.1051/0004-6361/201117576) — 정규화 DEM 방법 논문
- Hannah & Kontar (2013), [Multi-thermal dynamics and energetics of a coronal mass ejection in the low solar atmosphere](https://doi.org/10.1051/0004-6361/201219727) — 이 방법의 관측 적용 사례
