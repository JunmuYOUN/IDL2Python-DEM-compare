# aschwanden aia_teem 코드 설명

Aschwanden 2011/2013 — AIA 단일 Gaussian T/EM 적합

SpaceAI-DEM / `work/aschwanden_aia_teem_parity/`  · 2026-07-27

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 한눈에 — 역산이 아니라 순방향 적합
2. DEM 모양을 Gaussian 하나로 가정한다
3. `build_flux_table`  — 답안지 미리 만들기
4. `fit`  — 브루트포스 격자 탐색
5. float32 — 이 변환에서 가장 중요한 발견
6. IDL ↔ Python 대조
7. 검증 결과와 범위

## 1. 한눈에 — 역산이 아니라 순방향 적합

지금까지의 solver들은 관측에서 DEM을 **거꾸로 풀어냈습니다**(역산). Aschwanden의 방법은 접근 자체가 다릅니다.

**DEM의 모양을 미리 정해놓고, 파라미터만 맞춘다.**

"DEM은 log T 공간에서 **Gaussian 하나**일 것"이라 가정하고, 그 Gaussian의 **중심·폭·높이 3개 값**만 찾습니다. 역행렬도, 정규화도, 반복 수렴도 필요 없습니다 — **가능한 조합을 전부 대입해보고 제일 잘 맞는 걸** **고릅니다.**

결과적으로 문제의 성격이 뒤집힙니다:

|  | 역산 계열 (simple_reg_dem 등) | aia_teem |
|---|---|---|
| 미지수 | 41개 (온도격자) | 3개 (Te, σ, EM) |
| 관측 | 6개 | 6개 |
| 문제 성격 | under-determined → ill-posed | over-determined → 정상 적합 문제 |
| 정규화 | 필수 | 불필요 — 모양 가정 자체가 강한 제약 |
| χ² | 목표값(1)을 겨냥 | 진짜 최소화 대상. 자유도 = 6−3 = 3 |

대가는 명확합니다. **실제 DEM이 Gaussian 하나가 아니면 결과는 의미가 없습니다.** 다봉 분포나 긴 꼬리를 가진 플라즈마는 표현할 수 없습니다.

## 2. DEM 모양을 Gaussian 하나로 가정한다

```
DEM(logT) = EM · exp( −(logT − Te)² / (2σ²) )
                       └ 중심 ┘      └ 폭 ┘
            └ 높이 ┘
```

| 파라미터 | 뜻 | 탐색 범위 (검증 세트) |
|---|---|---|
| Te | 봉우리 중심 온도 (log T) | 29개 격자점 |
| σ | 온도 폭 | 20개 격자점 |
| EM | 총 emission measure (높이) | 연속 — 해석적으로 즉시 계산 |

핵심 요령은 `EM` 을 탐색 대상에서 뺀 것입니다. `(Te, σ)` 가 정해지면 **모양이 확정**되고, 그때 남은 `EM` 은 단순한 크기 배율이라 **나눗셈 한 번**으로 구해집니다. 덕분에 탐색 공간이 3차원이 아니라 **2차원(29×20 = 580조합)**으로 줄어듭니다.

## 3. build_flux_table — 답안지 미리 만들기

```
def build_flux_table(resp, telog, dte, tsig):
    """flux[nte,nsig,nwave] for unity-EM single-Gaussian DEMs. resp is [nte,nwave]."""
    nte, nwave = resp.shape
    nsig = tsig.size
    flux = np.zeros((nte, nsig, nwave))
    for i in range(nte):                                            # 중심 Te
        for j in range(nsig):                                       # 폭 σ
            em_kelvin = np.exp(-(telog - telog[i])**2 / (2.0 * tsig[j]**2))   # ①
            for iw in range(nwave):                                 # 채널
                flux[i, j, iw] = np.sum(resp[:, iw] * em_kelvin * dte)        # ②
    return flux
```

**① Gaussian 하나 만들기**

`telog[i]` 를 중심으로, 폭 `tsig[j]` 인 Gaussian을 온도격자 위에 그립니다. **높이는 1**(unity EM)입니다.

**② 그 DEM이 각 채널에 만드는 밝기 계산**

```
flux[i,j,iw] = Σ_t ( 응답[t,iw] × Gaussian[t] × 온도폭[t] )
```

순방향 적분 `∫ K(T)·DEM(T) dT` 를 이산합으로 계산한 것입니다.

**결과** `flux[i,j,:]` **= "중심** `telog[i]` **, 폭** `tsig[j]` **인 Gaussian DEM(높이 1)이 6채널에 만드는 밝기 6개".**

이걸 580개 조합에 대해 **미리 전부 계산해 표로 저장**합니다. 픽셀 수백만 개를 처리할 때 매번 적분하지 않고 표에서 꺼내 쓰기 위한 **룩업 테이블**입니다.

이 함수는 `float64` 로 계산합니다 (다음 장의 `fit` 과 대조됩니다). 검증에서 `01_flux`  체크포인트가 `1e-5`  기준을 통과했습니다.

## 4. fit — 브루트포스 격자 탐색

```
for j in range(nyy):
  for i in range(nxx):                                  # 픽셀마다
    flux_obs = images[i, j, :]                          # 관측 6개
    counts   = flux_obs * texp_                          # ① 총 광자수
    noise    = np.sqrt(counts) / texp_                   # ② Poisson 잡음
    chimin   = 9999.0
    for k in range(nte):                                # ③ 580조합 전수 탐색
      for l in range(nsig):
        flux_dem1 = flux[k, l, :]                       #    표에서 꺼냄
        em1       = np.sum(flux_obs) / np.sum(flux_dem1) # ④ EM 즉시 결정
        flux_dem  = flux_dem1 * em1
        chi = np.sqrt(np.sum((flux_obs - flux_dem)**2 / noise**2) / (nwave - nfree))  # ⑤
        if chi <= chimin:                               # ⑥ 최소값 갱신
            chimin = chi
            em_best = np.log10(em1); te_best = telog[k]; sig_best = tsig[l]
```

### ① ② Poisson 잡음 모델 54–55

```
counts = flux_obs * texp_          # DN/s × 초 = 총 카운트
noise  = np.sqrt(counts) / texp_   # √N 잡음을 다시 DN/s로 환산
```

광자 계수의 통계적 요동은 **√N**입니다. 노출시간을 곱해 총 카운트를 얻고, 잡음을 구한 뒤 다시 나눠 원래 단위로 되돌립니다.

> **⚠ 순수 photon noise만 고려합니다.** 읽기잡음·다크·보정 불확실도는 들어있지 않습니다. simple_reg_dem이 `errors`  배열을 외부에서 받는 것과 대조적으로, 이 코드는 **오차모델을 내부에 하드코딩**합니다.

### ④ EM을 해석적으로 결정 61

```
em1 = np.sum(flux_obs) / np.sum(flux_dem1)
```

관측 6개의 합을 모델 6개의 합으로 나눕니다. 모양이 맞다면 이 비율이 곧 크기 배율입니다.

> **⚠ 이것은 최소제곱 해가 아닙니다.** 오차 가중 최소제곱이라면 `Σ(obs·model/σ²) / Σ(model²/σ²)` 가 되어야 하는데, 여기서는 **단순 총합 비율**을 씁니다. 즉 **모든 채널을 동등하게 취급**해 EM을 정하고, 잡음 가중은 그 뒤 χ² 계산에서만 반영됩니다. 원본의 설계이며 변환본은 이를 그대로 재현했습니다.

### ⑤ reduced chi 63

```
chi = np.sqrt( Σ (flux_obs − flux_dem)² / noise² / (nwave − nfree) )
                                              nwave=6, nfree=3 → 자유도 3
```

**이름에 주의.** 바깥에 `sqrt` 가 있으므로 이 값은 reduced χ**²**가 아니라 **reduced χ 자체**입니다. 값이 1이면 χ²=1과 같은 뜻이지만, 다른 solver의 χ² 값과 직접 비교할 때는 제곱 관계를 감안해야 합니다.

그리고 **자유도가 진짜로 양수(6−3=3)**입니다. simple_reg_dem의 "reduced χ²"가 엄밀히는 자유도 미정이었던 것과 달리, 여기서는 통계적으로 정당한 reduced χ입니다.

### ⑥ 동점 처리 — 변환에서 놓치기 쉬운 부분 64

```
if chi <= chimin:              # IDL 'le': ties keep the later (k,l)
```

> **⚠** `<` **가 아니라** `<=` **입니다.** IDL 원본이 `le` 를 쓰기 때문입니다.

같은 `chi`  값이 두 조합에서 나오면 **나중에 나온** `(k,l)` **이 채택**됩니다. `<` 로 바꾸면 앞의 것이 남아 **결과가 달라집니다.** 브루트포스 탐색에서 동점은 드물지 않으므로, 이 한 글자가 parity를 좌우합니다.

### 왜 이 방식이 정확한 parity를 낼 수 있나

이 탐색은 **이산 격자 위의 argmin**입니다. 최적화 알고리즘이 아니라 **중첩 for문**이죠.

경로 의존성이 없다는 뜻입니다. Powell·MPFIT·Brent 같은 옵티마이저는 부동소수점 차이가 **탐색 경로 자체**를 바꿔 bit-parity가 불가능한데, 브루트포스는 순서가 고정되어 있어 **같은 입력이면 반드시 같은 출력**이 나옵니다.

이 저장소에서 **xrt_iterative(MPFIT)**와 **vdem(GCV/Brent)**이 최적화 루프를 "경계"로 문서화한 반면, aia_teem이 waiver 없이 ALL PASS를 받은 이유가 이것입니다.

## 5. float32 — 이 변환에서 가장 중요한 발견

```
images = np.asarray(images, dtype=np.float32)
flux   = np.asarray(flux,   dtype=np.float32)
telog  = np.asarray(telog,  dtype=np.float32)
...
```

`fit`  전체가 **의도적으로 float32**입니다. IDL의 `fltarr` 가 float32이기 때문인데, 문제는 **float64로 하면 오히려** **틀린다**는 점입니다.

### 왜 정밀도를 높이면 어긋나는가

최적점 근처에서 `flux_obs` 와 `flux_dem` 은 거의 같은 값입니다. 그 차를 빼면 **큰 수에서 큰 수를 빼 작은 수를** **얻는** 상황이 됩니다 — **catastrophic cancellation**입니다.

```
flux_obs  =  1234.5678
flux_dem  =  1234.5661
차        =     0.0017     ← 유효숫자가 급격히 줄어듦
```

이 지점에서 float32와 float64는 **서로 다른 값**을 냅니다. IDL이 float32로 계산했으므로, Python이 float64를 쓰면 `chi` 가 약 **1e-3** 어긋나고, 그 차이가 `<=`  비교를 뒤집어 **다른 (k,l)을 최적으로 선택**하게 만듭니다.

**실제로 겪은 순서:** run_01에서 float64로 시도 → `chi`  발산 → 원인 분석(파국적 상쇄) → **float32로 변경** → run_02 ALL PASS.

"정밀도는 높을수록 좋다"는 직관이 **parity 검증에서는 틀립니다.** 목표가 "더 정확한 답"이 아니라 **"IDL과** **같은 답"**이기 때문입니다.

같은 교훈이 이 저장소의 다른 패키지에서도 반복됐습니다 — vdem의 `!pi` (float32), dem_sites의 커널 (float32), simple_reg_dem의 누적 순서. `work/CONVERSION_SUMMARY.md` 에 **"float32 충실성이 결정적"**으로 정리돼 있습니다.

## 6. IDL ↔ Python 대조

| IDL | Python |
|---|---|
| fltarr  (float32) | dtype=np.float32  — 반드시 유지 |
| chi le chimin | chi <= chimin  — <  아님 |
| total(resp[*,iw]*em*dte) | np.sum(resp[:,iw]*em*dte) |
| alog10 | np.log10 |
| 원소별 >  / < | np.maximum  / np.minimum |
| 중첩 for  루프 | 동일 구조 유지 (벡터화하지 않음 — 순서 보존) |

**왜 벡터화하지 않았나.** numpy로 580조합을 한 번에 계산하면 훨씬 빠르지만, **합산 순서가 바뀌어** 부동소수점 결과가 달라지고 동점 처리 순서도 무너집니다. parity가 목적이므로 **루프 구조를 그대로 유지**했습니다.

## 7. 검증 결과와 범위

### 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일세트 — 29온도 × 20σ × 6채널, 9픽셀 (off-grid Gaussian 순방향모델) |
| 판정 | ALL PASS — 17 pass / 0 fail |
| 계측 무결성 | te/em/sig/chi  전부 max\|· − base\| = 0 |
| 회귀 테스트 | 3/3 |

**"off-grid Gaussian"이 검증 설계의 핵심입니다.** 정답 Gaussian의 중심·폭을 탐색 격자점과 **일부러 어긋나게** 놓았습니다. 격자점에 딱 맞으면 `chi` 가 0이 되어 비교가 무의미해지지만, 어긋나 있으면 실제로 **580조합을 비교해 최적점을 고르는 과정**이 검증됩니다.

추가로, 오라클이 룩업테이블을 인라인으로 만든 뒤 **실제** `aia_teem_table` **루틴과 max diff = 0**임을 확인했습니다. 즉 테이블 생성 자체도 원본과 동일합니다.

### 변환 범위 (엄격 스코핑)

**변환한 것** — DEM 알고리즘 본체:

`build_flux_table`  = `aia_teem_table.pro`  전체

`fit`  = `aia_teem_map.pro`  124–158줄 (픽셀별 격자 탐색, `nx<4096`  전픽셀 분기)

**제외한 것** — I/O 스캐폴딩, DEM 알고리즘이 아님:

FITS 읽기, `aia_prep` , `make_map` , `plot_map` , `window` , `aia_lct`

24파일 패키지의 `nlfff_dem`  / `euvi_loopdem`  / STEREO 다중시점 변형 — **DATA-BLOCKED** (벡터 자기장·loop 좌표·STEREO 데이터 필요)

### 알고 써야 할 것

**DEM이 단일 Gaussian이라는 가정이 전부입니다.** 다봉·비대칭 분포는 원리적으로 복원 불가

결과는 `(Te, σ, EM)`  3개 맵입니다. 온도격자 위 DEM 곡선이 필요하면 이 3개로 다시 그려야 합니다

**탐색 격자 해상도가 결과 해상도의 상한**입니다. `Te` 는 29단계로만 나옵니다 — 연속값이 아닙니다

오차모델이 photon noise로 고정돼 있습니다

검증은 **합성 9픽셀** 규모입니다. 실관측 맵과 디스크 limb 마스킹 분기( `nx=4096` )는 미검증

**출처 표시.** 3~6장의 코드 설명은 `inbox/aia_teem_table.pro` , `inbox/aia_teem_map.pro` , `converted/aia_teem.py` 에서, 7장은 `reports/09_parity_certificate.md` 에서 직접 확인했습니다. 1·2장의 방법 분류와 다른 solver와의 대비는 `work/DEM_SOLVER_METHODS.md` 와 코드에 근거한 해석이며, Aschwanden 2011/2013 원논문 대조는 수행하지 않았습니다.
