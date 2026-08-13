# aschwanden_aia_teem — AIA 관측에서 T와 EM 찾기

이 코드는 DEM이 **log T 축에서 Gaussian 하나**라고 가정하고, 그 모양을 가장 잘 설명하는
중심 온도 `Te`, 폭 `sigma`, 크기 `EM`을 찾는다.

## 먼저 읽기

1. 여러 `Te`와 `sigma` 조합에서 AIA 채널 밝기가 어떻게 나올지 미리 계산한다.
2. 실제 관측값과 가장 비슷한 조합을 찾는다.
3. 전체 밝기 비율로 `EM`을 정한다.

따라서 자유로운 모양의 DEM을 복원하는 코드는 아니다. 봉우리가 둘 이상이거나 한쪽으로 기운
DEM은 제대로 표현할 수 없다. 결과도 온도축 전체의 DEM 배열이 아니라 `Te`, `sigma`, `EM`
맵 세 개다. `Te`와 `sigma`는 미리 정한 격자에서 고르므로 결과의 세밀함도 그 격자보다 좋아질
수 없다.

Python 변환본은 합성 자료 9픽셀에서 17개 비교 항목과 회귀 테스트 3개를 통과했다. AIA
응답함수는 이 코드가 만들지 않으며 외부에서 넣어야 한다. 실제 관측 대형 맵과 태양 가장자리
마스크는 이번에 확인하지 않았다.

- IDL 원본: `idl/`
- Python 변환본: `python/`

아래에는 수식, 계산 순서, float32를 유지한 이유를 자세히 남겼다.

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 방법 요약 — 단일 Gaussian forward fitting
2. DEM 모양을 Gaussian 하나로 가정한다
3. `build_flux_table` — flux lookup table 만들기
4. `fit` — 전수 격자 탐색
5. float32 — 이 변환에서 가장 중요한 발견
6. IDL ↔ Python 대조
7. 검증 결과와 범위

## 1. 방법 요약 — 단일 Gaussian forward fitting

일반적인 DEM inversion과 달리, 이 방법은 DEM의 형태를 먼저 정한 뒤 관측값에 맞는
parameter를 찾습니다.

**DEM의 모양을 미리 정해놓고, 파라미터만 맞춘다.**

DEM을 log T 공간의 **단일 Gaussian**으로 가정하고 중심, 폭, 높이 세 parameter만 찾습니다.
정규화된 inverse를 푸는 대신 미리 정한 grid의 모든 조합을 평가해 최솟값을 선택합니다.

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

핵심 요령은 `EM`을 탐색 대상에서 뺀 것입니다. `(Te, σ)`가 정해지면 **모양이 확정**되고, 그때 남은 `EM`은 단순한 크기 배율이라 **나눗셈 한 번**으로 구해집니다. 덕분에 탐색 공간이 3차원이 아니라 **2차원(29×20 = 580조합)**으로 줄어듭니다.

## 3. `build_flux_table` — flux lookup table 만들기

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

`telog[i]`를 중심으로, 폭 `tsig[j]` 인 Gaussian을 온도격자 위에 그립니다. **높이는 1**(unity EM)입니다.

**② 그 DEM이 각 채널에 만드는 밝기 계산**

```
flux[i,j,iw] = Σ_t ( 응답[t,iw] × Gaussian[t] × 온도폭[t] )
```

순방향 적분 `∫ K(T)·DEM(T) dT`를 이산합으로 계산한 것입니다.

**결과** `flux[i,j,:]` **= "중심** `telog[i]` **, 폭** `tsig[j]` **인 Gaussian DEM(높이 1)이 6채널에 만드는 밝기 6개".**

이걸 580개 조합에 대해 **미리 전부 계산해 표로 저장**합니다. 픽셀 수백만 개를 처리할 때 매번 적분하지 않고 표에서 꺼내 쓰기 위한 **룩업 테이블**입니다.

이 함수는 `float64`로 계산합니다 (다음 장의 `fit`과 대조됩니다). 검증에서 `01_flux` 체크포인트가 `1e-5` 기준을 통과했습니다.

## 4. fit — 전수 격자 탐색

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

> **⚠ 순수 광자 잡음만 고려합니다.** 읽기잡음·다크·보정 불확실도는 들어 있지 않습니다. `simple_reg_dem`이 `errors` 배열을 외부에서 받는 것과 달리, 이 코드는 **오차모델을 내부에 고정**합니다.

### ④ EM을 해석적으로 결정 61

```
em1 = np.sum(flux_obs) / np.sum(flux_dem1)
```

관측 6개의 합을 모델 6개의 합으로 나눕니다. 모양이 맞다면 이 비율이 곧 크기 배율입니다.

> **⚠ 이것은 최소제곱 해가 아닙니다.** 오차 가중 최소제곱이라면 `Σ(obs·model/σ²) / Σ(model²/σ²)`가 되어야 하는데, 여기서는 **단순 총합 비율**을 씁니다. 즉 **모든 채널을 동등하게 취급**해 EM을 정하고, 잡음 가중은 그 뒤 χ² 계산에서만 반영됩니다. 원본의 설계이며 변환본은 이를 그대로 재현했습니다.

### ⑤ reduced chi 63

```
chi = np.sqrt( Σ (flux_obs − flux_dem)² / noise² / (nwave − nfree) )
                                              nwave=6, nfree=3 → 자유도 3
```

**이름에 주의.** 바깥에 `sqrt`가 있으므로 이 값은 reduced χ**²**가 아니라 **reduced χ 자체**입니다. 값이 1이면 χ²=1과 같은 뜻이지만, 다른 solver의 χ² 값과 직접 비교할 때는 제곱 관계를 감안해야 합니다.

그리고 **자유도가 진짜로 양수(6−3=3)**입니다. simple_reg_dem의 "reduced χ²"가 엄밀히는 자유도 미정이었던 것과 달리, 여기서는 통계적으로 정당한 reduced χ입니다.

### ⑥ 동점 처리 — 변환에서 놓치기 쉬운 부분 64

```
if chi <= chimin:              # IDL 'le': ties keep the later (k,l)
```

> **⚠** `<` **가 아니라** `<=` **입니다.** IDL 원본이 `le`를 쓰기 때문입니다.

같은 `chi` 값이 두 조합에서 나오면 **나중에 나온** `(k,l)` **이 채택**됩니다. `<`로 바꾸면 앞의 것이 남아 **결과가 달라집니다.** 전수 탐색에서 동점은 드물지 않으므로, 이 한 글자가 결과 일치를 좌우합니다.

### Grid search의 재현성

이 탐색은 **이산 격자 위의 argmin**이며 중첩된 `for` 문으로 모든 후보를 같은 순서로 확인합니다.

경로 의존성이 없다는 뜻입니다. Powell·MPFIT·Brent 같은 최적화 방법은 부동소수점 차이가 **탐색 경로 자체**를 바꿔 비트 단위 일치가 불가능한데, 이 방식은 순서가 고정되어 있어 **같은 입력이면 반드시 같은 출력**이 나옵니다.

`xrt_iterative`의 MPFIT이나 `vdem`의 GCV/Brent와 달리 탐색 순서가 고정되어 있어, 이 코드에서는 모든 비교 항목이 허용한 차이 없이 일치했습니다.

## 5. float32 — 이 변환에서 가장 중요한 발견

```
images = np.asarray(images, dtype=np.float32)
flux   = np.asarray(flux,   dtype=np.float32)
telog  = np.asarray(telog,  dtype=np.float32)
...
```

`fit` 전체가 **의도적으로 float32**입니다. IDL의 `fltarr`가 float32이기 때문인데, 문제는 **float64로 하면 오히려** **틀린다**는 점입니다.

### float64에서 결과가 달라지는 이유

최적점 근처에서 `flux_obs`와 `flux_dem`은 거의 같은 값입니다. 그 차를 빼면 **큰 수에서 큰 수를 빼 작은 수를** **얻는** 상황이 됩니다 — **catastrophic cancellation**입니다.

```
flux_obs  =  1234.5678
flux_dem  =  1234.5661
차        =     0.0017     ← 유효숫자가 급격히 줄어듦
```

이 지점에서 float32와 float64는 **서로 다른 값**을 냅니다. IDL이 float32로 계산했으므로, Python이 float64를 쓰면 `chi`가 약 **1e-3** 어긋나고, 그 차이가 `<=` 비교를 뒤집어 **다른 (k,l)을 최적으로 선택**하게 만듭니다.

**실제로 겪은 순서:** run_01에서 float64로 시도 → `chi` 발산 → 원인 분석(파국적 상쇄) → **float32로 변경** → run_02 모두 통과.

"정밀도는 높을수록 좋다"는 직관이 **결과 일치 검증에서는 틀립니다.** 목표가 "더 정확한 답"이 아니라 **"IDL과** **같은 답"**이기 때문입니다.

같은 현상은 다른 패키지에서도 확인됐습니다. 예를 들어 vdem의 `!pi`, dem_sites의 계산 커널,
simple_reg_dem의 누적 순서에서도 float32를 원본과 맞추는 일이 중요했습니다.

## 6. IDL ↔ Python 대조

| IDL | Python |
|---|---|
| fltarr  (float32) | dtype=np.float32  — 반드시 유지 |
| chi le chimin | chi <= chimin  — <  아님 |
| total(resp[*,iw]*em*dte) | np.sum(resp[:,iw]*em*dte) |
| alog10 | np.log10 |
| 원소별 >  / < | np.maximum  / np.minimum |
| 중첩 for  루프 | 동일 구조 유지 (벡터화하지 않음 — 순서 보존) |

NumPy로 580개 조합을 한 번에 vectorize하면 더 빠르지만 합산 순서와 동점 처리 순서가
달라집니다. IDL과 같은 연산 순서를 유지하기 위해 원본의 loop 구조를 보존했습니다.

## 7. 검증 결과와 범위

### 결과

| 항목 | 값 |
|---|---|
| 세트 | 단일세트 — 29온도 × 20σ × 6채널, 9픽셀 (격자 사이 값 Gaussian forward model) |
| 판정 | 모두 통과 — 17 통과 / 0 실패 |
| 중간값 저장 전후 결과 | te/em/sig/chi  전부 max\|· − base\| = 0 |
| 회귀 테스트 | 3/3 |

**"격자 사이 값 Gaussian"이 검증 설계의 핵심입니다.** 정답 Gaussian의 중심·폭을 탐색 격자점과 **일부러 어긋나게** 놓았습니다. 격자점에 딱 맞으면 `chi`가 0이 되어 비교가 무의미해지지만, 어긋나 있으면 실제로 **580조합을 비교해 최적점을 고르는 과정**이 검증됩니다.

추가로, IDL 기준값이 룩업테이블을 인라인으로 만든 뒤 **실제** `aia_teem_table` **루틴과 max diff = 0**임을 확인했습니다. 즉 테이블 생성 자체도 원본과 동일합니다.

### 이번에 옮기고 확인한 범위

**변환한 것** — DEM 알고리즘 본체:

`build_flux_table` = `aia_teem_table.pro` 전체

`fit` = `aia_teem_map.pro` 124–158줄 (픽셀별 격자 탐색, `nx<4096` 전픽셀 분기)

**제외한 것** — I/O 입력 준비 코드, DEM 알고리즘이 아님:

FITS 읽기, `aia_prep`, `make_map`, `plot_map`, `window`, `aia_lct`

24파일 패키지의 `nlfff_dem` / `euvi_loopdem` / STEREO 다중시점 변형 — **DATA-BLOCKED** (벡터 자기장·loop 좌표·STEREO 데이터 필요)

### 사용 시 주의사항

**DEM을 단일 Gaussian으로 가정합니다.** 따라서 다봉 또는 비대칭 분포는 복원할 수 없습니다.

결과는 `(Te, σ, EM)` 맵 세 개입니다. 온도 격자 위의 DEM이 필요하면 이 세 값으로 다시 계산해야 합니다.

**결과의 해상도는 탐색 격자보다 높을 수 없습니다.** `Te`는 연속값이 아니라 29개 후보 중 하나입니다.

오차모델은 광자 잡음으로 고정되어 있습니다.

검증은 **합성 9픽셀** 규모입니다. 실관측 맵과 태양 가장자리 마스킹 분기(`nx=4096`)는 확인하지 않았습니다.

**기록의 근거.** 코드 설명은 `idl/aia_teem_table.pro`, `idl/aia_teem_map.pro`,
`python/aia_teem.py`를 바탕으로 작성했다. 검증 수치는 변환 작업 당시의 비교 보고서를 따른다.

</details>

## 참고 논문

- Aschwanden, Boerner, Schrijver & Malanushenko (2013), [Automated Temperature and Emission Measure Analysis of Coronal Loops and Active Regions Observed with SDO/AIA](https://doi.org/10.1007/s11207-011-9876-5)
- Aschwanden & Boerner (2011), [Solar Corona Loop Studies with AIA. I. Cross-Sectional Temperature Structure](https://arxiv.org/abs/1103.0228) — AIA 루프의 온도 구조를 이해하기 위한 배경 자료
