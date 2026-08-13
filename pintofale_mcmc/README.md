# pintofale_mcmc — MCMC로 DEM 분포 살펴보기

이 코드는 DEM 후보를 조금씩 바꾸어 보면서 관측자료와 잘 맞는 후보를 모은다. 한 개의 최적해만
내놓기보다, 허용 가능한 DEM의 분포를 MCMC 체인으로 탐색하는 방법이다.

## 먼저 읽기

Python 변환에서 가장 어려운 부분은 계산식보다 **IDL과 같은 난수 흐름을 재현하는 일**이었다.
IDL의 MT19937 난수 순서, float32와 float64의 반올림, 수학 라이브러리 함수, 인수 전달에 따른
부작용까지 맞췄다. 그 결과 세 시험 세트의 짧은 체인 822단계가 비트 단위까지 같았다.

이 결과는 모든 환경에 그대로 적용되지 않는다. 비트 단위 일치는 시험한 Linux x86_64 환경과
그 환경의 수학 라이브러리를 기준으로 하며, 체인 길이도 `nsim=10~15`로 짧다. 실제 분석에서
쓰는 긴 체인의 통계적 수렴은 별도로 확인해야 한다.

`nosrch`와 격자 시작, `chi2`/`rchi`, `varsmooth`, 고정·자유 abundance 등 시험한 경로만
지원한다. LOOPY, SPLINY, ONLYRAT, MIXIE 등은 구현하지 않았고 사용하면 명시적으로 오류가
난다. IDL 원본의 동작과 맞추기 위해 원본 결함 일부도 그대로 재현했으므로 “IDL과 같음”이 곧
“모든 입력에서 안전함”을 뜻하지는 않는다.

- IDL 원본: `idl/`
- Python 변환본: `python/`

<details>
<summary>자세한 구현 및 검증 기록 열기</summary>

**차례**

1. 검증 결과 요약
2. 무엇을 푸는가 — Metropolis MCMC
3. 난수 스트림을 비트 단위로 재현하기
4. 초월함수 — 함수마다 정밀도 규약이 다르다
5. float32 / float64 승격 규칙
6. IDL 참조 전달과 `xfrac`
7. 충실히 재현한 원 배포본 결함 4건
8. 검증 결과와 변환 범위
9. 부록 — 왜 SSW가 아니라 원 배포 압축 파일인가

## 1. 검증 결과 요약

변환 작업 초기에는 MCMC라서 통계적인 비교만 가능할 것으로 예상했습니다. 그러나 난수 흐름과
계산 정밀도를 맞춘 뒤에는 시험용 체인을 비트 단위로 비교할 수 있었습니다.

**그 판정은 틀렸고, 이 패키지가 그것을 증명했습니다.**

**3개 세트 × 274 체크포인트 = 822/822 비트 단위까지 같음. 회귀테스트 36/36. 허용한 차이 0건.**

IDL 8.6의 기본 난수생성기는 **MT19937**입니다. 시험에서는 시드와 난수 소비 순서를 맞추면
NumPy의 기존 `RandomState`와 같은 난수 흐름이 나온다는 점도 확인했습니다.

MCMC는 random sampling을 사용하지만 pseudo-random number generator 자체는 deterministic합니다.
같은 seed에서 같은 난수열을 사용하고 이후 계산도 같으면 체인 전체를 재현할 수 있습니다.

> 단, 난수열만 같아서는 충분하지 않습니다. 각 단계의 likelihood가 1 ULP만 달라도 수락/거절
> 판정이 달라질 수 있습니다. 그래서 난수, 초월함수, 부동소수점 승격 규칙을 모두 IDL과
> 맞췄습니다. 자세한 내용은 3~5장을 참고하십시오.

## 2. 무엇을 푸는가 — Metropolis MCMC

관측은 **스펙트럼 선 flux**이고, 미지수는 `log10(DEM)`과 (선택적으로) `log10(원소비)` 입니다.

다른 solver들이 **하나의 최적해**를 찾는다면, MCMC는 **사후분포에서 표본을 뽑습니다.** 그 표본들의 분포가 곧 불확실도입니다.

### Metropolis 루프의 핵심 mcmc_dem.py 268–308

```
rn = rng.randomn(npar)              # ① 제안용 정규난수
ru = rng.randomu(2*npar + 1)        # ② 파라미터 선택 + 수락 판정용 균등난수
    tmpar = allpar + rn * sigpar    #   proposal: 현재 상태를 perturb
    tmpff, prb, tmpdem, tmpabnd = mcmc_dem_only(...)   # ③ 제안의 우도
    accept_fn   = min(oldprb - prb, 0.0)               # ④
    metro_check = alog_f32(ru[jpar + npar])            # ⑤ log(u)
    accepted    = metro_check < accept_fn              # ⑥
    if accepted:
        allpar, dem, oldprb, abnd = tmpar, tmpdem, prb, tmpabnd
```

**④~⑥ 수락 판정**

`prb`는 χ² 계열 값이라 **작을수록 좋습니다.**

```
제안이 더 좋으면   oldprb − prb > 0  →  accept_fn = 0   →  log(u) < 0 은 항상 참  →  무조건 수락
제안이 더 나쁘면   oldprb − prb < 0  →  accept_fn < 0   →  확률 exp(accept_fn) 로 수락
```

표준 Metropolis 규칙입니다. 나빠지는 방향도 **가끔 받아들이기** 때문에 지역 최소점에 갇히지 않고 분포 전체를 탐색합니다.

⑥의 부등식은 체인 일치 여부를 결정하는 핵심입니다. `metro_check`나 `accept_fn`의 마지막 비트가
달라 수락/거절 판정이 한 번만 바뀌어도, 그다음부터 IDL과 Python 체인은 서로 다른 상태를
따라갑니다. 평균값이 비슷한지만 확인해서는 같은 체인을 재현했다고 볼 수 없습니다.

## 3. 난수 스트림을 비트 단위로 재현하기

이 대응 규칙은 통제된 시험으로 확인했으며 `python/idl_random.py`에 구현했습니다.

| IDL 호출 | Python 등가 | 워드 소비 |
|---|---|---|
| randomu(s,n,/double) | rs.random_sample(n) | 2 워드/값 (53-bit) |
| randomu(s,n) | float32(word / 2**32) | 1 워드/값 |
| randomu(s,n,/long) | word >> 1 | 1 워드/값 |
| randomn(s,n,/double) | rs.standard_normal(n)  + 캐시 클리어 | polar over 53-bit |
| randomn(s,n) | polar over 1-워드 균등난수 + 특정 반올림 배치 | 1 워드/균등난수 |

시드는 INT·LONG·DOUBLE 중 무엇으로 줘도 스트림이 같습니다.

### 차이 ① — float 정규난수의 반올림 위치

`float32(double 결과)`가 **아닙니다.** 전수 탐색(반올림 배치 2⁷ × 단정밀도 커널 24종)으로 **유일해**를 찾았습니다:

```
u1, u2 = word / 2**32                       ← float64 유지 (여기서 f32로 반올림하지 않음)
v1, v2 = float32(2*u - 1)                   ← 여기서 float32
rsq    = v1*v1 + v2*v2                      ← float32 연산
fac    = float32( sqrt(-2*log(rsq)/rsq) )   ← float64로 계산 후 마지막에 1회 반올림
out    = float32(v2 * fac),  spare = float32(v1 * fac)
```

중간 단계를 **전부 float32로 하거나 전부 float64로 하면** 100개 중 9~77개가 **±1 ulp** 어긋납니다. "어느 지점에서 반올림하는가"가 정확히 한 가지로 정해져 있습니다.

### 차이 ② — 남은 편차(spare deviate)를 처리하는 순서

polar 방법은 정규난수를 **2개씩** 만듭니다. 홀수 개를 요청하면 1개가 남는데:

|  | 호출 경계에서 남은 1개 | numpy 기본 |
|---|---|---|
| randomn(s,n)  (float) | 유지 — 중간에 randomu 가 끼어도 살아남음 | 유지 ✔ |
| randomn(s,n,/double) | 폐기 | 유지 ✘ → 클리어 필요 |

근거 (각 4개 중간값로 교차확인):

**float**: `randomn(s,1)` ×6 == `randomn(s,6)`. `randomn(s,2)` ×3, `randomn(s,3)` ×2도 동일

**double**: `randomn(s,1)` ×6 → 단일 호출의 **0, 2, 4, 6, 8, 10**번째. `randomn(s,3)` ×2 → **0,1,2** 그리고 **4,5,6** — **3번 인덱스가 버려짐**

`mcmc_dem.pro`는 `randomn(seed, mno)` 처럼 **가변·홀수 길이**로 호출합니다. 이 규칙을 틀리면 **체인이 즉시** **분기**합니다. 실사용 RNG 호출 4곳(505·506·524·529줄)은 전부 기본 생성기 + **float** 경로라 위 검증 범위 안에 있습니다.

## 4. 초월함수 — 함수마다 정밀도 규약이 다르다

이것도 예상 밖의 발견이었습니다. IDL은 **함수마다 다른 방식**으로 계산합니다.

| IDL | 실제 규약 | 근거 |
|---|---|---|
| ALOG10(float) | 플랫폼 log10f  (단정밀도 libm) | libm 30/30 일치. f32(log10(f64)) 는 3개 틀림 |
| ALOG(float) | 플랫폼 logf | metro_check  전량 일치 |
| EXP(float) | double 계산 후 float32 반올림 | 202/202 일치. expf 는 프로파일당 1개 틀림 |
| 10.D^x | 플랫폼 pow(10,x) | 30/30 일치. numpy ** 는 7개 틀림 |

> **⚠** `ALOG` **·** `ALOG10` **은 단정밀도 libm인데** `EXP` **는 아닙니다.** 같은 IDL 안에서도 규약이 통일돼 있지 않으므로, 함수별로 실측해서 맞춰야 합니다.

### 해결 — ctypes로 libm 직접 호출

```
_libm = ctypes.CDLL(ctypes.util.find_library("m") or "libm.so.6")
for _n in ("log10f", "logf"):
    _f = getattr(_libm, _n)
    _f.restype  = ctypes.c_float
    _f.argtypes = [ctypes.c_float]
```

numpy의 float32 커널과 `**`는 **어느 쪽과도 다르므로**, 시스템 libm을 그대로 불러 씁니다

이 처리는 `python/idl_compat.py`에 구현했습니다.

**환경에 따른 제한:** 이 코드는 운영체제의 `libm`에 의존합니다. 다른 OS나 C 라이브러리에서는
결과가 달라질 수 있습니다. 비트 단위 일치는 IDL 기준값을 만든 것과 같은 Linux x86_64
환경에서만 확인했습니다.

## 5. float32 / float64 승격 규칙

| 지점 | 내용 |
|---|---|
| alog10(abnd) | abnd 가 FLTARR → 결과가 float32. double dem  로그와 concat되며 승격되지만 값은 이미 float32로 반올림된 뒤 |
| xdem=findgen(n)*0.1+mindem | findgen*0.1 은 float32, +mindem (double)에서 double 승격. numpy는 스칼라 덧셈에서 float32를 유지하므로 명시적 승격 필요 (4.8e-8 오차의 원인) |
| ebound=0.9 | BOUND 미지정 시 FLOAT 리터럴. ebound*mou 가 float32 곱이라 0.9f × 10 → 9.0 으로 반올림. int64와 곱해 float64로 승격시키면 신뢰구간이 2.6e-8 어긋남 |
| aberr  vs demerr | fltarr(nab,2) (float32) vs dblarr(nt,2) (double) — 같은 함수의 출력인데 정밀도가 다름 |
| TOTAL | IDL은 순차 누적, numpy sum 은 pairwise → 순차 루프로 대체 |
| INTERPOL | X가 감소 배열일 수 있음( reverse(findgen(mol)) ). np.searchsorted  사용 불가 |

| 0.9f × 10 → |  |
|---|---|

|  | 9.0 |
|---|---|

`0.9f × 10 → 9.0` 같은 항목이 이 작업의 성격을 보여줍니다. **"대략 9"가 아니라 "정확히 어느 비트로 9인가"**를 맞춰야 체인이 갈라지지 않습니다.

## 6. IDL 참조 전달과 `xfrac`

비교 보고서가 **"가장 까다로웠던 지점"**으로 꼽은 부분입니다.

### 일어나는 일

```
mcmc_dem:       xfrac = 0*wvl*1.0          ← 0으로 채운 배열 (!)
                ff = lineflx(...) / xfrac  ← 0으로 나눔 → ff = +Inf
                oldprb = bestprb = +Inf
첫 스텝:        accept_fn = min(Inf − prb, 0) = 0
                → metro_check = log(u) < 0 은 항상 참
                → 첫 제안은 무조건 수락
```

여기까지는 단순한 결함입니다. 문제는 다음입니다:

```
mcmc_dem_only:  xfrac 를 키워드로 받아  xfrac = 0*wvl + 1  로 덮어씀
```

> **⚠ IDL 키워드는 참조전달(참조 전달)입니다.** 따라서 `mcmc_dem_only`가 `xfrac`를 1로 바꾸면 **호출자인** `mcmc_dem` **의** `xfrac` **도 1이 됩니다.** Python이었다면 지역 변수 재할당으로 끝났을 일이 IDL에서는 호출자 상태를 바꿉니다.

### 결과 — 시작 모드에 따라 같은 코드가 다르게 동작한다

| 모드 | mcmc_dem_only  사전 호출 | 루프 진입 시 ff |
|---|---|---|
| /nosrch  (set1, set3) | 없음 → xfrac 는 0 | +Inf → 첫 제안 무조건 수락 |
| 그리드 탐색 ON (set2) | 탐색이 호출 → xfrac 가 1로 바뀜 | 유한값 → 정상 판정 |

`set2`는 이 차이를 확인하기 위해 grid-search 시작 경로를 사용했습니다. `set1`만 시험했다면
`xfrac = 0` 경로만 실행되어 이 참조 전달 효과를 확인할 수 없었습니다.

## 7. 충실히 재현한 원 배포본 결함 4건

변환의 원칙은 **"고치지 않고 그대로 재현한다"**입니다. 목표가 "더 나은 코드"가 아니라 **"IDL과 같은 답"**이기 때문입니다.

| # | 결함 | 영향 |
|---|---|---|
| 1 | xfrac = 0*wvl*1.0  — 다른 곳은 전부 0*wvl+1  → 0으로 나눔 | 초기 ff = Inf, 첫 제안 무조건 수락 (6장) |
| 2 | `mcmc_dem_only`에 `pardem=lscal`을 넘기지만 해당 범위에 `lscal`이 없음(인자 이름은 `scale`) | 단계별 평활화 크기가 상수 3으로 고정되고 `findscale` 결과가 반복문 안에서 무시됨 |
| 3 | lineflx 의 nh_ne(oz)  — nh_ne 는 온도 크기인데 파장 인덱스로 접근 | Z==1  라인이 있을 때만 발현. 변환본은 그 경로를 추측하지 않고 NotImplementedError |
| 4 | `storpar(istor)` 범위 검사가 바깥쪽 반복 사이에만 있음 | `nburn*nbatch > nsim`이면 IDL이 중단됨. 시험 입력을 만들 때 이 조건을 피했음 |

**결함 2가 특히 중요합니다.** `findscale`로 애써 계산한 평활화 스케일이 MCMC 루프 안에서는 쓰이지 않고 상수 3이 대신 들어갑니다. **이 코드로 나온 기존 연구 결과들에도 그대로 적용되는 사실**입니다.

**결함 3의 처리 방식이 모범적입니다.** "아마 이런 의도였겠지"라고 고쳐 쓰면 IDL과 다른 답이 나오고, 그 차이를 나중에 아무도 추적할 수 없습니다. **명시적 예외를 던져 "여기는 검증되지 않았다"를 실행 시점에** **알리는** 편이 안전합니다.

## 8. 검증 결과와 변환 범위

### 세트 구성

| 세트 | 구성 | npar | 체크포인트 |
|---|---|---|---|
| set1 | /nosrch, 원소비 고정, nsim=10 · nbatch=5 · nburn=2 | 12 | 274/274 |
| set2 | 그리드 탐색 ON, 원소비 고정 | 12 | 274/274 |
| set3 | /nosrch, ABRNG로 원소비 해동, nsim=15 · nbatch=5 · nburn=3 | 42 | 274/274 |

**합계 822/822 비트 단위까지 같음, 허용한 차이 0건, 회귀테스트 36/36.**

### 비교 대상

표면적인 최종 출력만이 아니라 **체인 내부까지** 대조했습니다:

**setup 전량** — `allpar`, `ipar`, `opar`, `aparmx/mn`, `sigpar`, `emenv`, `cenv`, `xdem`, `xab`, `rngD`, `rngA` **배치별 난수 20세트** — `rn`, `ru`

**단계별 상태 20개** — `tmpar`, `tmpdem`, `tmpabnd`, `tmpff`, `prb`, `accept_fn`, `metro_check`, `allpar`

**최종 출력 전량** — `dem`, `demerr`, `abund`, `aberr`, `simprb`, `simdem`, `simabn`, `simflx`, `storpar`, `storidx`, `stordem`, `bestpar`

`accept_fn` **과** `metro_check` **를 단계별로 비교한 것**이 이 검증의 핵심입니다. 최종 DEM만 맞춰서는 "우연히 비슷한 곳에 도달"과 구분되지 않습니다. **매 단계의 수락/거절 판정이 같았다**는 것이 진짜 증거입니다.

IDL 기준값 결정론 확인: 동일 배치를 **2회 실행** → 60 중간값 / 427 변수 전부 비트 동일.

### 검증을 위한 변경 — 알고리즘 변경 아님

변환 당시 사용한 계측 스크립트는 원 배포본의 계산식을 바꾸지 않고 다음 기능만 덧붙였습니다.

1. **시드 고정** — 원 배포본은 `seed`를 초기화하지 않아 실행 시각에 따라 난수열이 달라집니다.
   중간값 저장용 실행본은 `MCMC_SEED`(기본 42)를 사용했습니다. 수치 결과에 영향을 주는 변경은 이것뿐이며,
   시드를 고정하지 않으면 IDL끼리도 같은 결과를 재현할 수 없습니다.

2. `chk_dump` 중간값 삽입 (문장 경계에만)

3. 실행 환경: `verbose=0`과 `set_plot,'z'`를 사용해 화면을 띄우지 않고 실행했습니다. IDL
   기준값은 원 배포본의 계산식을 수정하지 않은 상태에서 만들었습니다.

### IDL 기준값 입력

`abnd` ← `GETABUND('anders & grevesse')`

`lscal` ← `FINDSCALE(line)`

테스트 입력 `emis`·`truedem` — IDL의 `libm`과 NumPy의 `exp` 결과가 최종 1 ULP(2.3e-16)만큼 달라, 입력을 만드는 코드의 차이가 알고리즘 검증에 섞이지 않도록 IDL에서 계산한 값을 사용

### 변환 범위

| 구분 | 내용 |
|---|---|
| 구현 | /nosrch  및 그리드 탐색 시작, /chi2 (+ /rchi ) 우도, varsmooth 평활화, 원소비 고정/해 동, burn-in 시그마 재추정, 신뢰구간 산출 |
| 미구현 (호출 시 명시적 예 외) | LOOPY, SPLINY, ONLYRAT, MIXIE  디블렌딩, SYSDEV  emissivity jitter, SOFTLIM, SAMPENV, upper limits, savfil |
| 생략 | PRED_FLX (진행표시용, 체인 상태 무영향), ADJUSTIE (TIES 미설정 시 no-op) |

### 재현

```
python python/run_parity.py all        # 시험 자료가 있을 때 822/822 비교
python -m pytest tests/ -q             # 36 passed
```

## 9. 부록 — 왜 SSW가 아니라 원 배포 압축 파일인가

다른 11개 패키지는 서버의 SSW 설치본을 썼는데, 이 패키지만 **원 배포 압축 파일**( `PoA_pro_current.tar.gz`, sha256 `d205c2b8…` )을 씁니다.

### 이유 ① — SSW 설치본을 읽을 수 없음

서버의 SSW `poa` 설치본은 `idl/stat/` · `idl/fitting/` 디렉터리가 `drwxr--r--` **(실행권한 없음)**이라 `likeli.pro`를 **읽을 수도 컴파일할 수도 없습니다.** IDL도 같은 사용자로 실행되므로 **IDL 기준값 자체가 불가능**합니다. 서버 전체를 검색해도 읽을 수 있는 사본이 없었습니다.

### 이유 ② — 서로 다른 버전이 섞이는 것을 막기 위해

SSW 본체(2004)에 tarball의 `likeli` (이력 2008까지)를 섞으면 **실존한 적 없는 조합**이 됩니다. 특히 **"Apr08:** `/chi2` **출력 보정"**은 **Metropolis를 구동하는 바로 그 반환값**을 바꾸는 변경입니다.

그래서 **일관된 원 배포본 트리 하나**로 IDL 기준값과 변환 대상을 통일했습니다.

### 버전 델타

`solarsoft_dem/pintofale_mcmc/mcmc_dem.pro` (SSW 2004 스냅샷)와 **352줄** 차이가 있습니다. 주요 기능 델타:

`verbose` 키워드 (플로팅 차단)

`sysdev` / `sysgrp` (emissivity jitter)

`xnsig`, `stordem`

`hitpar(opar)` 버그 수정

burn-in 시그마 추정을 **히스토그램 → 분위수** 방식으로 교체 (VK Apr2012)

따라서 이 검증은 당시 받은 원 배포본의 현재 버전을 대상으로 합니다. SSW의 2004년 보관본을
쓰는 기존 연구 코드와는 동작이 다를 수 있습니다.

### 사용 시 주의사항

**원 배포본은 시드를 초기화하지 않습니다.** 재현이 필요하면 **반드시 시드를 명시**해야 합니다 — IDL에서도 마찬가지입니다

원본 결함 2 때문에 `findscale`의 smoothing scale은 MCMC loop에서 무시되고 상수 3이 사용됩니다.

원본에서는 `nburn*nbatch > nsim`이면 array 범위를 벗어나 IDL 실행이 중단됩니다(결함 4).

구현하지 않은 keyword를 사용하면 명시적인 exception이 발생합니다.

비트 단위 일치는 시험에 사용한 **Linux x86_64와 해당 `libm`** 조합에서만 확인했습니다.

검증에 사용한 체인은 `nsim=10~15`로 짧습니다. 수천~수만 단계의 실사용 체인은 검증하지 않았습니다.

**기록의 근거.** 계산 규칙은 `idl/` 원본과 `python/idl_random.py`, `python/idl_compat.py`,
`python/mcmc_dem.py`를 바탕으로 작성했다. 검증 수치는 변환 작업 당시의 비교 보고서를 따른다.

</details>

## 참고 논문

- Kashyap & Drake (1998), [Markov-Chain Monte Carlo Reconstruction of Emission Measure Distributions: Application to Solar Extreme-Ultraviolet Spectra](https://ui.adsabs.harvard.edu/abs/1998ApJ...503..450K/abstract) ([저자 보관 PDF](https://hea-www.cfa.harvard.edu/AstroStat/etc/mcmc_dem_kashyap%2Bdrake_1998.pdf))
