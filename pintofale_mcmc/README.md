# pintofale_mcmc 코드 설명

Kashyap & Drake 1998 — PINTofALE Metropolis MCMC DEM

SpaceAI-DEM / `work/pintofale_mcmc_parity/`  · 2026-07-28

> 이 문서는 변환 하네스가 만든 walkthrough PDF를 마크다운으로 옮긴 것이다. 본문의 경로 표기는
> 하네스 작업트리 기준이라 이 저장소와 다르다 — `converted/`는 이 폴더의 `python/`,
> `inbox/`는 `idl/`에 해당한다. `probes/`·`reports/`·`staging/`은 parity 검증 산출물이라
> 이 저장소에 포함하지 않았다(무엇을 왜 뺐는지는 [상위 README](../README.md) 참조).

**차례**

1. 헤드라인 — "불가능" 판정을 뒤집다
2. 무엇을 푸는가 — Metropolis MCMC
3. 난수 스트림을 비트 단위로 재현하기
4. 초월함수 — 함수마다 정밀도 규약이 다르다
5. float32 / float64 승격 함정
6. IDL 참조전달 부작용 — `xfrac`  이야기
7. 충실히 재현한 upstream 결함 4건
8. 검증 결과와 변환 범위
9. 부록 — 왜 SSW가 아니라 upstream tarball인가

## 1. 헤드라인 — "불가능" 판정을 뒤집다

이 저장소의 `work/SPECIAL_CASES.md` 는 pintofale_mcmc를 **"MCMC 확률적 → 통계검증만"**으로 분류하고 bit- parity 대상에서 제외했었습니다.

**그 판정은 틀렸고, 이 패키지가 그것을 증명했습니다.**

**3개 세트 × 274 체크포인트 = 822/822 bit-identical. 회귀테스트 36/36. waiver 0건.**

근거: IDL 8.6의 기본 난수생성기가 **MT19937**이고, 시드·소비 규칙까지 numpy의 legacy `RandomState` 와**비트 단위로 동일**하다는 사실을 통제실험으로 확인했습니다 ( `reports/00_rng_feasibility.md` ).

"확률적이니 재현 불가"라는 직관이 왜 틀렸는지가 핵심입니다. **MCMC는 확률적으로 보이지만 실제로는 결정** **론적 프로그램**입니다 — 난수생성기가 결정론적이기 때문입니다. 같은 시드에서 같은 난수열이 나오면, 체인 전체가 같은 경로를 밟습니다.

> **⚠ 단, 조건이 있습니다.** 난수열이 맞는 것만으로는 부족합니다. 매 스텝의 우도 계산이 **1 ulp라도** 갈리면 수락/거절 판정이 뒤집히고, 그 순간부터 체인이 완전히 분기합니다. 그래서 이 변환은 **난수 + 초월함수** **+ 부동소수점 승격**을 전부 IDL과 맞춰야 했습니다. 3~5장이 그 이야기입니다.

## 2. 무엇을 푸는가 — Metropolis MCMC

관측은 **스펙트럼 선 flux**이고, 미지수는 `log10(DEM)` 과 (선택적으로) `log10(원소비)` 입니다.

다른 solver들이 **하나의 최적해**를 찾는다면, MCMC는 **사후분포에서 표본을 뽑습니다.** 그 표본들의 분포가 곧 불확실도입니다.

### Metropolis 루프의 핵심 mcmc_dem.py 268–308

```
rn = rng.randomn(npar)              # ① 제안용 정규난수
ru = rng.randomu(2*npar + 1)        # ② 파라미터 선택 + 수락 판정용 균등난수
    tmpar = allpar + rn * sigpar    #   제안: 현재 상태를 흔든다
    tmpff, prb, tmpdem, tmpabnd = mcmc_dem_only(...)   # ③ 제안의 우도
    accept_fn   = min(oldprb - prb, 0.0)               # ④
    metro_check = alog_f32(ru[jpar + npar])            # ⑤ log(u)
    accepted    = metro_check < accept_fn              # ⑥
    if accepted:
        allpar, dem, oldprb, abnd = tmpar, tmpdem, prb, tmpabnd
```

**④~⑥ 수락 판정**

`prb` 는 χ² 계열 값이라 **작을수록 좋습니다.**

```
제안이 더 좋으면   oldprb − prb > 0  →  accept_fn = 0   →  log(u) < 0 은 항상 참  →  무조건 수락
제안이 더 나쁘면   oldprb − prb < 0  →  accept_fn < 0   →  확률 exp(accept_fn) 로 수락
```

표준 Metropolis 규칙입니다. 나빠지는 방향도 **가끔 받아들이기** 때문에 지역 최소점에 갇히지 않고 분포 전체를 탐색합니다.

**여기가 parity의 급소입니다.** ⑥은 **부등식 하나**입니다. `metro_check` 나 `accept_fn` 이 마지막 비트 하나만 달라도 판정이 뒤집히고, 그 스텝부터 IDL과 Python의 체인이 **완전히 다른 곳으로 갑니다.** 평균적으로 비슷한 정도로는 parity가 성립하지 않습니다.

## 3. 난수 스트림을 비트 단위로 재현하기

통제실험 결과 확정된 대응 규칙입니다 ( `converted/idl_random.py` ).

| IDL 호출 | Python 등가 | 워드 소비 |
|---|---|---|
| randomu(s,n,/double) | rs.random_sample(n) | 2 워드/값 (53-bit) |
| randomu(s,n) | float32(word / 2**32) | 1 워드/값 |
| randomu(s,n,/long) | word >> 1 | 1 워드/값 |
| randomn(s,n,/double) | rs.standard_normal(n)  + 캐시 클리어 | polar over 53-bit |
| randomn(s,n) | polar over 1-워드 균등난수 + 특정 반올림 배치 | 1 워드/균등난수 |

시드는 INT·LONG·DOUBLE 중 무엇으로 줘도 스트림이 같습니다.

### 함정 ① — float 정규난수의 반올림 배치

`float32(double 결과)` 가 **아닙니다.** 전수 탐색(반올림 배치 2⁷ × 단정밀도 커널 24종)으로 **유일해**를 찾았습니다:

```
u1, u2 = word / 2**32                       ← float64 유지 (여기서 f32로 반올림하지 않음)
v1, v2 = float32(2*u - 1)                   ← 여기서 float32
rsq    = v1*v1 + v2*v2                      ← float32 연산
fac    = float32( sqrt(-2*log(rsq)/rsq) )   ← float64로 계산 후 마지막에 1회 반올림
out    = float32(v2 * fac),  spare = float32(v1 * fac)
```

중간 단계를 **전부 float32로 하거나 전부 float64로 하면** 100개 중 9~77개가 **±1 ulp** 어긋납니다. "어느 지점에서 반올림하는가"가 정확히 한 가지로 정해져 있습니다.

### 함정 ② — 여분 편차(spare deviate)의 처리, float과 double이 반대

polar 방법은 정규난수를 **2개씩** 만듭니다. 홀수 개를 요청하면 1개가 남는데:

|  | 호출 경계에서 남은 1개 | numpy 기본 |
|---|---|---|
| randomn(s,n)  (float) | 유지 — 중간에 randomu 가 끼어도 살아남음 | 유지 ✔ |
| randomn(s,n,/double) | 폐기 | 유지 ✘ → 클리어 필요 |

근거 (각 4개 프로브로 교차확인):

**float**: `randomn(s,1)` ×6 == `randomn(s,6)` . `randomn(s,2)` ×3, `randomn(s,3)` ×2도 동일

**double**: `randomn(s,1)` ×6 → 단일 호출의 **0, 2, 4, 6, 8, 10**번째. `randomn(s,3)` ×2 → **0,1,2** 그리고 **4,5,6** — **3번 인덱스가 버려짐**

`mcmc_dem.pro` 는 `randomn(seed, mno)` 처럼 **가변·홀수 길이**로 호출합니다. 이 규칙을 틀리면 **체인이 즉시** **분기**합니다. 실사용 RNG 호출 4곳(505·506·524·529줄)은 전부 기본 생성기 + **float** 경로라 위 검증 범위 안에 있습니다.

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

numpy의 float32 커널과 `**` 는 **어느 쪽과도 다르므로**, 시스템 libm을 그대로 불러 씁니다

( `converted/idl_compat.py` ).

**이식성 대가:** 이 코드는 플랫폼 libm에 의존합니다. 다른 OS·libc에서는 결과가 달라질 수 있습니다. IDL 오라클이 돌던 **같은 linux x86_64**에서만 bit-parity가 보장됩니다.

## 5. float32 / float64 승격 함정

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

`0.9f × 10 → 9.0`  같은 항목이 이 작업의 성격을 보여줍니다. **"대략 9"가 아니라 "정확히 어느 비트로 9인가"**를 맞춰야 체인이 갈라지지 않습니다.

## 6. IDL 참조전달 부작용 — xfrac 이야기

인증서가 **"가장 까다로웠던 지점"**으로 꼽은 부분입니다.

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

> **⚠ IDL 키워드는 참조전달(pass-by-reference)입니다.** 따라서 `mcmc_dem_only` 가 `xfrac` 를 1로 바꾸면 **호출자인** `mcmc_dem` **의** `xfrac` **도 1이 됩니다.** Python이었다면 지역 변수 재할당으로 끝났을 일이 IDL에서는 호출자 상태를 바꿉니다.

### 결과 — 시작 모드에 따라 같은 코드가 다르게 동작한다

| 모드 | mcmc_dem_only  사전 호출 | 루프 진입 시 ff |
|---|---|---|
| /nosrch  (set1, set3) | 없음 → xfrac 는 0 | +Inf → 첫 제안 무조건 수락 |
| 그리드 탐색 ON (set2) | 탐색이 호출 → xfrac 가 1로 바뀜 | 유한값 → 정상 판정 |

**set2가 이걸 잡아냈습니다.** set1만 검증했다면 `xfrac = 0`  경로만 통과하고 탐색 경로에서 조용히 틀렸을 겁니다. **서로 다른 시작 모드를 세트로 나눈 설계**가 이 부작용을 드러냈습니다.

## 7. 충실히 재현한 upstream 결함 4건

변환의 원칙은 **"고치지 않고 그대로 재현한다"**입니다. 목표가 "더 나은 코드"가 아니라 **"IDL과 같은 답"**이기 때문입니다.

| # | 결함 | 영향 |
|---|---|---|
| 1 | xfrac = 0*wvl*1.0  — 다른 곳은 전부 0*wvl+1  → 0으로 나눔 | 초기 ff = Inf , 첫 제안 무조건 수락 (6장) |
| 2 | mcmc_dem_only 에 pardem=lscal 을 넘기 는데 그 스코프에 lscal 이 없음 (인자명은 scale ) | 스텝별 평활화 스케일이 상수 3으로 고정. findscale  결과 가 루프 안에서 완전히 무시됨 |
| 3 | lineflx 의 nh_ne(oz)  — nh_ne 는 온도 크기인데 파장 인덱스로 접근 | Z==1  라인이 있을 때만 발현. 변환본은 그 경로를 추측하지 않고 NotImplementedError |
| 4 | storpar(istor)  경계검사가 outer 반복 사 이에만 있음 | nburn*nbatch > nsim 이면 IDL이 크래시. 세트 구성 시 제 약으로 반영 |

**결함 2가 특히 중요합니다.** `findscale` 로 애써 계산한 평활화 스케일이 MCMC 루프 안에서는 쓰이지 않고 상수 3이 대신 들어갑니다. **이 코드로 나온 기존 연구 결과들에도 그대로 적용되는 사실**입니다.

**결함 3의 처리 방식이 모범적입니다.** "아마 이런 의도였겠지"라고 고쳐 쓰면 IDL과 다른 답이 나오고, 그 차이를 나중에 아무도 추적할 수 없습니다. **명시적 예외를 던져 "여기는 검증되지 않았다"를 실행 시점에** **알리는** 편이 안전합니다.

## 8. 검증 결과와 변환 범위

### 세트 구성

| 세트 | 구성 | npar | 체크포인트 |
|---|---|---|---|
| set1 | /nosrch , 원소비 고정, nsim=10 · nbatch=5 · nburn=2 | 12 | 274/274 |
| set2 | 그리드 탐색 ON, 원소비 고정 | 12 | 274/274 |
| set3 | /nosrch , ABRNG로 원소비 해동, nsim=15 · nbatch=5 · nburn=3 | 42 | 274/274 |

**합계 822/822 bit-identical, waiver 0건, 회귀테스트 36/36.**

### 비교 대상

표면적인 최종 출력만이 아니라 **체인 내부까지** 대조했습니다:

**setup 전량** — `allpar` , `ipar` , `opar` , `aparmx/mn` , `sigpar` , `emenv` , `cenv` , `xdem` , `xab` , `rngD` , `rngA` **배치별 난수 20세트** — `rn` , `ru`

**스텝별 상태 20개** — `tmpar` , `tmpdem` , `tmpabnd` , `tmpff` , `prb` , `accept_fn` , `metro_check` , `allpar`

**최종 출력 전량** — `dem` , `demerr` , `abund` , `aberr` , `simprb` , `simdem` , `simabn` , `simflx` , `storpar` , `storidx` , `stordem` , `bestpar`

`accept_fn` **과** `metro_check` **를 스텝별로 비교한 것**이 이 검증의 핵심입니다. 최종 DEM만 맞춰서는 "우연히 비슷한 곳에 도달"과 구분되지 않습니다. **매 스텝의 수락/거절 판정이 같았다**는 것이 진짜 증거입니다.

오라클 결정론 확인: 동일 배치를 **2회 실행** → 60 probes / 427 변수 전부 비트 동일.

### 계측 델타 — 알고리즘 변경 아님

`staging/make_staging.py` 가 upstream 소스에 **추가만** 합니다 (기존 문장 수정·삭제 없음):

1. **시드 고정** — upstream은 `seed` 를 초기화하지 않아 **시계 시딩**됩니다. 즉 **IDL에서조차 재현 불가**입니다. 계측본은 `MCMC_SEED` (기본 42)로 고정합니다. **유일한 수치 영향 델타**이며, 이것 없이는 어떤 방식의 parity도 성립하지 않습니다

2. `chk_dump`  프로브 삽입 (문장 경계에만)

3. 실행 환경: `verbose=0`  + `set_plot,'z'` 로 헤드리스 실행 — **소스 패치 없이**, 오라클은 pristine upstream 코드를 실행

### 오라클 주입

`abnd`  ← `GETABUND('anders & grevesse')`

`lscal`  ← `FINDSCALE(line)`

테스트 입력 `emis` · `truedem`  — IDL libm과 numpy의 `exp` 가 최종 **1 ulp (2.3e-16)** 다르므로, 드라이버 특성이 알고리즘 검증을 오염시키지 않도록 오라클 값 사용

### 변환 범위

| 구분 | 내용 |
|---|---|
| 구현 | /nosrch  및 그리드 탐색 시작, /chi2 (+ /rchi ) 우도, varsmooth 평활화, 원소비 고정/해 동, burn-in 시그마 재추정, 신뢰구간 산출 |
| 미구현 (호출 시 명시적 예 외) | LOOPY , SPLINY , ONLYRAT , MIXIE  디블렌딩, SYSDEV  emissivity jitter, SOFTLIM , SAMPENV , upper limits, savfil |
| 생략 | PRED_FLX (진행표시용, 체인 상태 무영향), ADJUSTIE (TIES 미설정 시 no-op) |

### 재현

```
python converted/run_parity.py all     # 822/822 bit-identical
python -m pytest tests/ -q             # 36 passed
```

## 9. 부록 — 왜 SSW가 아니라 upstream tarball인가

다른 11개 패키지는 서버의 SSW 설치본을 썼는데, 이 패키지만 **upstream tarball**( `PoA_pro_current.tar.gz` , sha256 `d205c2b8…` )을 씁니다.

### 이유 ① — SSW 설치본을 읽을 수 없음

서버의 SSW `poa`  설치본은 `idl/stat/` · `idl/fitting/`  디렉터리가 `drwxr--r--` **(실행권한 없음)**이라 `likeli.pro` 를 **읽을 수도 컴파일할 수도 없습니다.** IDL도 같은 사용자로 실행되므로 **오라클 자체가 불가능**합니다. 서버 전체를 검색해도 읽을 수 있는 사본이 없었습니다.

### 이유 ② — 버전 키메라를 피하기 위해

SSW 본체(2004)에 tarball의 `likeli` (이력 2008까지)를 섞으면 **실존한 적 없는 조합**이 됩니다. 특히 **"Apr08:** `/chi2` **출력 보정"**은 **Metropolis를 구동하는 바로 그 반환값**을 바꾸는 변경입니다.

그래서 **일관된 upstream 트리 하나**로 오라클과 변환 대상을 통일했습니다.

### 버전 델타

`solarsoft_dem/pintofale_mcmc/mcmc_dem.pro` (SSW 2004 스냅샷)와 **352줄** 차이가 있습니다. 주요 기능 델타:

`verbose`  키워드 (플로팅 차단)

`sysdev` / `sysgrp`  (emissivity jitter)

`xnsig` , `stordem`

`hitpar(opar)`  버그 수정

burn-in 시그마 추정을 **히스토그램 → 분위수** 방식으로 교체 (VK Apr2012)

**즉 이 인증은 upstream current release에 대한 것**이며, SSW 2004 스냅샷을 쓰는 기존 연구 코드와는 **동작이 다를 수 있습니다.**

### 알고 써야 할 것

**upstream은 시드를 초기화하지 않습니다.** 재현이 필요하면 **반드시 시드를 명시**해야 합니다 — IDL에서도 마찬가지입니다

**결함 2 때문에** `findscale` **의 평활화 스케일이 MCMC 루프에서 무시**되고 상수 3이 쓰입니다

`nburn*nbatch > nsim` 이면 **IDL이 크래시**합니다 (결함 4)

미구현 키워드를 쓰면 **명시적 예외**가 납니다 — 조용히 틀린 답이 나오지 않습니다

bit-parity는 **linux x86_64 + 그 libm**에서만 보장됩니다 (4장)

검증 규모는 **nsim 10~15**의 짧은 체인입니다. 수천~수만 스텝의 실사용 체인에서 끝까지 일치가 유지되는지는 **미검증**입니다

**출처 표시.** 이 문서의 모든 수치·판정·규칙은 `reports/09_parity_certificate.md` , `reports/00_rng_feasibility.md` , `converted/idl_random.py` , `converted/idl_compat.py` , `converted/mcmc_dem.py` 에서 직접 확인했습니다. 1·2장의 설명 구성과 결함 재현 원칙에 대한 평가는 코드와 인증서에 근거한 해석이며, Kashyap & Drake 1998 원논문 대조는 수행하지 않았습니다.

**참고:** `work/CONVERSION_SUMMARY.md` 와 `work/SPECIAL_CASES.md` 는 아직 pintofale_mcmc를 "bit-parity 근본 불가"로 기록하고 있어 이 인증서와 상충합니다.
