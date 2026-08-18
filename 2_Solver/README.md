# IDL2Python-DEM — DEM 코드 12종의 IDL 원본과 Python 변환본

태양 코로나의 DEM을 계산하는 기존 IDL 코드 12종을 Python으로 옮긴 저장소다. 각 폴더에는
원본 `.pro`, 변환한 `.py`, 그리고 코드의 목적과 한계를 설명하는 README가 들어 있다.

Python 코드가 제대로 옮겨졌는지는 다음과 같이 확인했다. 먼저 IDL 원본을 실행해 기준값을
만들고, 같은 입력을 Python 코드에 넣어 두 결과를 수치로 비교했다. 정해 둔 비교 항목은 12종
모두 통과했고, 회귀 테스트 84개도 통과했다. `pintofale_mcmc`는 짧은 시험용 MCMC 체인
822단계가 비트 단위까지 같았다.

여기서 `ALL PASS`는 **문서에 적힌 입력과 비교 항목을 모두 통과했다**는 뜻이다. 모든 관측자료,
모든 선택 기능, 모든 컴퓨터 환경에서 결과가 같다는 뜻은 아니다. 각 README의 “이번에 확인한
범위”와 “한계”를 함께 읽어야 한다.

이 저장소는 주로 다음 두 가지 용도로 쓸 수 있다.

- 새 DEM 방법의 결과를 기존 방법과 비교할 때 기준으로 사용
- 응답함수 적용, χ² 계산, 정규화 방식 등을 IDL 원본과 Python 코드에서 함께 확인

## 이 저장소에서 쓰는 검증 용어

- **결과 일치 검증**: 같은 입력을 IDL과 Python에 넣고 결과가 같은지 확인하는 작업
- **IDL 기준값**: IDL 원본을 실행해서 미리 저장한 비교용 결과
- **허용한 차이**: 원인이 확인됐고 최종 결과에는 문제가 없어서 예외로 기록한 작은 차이
- **이번에 확인한 범위**: 실제로 실행하고 비교한 함수, 입력, 선택 기능의 범위

`ill-posed`, inversion, regularization, response, forward model, Jacobian처럼 분야에서 널리 쓰는
용어는 억지로 풀어 쓰지 않고 그대로 사용한다.

## 디렉터리 구조

```
<패키지>/
├── README.md  쉬운 방법 설명, 코드 해설, 확인 결과와 한계
├── idl/       원본 .pro (변경 없음. 원저자 배포 그대로)
└── python/    변환된 Python
```

패키지별 `README.md`에는 각 코드가 **무엇을 계산하는지**, **어떤 가정을 쓰는지**, IDL에서
Python으로 옮길 때 무엇을 조심했는지, 그리고 **어디까지 확인했는지**를 적었다. 아래 표의
패키지 이름을 누르면 바로 열 수 있다.

`demreg`만 예외로 `python_ref/`가 하나 더 있다 — 아래 "demreg의 python_ref" 참조.

## 패키지 목록

아래 수치는 변환 당시 저장한 검증 기록을 기준으로 한다. 마지막 열에는 아직 실행해 보지 않았거나
이번 Python 변환에 넣지 않은 부분을 적었다.

| 패키지 | 방법 | IDL과 Python 결과 대조 | 회귀 테스트 | 이번에 확인하지 않은 부분 |
|---|---|---|---|---|
| [`simple_reg_dem`](simple_reg_dem/README.md) | Plowman 정규화 DEM | 3세트 통과(24개 일치, 차이 2개 허용), 512×832 자료 확인 | 7/7 | float32 커널의 작은 차이 2개 |
| [`demreg`](demreg/README.md) | Hannah & Kontar GSVD 정규화 | 단일 픽셀과 2×2 맵 통과 | 5/5 | 초기 DEM을 자동으로 만드는 GSVD 경로와 `gloci` |
| [`dem_sites`](dem_sites/README.md) | SITES 반복 DEM | 19개 항목 통과 | 4/4 | Gaussian 커널은 외부 계산값 사용 |
| [`trace_dem`](trace_dem/README.md) | TRACE 3채널 선형 DEM | 13개 항목 통과 | 4/4 | TRACE 응답함수 생성 과정 |
| [`firdems`](firdems/README.md) | FIR DEM | 21개 항목 통과 | 5/5 | 최종 반복 단계인 `firdem_iterate` |
| [`cheung_sparse_em`](cheung_sparse_em/README.md) | Sparse-EM | 2세트에서 각각 15개 항목 통과 | 4/4 | IDL의 `simplex` 대신 SciPy `linprog` 사용 |
| [`aschwanden_aia_teem`](aschwanden_aia_teem/README.md) | 단일 Gaussian T/EM 적합 | 17개 항목 통과 | 3/3 | 실관측 대형 맵과 태양 가장자리 마스크 |
| [`xrt_iterative`](xrt_iterative/README.md) | XRT 관측값·χ² 계산 | 12개 항목 통과 | 3/3 | 전체 MPFIT 반복 계산은 옮기지 않음 |
| [`eit_dem`](eit_dem/README.md) | SOHO/EIT 기준 DEM 보정 | 25개 항목 통과 | 6/6 | `eit_line_map`을 이용한 관측값 계산 |
| [`chianti_dem`](chianti_dem/README.md) | CHIANTI 스플라인 LM 적합 | 15개 일치, 중간값 차이 2개 허용 | 3/3 | `dem_fit` 밖의 보조 래퍼 |
| [`vdem`](vdem/README.md) | 속도 DEM 정규화 | 23개 항목 통과 | 4/4 | GCV/Brent로 가중치를 자동 선택하는 부분 |
| [`pintofale_mcmc`](pintofale_mcmc/README.md) | PINTofALE MCMC DEM | 3세트, 822단계가 비트 단위까지 일치 | 36/36 | LOOPY, SPLINY, ONLYRAT, MIXIE 선택 기능 |

합계 **12종 / 회귀테스트 84개**.

## 변환하면서 여러 번 확인한 점

1. **응답함수까지 모두 옮긴 것은 아니다.** SSW·CHIANTI 함수의 출력값은 IDL에서 먼저 계산해
   Python에 입력했다. 따라서 주로 확인한 대상은 DEM 계산 부분이다.
2. **숫자 자료형이 중요했다.** IDL의 `FLTARR`와 `!pi`는 float32다. Python에서 무조건 더 높은
   정밀도를 쓰면 반올림 순서가 달라져 오히려 IDL 결과와 어긋날 수 있었다.
3. **반복 최적화는 중간 경로가 달라질 수 있다.** `firdem_iterate`, MPFIT, GCV/Brent는 작은
   반올림 차이 때문에 다음 반복 단계가 달라질 수 있다. 그래서 옮긴 범위와 그렇지 않은 범위를
   각 문서에 따로 적었다. MCMC 시험에서는 IDL과 NumPy의 난수 흐름을 맞춰 짧은 체인 전체가
   비트 단위까지 같았다.
4. **IDL 축·연산자 규칙** — IDL `#` ↔ NumPy `@`/`np.outer`, 원소별 `>`/`<`는 비교가 아니라
   max/min, `shift` ↔ `np.roll(axis=0)`, GSVD/SVD 부호 규약.
5. **기존 포팅 재검증의 가치** — `demreg` 커뮤니티 포팅에서 실제 버그 3건(logt 중점식·
   tresp 외삽·elogt 공식)을 발견해 수정했다.

## Python 코드 사용 시 주의

- 각 `python/` 디렉터리는 **다른 패키지와 분리되어 있다.** solver 모듈이 같은 디렉터리의
  `chk_dump` 등을 불러오므로, 해당 디렉터리를
  `sys.path`에 넣고 쓴다.
- 의존성은 대부분 **numpy + scipy**뿐이다(결과 일치 검증 환경: Python 3.10 / numpy 1.26.4 /
  scipy 1.11.4). 예외는 `demreg`로 `astropy`·`tqdm`·`threadpoolctl`을 추가로 쓴다.
  `run_parity_*.py`는 IDL `.sav`를 읽느라 `scipy.io.readsav`를 쓴다.
- `chk_dump.py`는 중간 계산값을 저장하는 검증용 파일이다. 환경변수 `CHK_DIR`이 없으면 아무
  일도 하지 않으므로 그대로 두어도 된다.
- `run_parity_*.py`는 IDL 기준값과 Python 결과를 비교하는 **검증용 실행 파일**이다.
  비교 데이터(`.sav`/`.npz`)는 용량 때문에 여기 포함하지 않았으므로 **현재 상태로는 실행되지
  않는다.** 어떤 입력을 어떤 방식으로 넣어 비교했는지가 코드에 적혀 있어, 변환의
  정확한 검증 조건을 읽기 위한 참조로 넣어 뒀다.

### demreg의 `python_ref`

`demreg/python_ref/`는 **기존에 공개되어 있던 Python 변환본**이다. 이번 작업에서
새로 만든 것이 아니라 **재검증 대상**이었고, 그 과정에서 버그 3건을 찾았다.
실제로 쓸 것은 `demreg/python/` 쪽이며, `python_ref/`는 비교 이력을 남기려고 보존한다.

## 원본 코드의 출처와 인용

`idl/` 아래 `.pro`는 **전부 제3자 코드**이며 원저자 배포본 그대로다. 우리가 작성한 것이 아니다.

| 패키지 | 원본 배포처 |
|---|---|
| `firdems` | `solar.physics.montana.edu/plowman/firdems.tgz` |
| `dem_sites` | SolarSoft `dem_sites/idl/` |
| `pintofale_mcmc` | PINTofALE 배포본 |
| `chianti_dem` | CHIANTI |
| 나머지 | SolarSoft (<https://sohoftp.nascom.nasa.gov/solarsoft/>) |

**라이선스 주의** — 원본 파일에는 LICENSE 파일도, 소스 내 copyright 표기도 없다(SolarSoft
코드의 일반적인 상태다). 명시적 재배포 허가가 문서화돼 있지 않다는 뜻이므로, 이 디렉터리는
연구용 참조로 다루고, 외부 공개·재배포 범위를 넓힐 때는 각 원저자에게 확인이 필요하다.

학술적으로 사용할 때는 각 방법의 원 논문을 인용한다. 특히 `dem_sites`의 원본
`idl/readme.txt`는 SITES와 Grid-SITES 논문 인용을 요청한다. 확인한 논문은 각 패키지 README에
정리했으며, 대표 자료는 다음과 같다.

- Plowman & Caspi (2020), [A Fast, Simple, Robust Algorithm for Coronal Temperature Reconstruction](https://arxiv.org/abs/2006.06828)
- Hannah & Kontar (2012), [Differential emission measures from the regularized inversion of Hinode and SDO data](https://doi.org/10.1051/0004-6361/201117576)
- Morgan & Pickering (2019), [SITES](https://doi.org/10.1007/s11207-019-1525-4), Pickering & Morgan (2019), [Grid-SITES](https://doi.org/10.1007/s11207-019-1526-3)
- Plowman, Kankelborg & Martens (2013), [Fast Differential Emission Measure Inversion of Solar Coronal Data](https://doi.org/10.1088/0004-637X/771/1/2)
- Cheung et al. (2015), [Thermal Diagnostics with AIA: A Validated Method for DEM Inversions](https://doi.org/10.1088/0004-637X/807/2/143)
- Aschwanden et al. (2013), [Automated Temperature and Emission Measure Analysis of Coronal Loops and Active Regions Observed with SDO/AIA](https://doi.org/10.1007/s11207-011-9876-5)
- Kashyap & Drake (1998), [Markov-Chain Monte Carlo Reconstruction of Emission Measure Distributions](https://ui.adsabs.harvard.edu/abs/1998ApJ...503..450K/abstract)
- Newton, Emslie & Mariska (1995), [The Velocity Differential Emission Measure](https://ui.adsabs.harvard.edu/abs/1995ApJ...447..915N/abstract)

## 포함하지 않은 것

원본·변환 코드만 담고, 검증 과정에서 생긴 파일과 용량이 큰 원 배포본은 제외했다.

- `poa_dist/` — PINTofALE 원 배포본 전체(44 MB, 압축 파일 포함). 변환 대상이 아니라 참고용이었다
- `firdems.tgz` — 같은 내용이 `firdems/idl/firdems/`에 풀려 있다
- 비교 데이터(`.sav`/`.npz`), 결과 일치 검증 보고서, 중간값 저장용 IDL, 실행 로그, `__pycache__`
