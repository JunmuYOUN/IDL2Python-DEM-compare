# IDL2Python-DEM — 고전 DEM solver의 IDL 원본과 Python 변환본

태양 코로나 **DEM(Differential Emission Measure) inversion 코드 12종**을 IDL에서 Python으로
옮긴 결과다. 각 패키지마다 원본 `.pro`와 변환된 `.py`를 나란히 두고, 무엇을 어떻게 옮겼는지를
설명한 문서를 함께 담았다.

변환은 눈으로 대조한 것이 아니다. **IDL 원본을 실제로 실행해 얻은 oracle 값과 Python 결과를
수치 대조(parity)** 해서 검증했다. 12종 전부 ALL PASS, 회귀테스트 84개, `pintofale_mcmc`는
MCMC 체인 822 스텝이 bit-identical하다.

쓸 곳은 두 가지다.

- **검증 기준선** — 새 DEM 방법(신경망 등)의 결과를 고전 inversion과 비교할 때의 레퍼런스
- **물리 구현 참조** — response 적용, χ² 정의, 정규화 항 같은 세부를 IDL 원본까지 거슬러 확인

## 디렉터리 구조

```
<패키지>/
├── README.md  변환 설명서 — 알고리즘·코드 한 줄씩 해설·검증 결과·경계
├── idl/       원본 .pro (변경 없음. 원저자 배포 그대로)
└── python/    변환된 Python
```

패키지별 `README.md`가 이 자산의 핵심이다. 각 solver가 **어떤 방법인지**, 변환된 코드가
**무슨 일을 하는지 줄 단위로**, IDL에서 무엇을 조심해야 했는지, 그리고 **어디까지 검증됐고
어디부터는 아닌지**를 담고 있다. 아래 표의 패키지 이름에서 바로 열 수 있다.

`demreg`만 예외로 `python_ref/`가 하나 더 있다 — 아래 "demreg의 python_ref" 참조.

## 패키지 목록

parity·테스트 수치는 변환 하네스가 발급한 인증서 기준이다(원 출처: `CONVERSION_SUMMARY.md`).
"경계"는 **bit-parity를 보장하지 않는 범위**를 뜻한다.

| 패키지 | 방법 (참고문헌) | parity | 회귀테스트 | 경계 |
|---|---|---|---|---|
| [`simple_reg_dem`](simple_reg_dem/README.md) | Plowman 정규화 DEM (arXiv:2006.06828) | 3세트 ALL PASS (24 pass / 2 waived), 전수 512×832 | 7/7 | float32 커널 waiver 2건 |
| [`demreg`](demreg/README.md) | Hannah & Kontar GSVD 정규화 | ALL PASS (단일 픽셀 + 2×2 맵) | 5/5 | 초기 guess GSVD 경로·gloci 미검증 |
| [`dem_sites`](dem_sites/README.md) | Morgan & Pickering SITES (2019 SoPh 294, 135 / 294, 136) | ALL PASS 19/0 | 4/4 | `gaussian_function` 커널 주입 |
| [`trace_dem`](trace_dem/README.md) | Kankelborg TRACE 선형 DEM (1998) | ALL PASS 13/0 | 4/4 | response 설정(`int_tabulated`+`trace_t_resp`) 주입 |
| [`firdems`](firdems/README.md) | Plowman FIR fast/iterative/regularized DEM | ALL PASS 21/0 | 5/5 | `firdem_iterate`(2000회 반복) 경로의존 |
| [`cheung_sparse_em`](cheung_sparse_em/README.md) | Cheung sparse-EM / basis pursuit | ALL PASS 15/0 ×2세트 | 4/4 | IDL simplex → `scipy.optimize.linprog`(HiGHS) |
| [`aschwanden_aia_teem`](aschwanden_aia_teem/README.md) | Aschwanden 단일 Gaussian T/EM | ALL PASS 17/0 | 3/3 | χ²를 float32로 맞춰야 정확 parity |
| [`xrt_iterative`](xrt_iterative/README.md) | Weber XRT 순방향 + χ² (`xrt_dem_iterative2`) | ALL PASS 12/0 | 3/3 | MPFIT spline-knot 최적화 루프 경로의존 |
| [`eit_dem`](eit_dem/README.md) | Newmark SOHO/EIT | ALL PASS 25/0 | 6/6 | CHIANTI 순방향(`eit_line_map`) 주입 |
| [`chianti_dem`](chianti_dem/README.md) | CHIANTI spline LM (`dem_fit`) | ALL PASS (15 pass / 2 waived) | 3/3 | IDL `SPLINE`(장력) 완전포팅. demreg/MCMC 래퍼는 제외 |
| [`vdem`](vdem/README.md) | Newton 속도 DEM + GCV 정규화 | ALL PASS 23/0 | 4/4 | **온도 DEM이 아니라 속도 DEM**(시선속도 분포). GCV/Brent 최적화 경계 |
| [`pintofale_mcmc`](pintofale_mcmc/README.md) | Kashyap & Drake MCMC DEM | **822/822 bit-identical** (3세트) | 36/36 | LOOPY/SPLINY/ONLYRAT/MIXIE 경로 |

합계 **12종 / 회귀테스트 84개**.

## 변환에서 반복 확인된 것

1. **oracle 입력 주입** — SSW·CHIANTI response나 라이브러리 함수 출력(`aia_get_response`,
   `eit_line_map`, `trace_t_resp`, `gaussian_function` 등)은 재구현하지 않고 IDL oracle이 뽑은
   값을 그대로 주입했다. 따라서 변환·검증 대상은 **DEM 알고리즘 자체**로 한정된다.
2. **float32 충실성이 결정적** — IDL `FLTARR`, `!pi`(IDL에서 float32다), 커널, χ² 잔차의
   catastrophic cancellation 지점에서는 numpy도 float32로 맞춰야 strict parity가 나온다.
3. **경로의존 최적화는 bit-parity 불가** — `firdem_iterate`, MPFIT, GCV/Brent는 optimizer가
   밟는 경로가 갈린다. 대신 solver가 최적화하는 **물리(순방향 모델 + χ²)** 를 검증하고
   최적화 루프는 경계로 문서화했다. 단 MCMC는 예외였다 — IDL 8.6의 기본 RNG가 MT19937이라
   numpy `RandomState`와 스트림이 일치, 체인 전체가 bit-identical하게 재현됐다.
4. **축·연산자 함정** — IDL `#` ↔ numpy `@`/`np.outer`, 원소별 `>`/`<`는 비교가 아니라
   max/min, `shift` ↔ `np.roll(axis=0)`, GSVD/SVD 부호 규약.
5. **기존 포팅 재검증의 가치** — `demreg` 커뮤니티 포팅에서 실제 버그 3건(logt 중점식·
   tresp 외삽·elogt 공식)을 발견해 수정했다.

## Python 코드 사용 시 주의

- 각 `python/` 디렉터리는 **self-contained**다. 패키지 간 import도, 레포 내 다른 모듈에 대한
  의존도 없다. solver 모듈이 같은 디렉터리의 `chk_dump` 등을 import하므로, 해당 디렉터리를
  `sys.path`에 넣고 쓴다.
- 의존성은 대부분 **numpy + scipy**뿐이다(parity 검증 환경: Python 3.10 / numpy 1.26.4 /
  scipy 1.11.4). 예외는 `demreg`로 `astropy`·`tqdm`·`threadpoolctl`을 추가로 쓴다.
  `run_parity_*.py`는 IDL `.sav`를 읽느라 `scipy.io.readsav`를 쓴다.
- `chk_dump.py`는 parity probe dumper다. 환경변수 `CHK_DIR`이 없으면 **no-op**이라 그대로
  둬도 무해하다.
- `run_parity_*.py`는 IDL oracle이 만든 probe 파일을 읽어 비교하는 **검증 드라이버**다.
  probe 데이터(`.sav`/`.npz`)는 용량 때문에 여기 포함하지 않았으므로 **그대로는 실행되지
  않는다.** 어떤 입력을 어떤 방식으로 주입해 비교했는지가 코드에 명시돼 있어, 변환의
  정확한 검증 조건을 읽기 위한 참조로 넣어 뒀다.

### demreg의 `python_ref`

`demreg/python_ref/`는 **커뮤니티가 먼저 만들어 배포하던 Python 포팅**이다. 이번 변환에서
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

학술적으로 사용할 때는 각 방법의 원 논문을 인용한다. 특히 `dem_sites`는 동봉된
`idl/readme.txt`에서 위 두 논문의 인용을 명시적으로 요구하고 있다.

## 포함하지 않은 것

원본·변환 코드만 담고, 검증 과정에서 생긴 산출물과 대용량 상류 배포본은 제외했다.

- `poa_dist/` — PINTofALE 상류 배포본 전체(44 MB, tarball 포함). 변환 대상이 아니라 참조용이었다
- `firdems.tgz` — 같은 내용이 `firdems/idl/firdems/`에 풀려 있다
- probe 데이터(`.sav`/`.npz`), parity 리포트, 계측용 IDL(`staging/`), 실행 로그, `__pycache__`
