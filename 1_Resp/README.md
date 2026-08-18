# AIA temperature response

`aia_temperature_response.asdf`는 SunPy/aiapy 기반의
[`synthesizAR`](https://github.com/wtbarnes/synthesizAR)에 포함된 SDO/AIA 온도응답 표를
원본 그대로 받은 파일이다. 업스트림 생성 코드는 SolarSoft
`aia_get_response(/temp, /dn)`의 결과를 ASDF로 저장한다.

- 채널: 94, 131, 171, 193, 211, 335 Å
- 원본 온도 격자: `log10(T/K) = 4.0 .. 9.0`, 101점
- 응답 단위: `DN cm^5 / (pixel s)`
- 원본 SHA-256: `aec75a29b2abe108a086932e5dd6a098cba628d69e8fc3ffa2e8ecfff4708740`
- 고정한 upstream commit: `77aab1767e25bc10a200e1b26270da1de20922c2`

`aia_temperature_response.npz`는 같은 값을 NumPy만으로 읽을 수 있게 변환한 파일이다.
`3_DEM_PoC/DEM_solver_PoC.ipynb`는 이 NPZ를 읽고 관측 채널 순서에 맞춰 보간한다.

재생성 명령:

```powershell
conda run -n ssw python 1_Resp/fetch_aia_temperature_response.py
```

현재 표는 upstream 생성 호출에 `/timedepend_date`, `/evenorm`, `/chiantifix`가 없으므로
관측시각별 성능저하, EVE 정규화, 94/131 Å CHIANTI 보정을 적용하지 않은 기준 응답이다.
정량 연구에서는 분석 조건에 맞는 SSW/CHIANTI 설정으로 다시 생성해야 한다.

참고:

- [SunPy GENX reader](https://docs.sunpy.org/en/stable/generated/api/sunpy.io.special.genx.read_genx.html)
- [aiapy response API](https://aiapy.readthedocs.io/en/stable/reference/response.html)
- [upstream ASDF generation example](https://gist.github.com/wtbarnes/4fb5ef1ee2003aa219a533e01a957356)
