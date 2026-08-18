쉬운 한국어 안내
================

이 폴더는 SITES와 Grid-SITES의 IDL 원본 배포본이다.

- SITES는 한 픽셀의 여러 채널 관측값에서 DEM을 반복 계산한다.
- Grid-SITES는 많은 픽셀이나 시간 자료를 더 빠르게 처리하기 위한 방법이다.
- 가장 간단한 시작 예제는 dem_example_sites.pro와 dem_example_gridsites.pro다.
- AIA 예제를 실행하려면 AIA/SDO SSW 환경이 필요하며, 날짜와 폴더 위치는 사용 환경에 맞게
  바꿔야 한다.
- 관측 오차와 온도 응답함수의 오차가 모두 필요하다. 온도 구간 수, 온도 범위, 부드럽게 만드는
  커널 폭은 사용자가 정한다.
- 원저자는 AIA만 쓸 때 약 0.6 MK 아래로 내려가지 말 것을 권하고, 플레어 수준의 높은 온도는
  충분히 시험하지 않았다고 밝힌다.

이 방법을 연구에 사용하면 다음 두 논문을 인용해야 한다.

[1] Morgan & Pickering (2019), SITES:
    https://doi.org/10.1007/s11207-019-1525-4
[2] Pickering & Morgan (2019), Grid-SITES:
    https://doi.org/10.1007/s11207-019-1526-3

아래는 원저자가 배포한 영어 안내문이며 내용은 그대로 보존했다.
----------------------------------------------------------------

This package implements the SITES DEM inversion method described in [1]. The package also includes the Grid-SITES method, that can greatly increase the efficiency of SITES (or other DEM inversions) described in [2].

[1] https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract
[2] https://ui.adsabs.harvard.edu/abs/2019SoPh..294..136P/abstract

If you use these methods in your work, please cite the papers.

QUICK-START
Perhaps the quickest way of using these codes for AIA is to use the example calls provided in dem_example_sites.pro and dem_example_gridsites.pro. Users will need to change some values in these codes such as dates and directory structure. Note that you must be in your AIA/SDO SSW environment for these examples to work.

SITES
The core inversion function is dem_sites.pro. This will invert a single multi-channel measurement (i.e. one pixel). See dem_example_sites.pro for an example of using SITES for AIA data, particularly for looping through multiple pixels to create a DEM map of larger regions.

GRID-SITES
The Grid-SITES function is dem_gridsites.pro. This will invert multi-channel measurements for multiple pixels. See dem_example_gridsites.pro for an example of using Grid-SITES for AIA data. Grid-SITES uses the dem_sites.pro function to invert the DEMs, but is far more efficient, particularly for larger images or time series. 

ERRORS, PREPARATION OF DATA, AND TEMPERATURE RESPONSES
SITES requires the estimated errors of both the measurement and the response functions. The procedure dem_openaiafiles_sites.pro provides the measurement data cube in multiple channels, plus the errors. It can also combine multiple consecutive observations to improve the signal to noise (and other convenient tasks, for example extracting a  region). See the dem_example_sites procedure for an example on using dem_openaiafiles_sites. 

The dem_aiainterpolresponse_sites.pro procedure takes the structure returned by aia_get_response (standard AIA SSW calibration/analysis software) and interpolates the response curves to the user-defined temperature bins. It also estimates the uncertainty in the responses. 

OTHER INSTRUMENTS
In principle, if one has multi-wavelength measurements, appropriate temperature response curves, and measurement/response uncertainties, SITES (and Grid-SITES) should work. Note that SITES has been described exclusively in the context of AIA in [1], and I would be very interested to see results applied to other instruments. See contact below.

INVERSION OPTIONS
There are only 3 user-set parameters that effect the core inversion routine. These are:
1) Number of temperature bins. I usually set this at around 45.
2) Range of temperatures. I usually set a range of around 0.5 to 8MK, on a regular logarithmic temperature scale. I have not tested the code on flare-like temperatures. See [1] for a discussion of low <1MK temperatures in the context of AIA (basically never go below 0.6MK with AIA only).
3) Width of smoothing kernel. The default setting has been found through trial and error - a compromise between stable solutions and over-smoothing the results. Users can provide their own smoothing kernel if they wish to experiment with this parameter.

CONTACT
Huw Morgan, Aberystwyth University, hmorgan@aber.ac.uk
