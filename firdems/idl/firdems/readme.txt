쉬운 한국어 안내
================

이 폴더는 FIR DEM의 2012년 예비 IDL 배포본이다. AIA와 EIS용 실행 도우미, 일반 기기용
firdem 함수, 결과 표시 함수, 시험 스크립트가 함께 들어 있다.

- AIA: aia_firdem_wrapper를 사용한다. 채널 이름, 노출시간, 채널 영상 배열이 기본 입력이다.
- EIS: 먼저 각 선의 방출률 정보를 파일로 만든 뒤 eis_firdem_wrapper를 사용한다.
- 일반 기기: firdem을 직접 호출한다. 관측값, 오차, 노출시간, 채널 이름, 온도 응답함수가
  필요하다.
- 출력 구조체에는 DEM, DEM 오차, 온도축, 관측값과 다시 계산한 값, 채널별 χ² 등이 들어 있다.
- test=1은 짧은 확인용이고 test=2는 시간이 더 오래 걸리는 상세 확인이다.

관련 논문:
Plowman, Kankelborg & Martens (2013),
"Fast Differential Emission Measure Inversion of Solar Coronal Data"
https://doi.org/10.1088/0004-637X/771/1/2

아래는 원저자가 배포한 영어 안내문이며 내용은 그대로 보존했다.
----------------------------------------------------------------

Preliminary version of J.E. Plowman et al. fast, iterative, regularized DEM code, with test scripts, EIS and AIA wrappers, and routines for displaying the results.

Computing DEMs:

	For AIA, no setup is required. The DEM is computed by calling 'aia_firdem_wrapper'. The minimum inputs are, in order, a list of AIA channels (named 'A171', etc, as in the output of aia_get_response), an array of exposure times, and an array (dimensions [nx,ny,nchannels]) containing the set of AIA channel images to be inverted. The inverted DEM is placed in the fourth argument on return. See the comments on the corresponding file for detailed usage instructions.

	For EIS, a structure containing the emissivities of each spectral line must first be created and written to disk. This is done by the calling 'eis_firdem_add_line_to_dem' for each spectral line that is to be used. The DEM can then be computed for a set of images (one for each spectral line) by calling 'eis_firdem_wrapper'. The minimum inputs are, in order, a list of EIS lines (named for instance, "fe_09_188.497"), an array (dimensions [nx,ny,nlines]) containing the set of EIS line images to be inverted, and an array (dimensions [nx,ny,nlines]) containing estimated errors in each element of the data array. The inverted DEM is placed in the fourth argument on return. See the comments on the corresponding files for detailed usage instructions. 
	
	DEMs for a generic instrument can be computed by calling 'firdem'; the EIS and AIA routines mentioned above are wrappers for this routine. It has the following required arguments:
	
		datain: Array (dimensions [nx,ny,nchan]) containing intensities in each of the spectral
				channels/lines from which the DEM will be computed. Note that negative data
				elements are not zeroed by this routine.
		errsin: The uncertainties corresponding to datain, with the same array dimensions.
		exptimes: Exposure times for each channel.
		data_channeltags: Names of each channel contained in datain.
		tresps: Array (dimensions [n_tresps,n_temps]) containing a set of temperature response
				functions which must include all those found in the input data.
		tresp_logt: Array (dimensions [n_temps]) giving the logt axis for the tresps array.
		tresp_channeltags: Array (dimensions [n_tresps]) labeling each channel contained
							in tresps. Labels must be a superset of data_channeltags.
	
		dem_out: Structure containing the output DEMs for each pixel in datain. Has the
			following elements:
				t(nt):	Array of temperatures at which the basis elements 
					have been evaluated.
				basis(nt,n_bases): The amplitudes of the basis elements, in
					units of differential emission measure (cm^-5 per unit
					alog10(T)). The number of basis elements used is set by
					the optional input 'nb2'.
				coffs(nx,ny,n_bases): Sets of basis element coefficients, which define the
					DEM for each pixel in the input data. May have small residual negative
					elements, which the user may safely zero.
				t0a(n_bases): Array listing the centers of each basis element

	
Using computed DEMs:

	The code returns a structure containing an array of temperatures, an array of the basis elements evaluated at each of those temperatures, and an array containing the basis coefficients which define the DEM corresponding to each pixel in the input image stack. I have included a few routines that operate on this structure to produce more directly interpretable output:
		
		firdem_image: returns an array containing the DEM at a single specified temperature at every pixel.
		
		firdem_get_pixel: returns the single DEM curve corresponding to a specified x and y pixel coordinate in the input image stack.
		
		firdem_emwtemp_image_color: Calculates a true-color image with hue indicating emission-weighted median temperature and brightness indicating total emission measure.

		firdem_vchan_image: ; Makes a true color image from a DEM computed by firdem using 3 'virtual' channels with specified temperature response functions.
		
Test scripts (These create EPS output in the eis_dems/plots and aia_dems/plots directory):

	aia_dems/aia_firdem_test_script_1.pro and eis_dems/eis_firdem_test_script_1.pro: Tests ability to recover a log-normal DEM of width 0.1 and total EM 5.0x10^28 at several input temperatures. Each log-normal distribution is tried repeatedly at several input temperatures, and the results plotted to get an idea of the statistical variance in the recovered DEMs.
		
	aia_dems/aia_firdem_test_script_2.pro and eis_dems/eis_firdem_test_script_2.pro:  computes and plots the response of our fast DEM method to  of log-normal distributions with specified width and total em, at temperatures ranging from 10^5.5 to 10^7.0 Kelvin. The DEMs are plotted in a single 2-D intensity plot, with each vertical slice corresponding to a DEM like that plotted by aia_test_dem_iterative_script.pro. 
	
	aia_dems/aia_firdem_test_script_4.pro: Tests ability to recover a set of test DEMS composed of the sums of 5 random log-normal distributions from AIA data. Each log-normal distribution is tried repeatedly at several input temperatures, and the results plotted to get an idea of the statistical variance in the recovered DEMs. aia_firdem_test_script_3.pro is a simplified version of this script which computes a single such DEM.
	
	eis_dems/eis_firdem_warren_comparison_script.pro: Computes DEMs using the EIS data and uncertainties reported in ApJ, 734, 90.
	
	aia_dems/aia_firdem_chianti_test_script.pro: Tests ability to recover the example quiet sun, active region, and flare DEMs provided with CHIANTI.

These scripts may all be run in IDL from their directory with, for instance, @aia_firdem_test_script_1.pro
	
Please email plowman@physics.montana.edu with any questions -- JEP 10-12-12
