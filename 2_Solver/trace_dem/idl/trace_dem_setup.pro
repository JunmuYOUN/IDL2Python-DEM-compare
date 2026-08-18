;+
;NAME:
;	TRACE_DEM_SETUP
;PURPOSE:
;	Calculate and invert the temperature response matrix for a 3-component
;	model of the differential emission measure. The resulting 3x3 matrix,
;	K_inv, can be used to produce differential emission measure curves
;	and temperature maps based on TRACE data in all three EUV wavelengths
;	(171, 195, 284). If the dn input parameter is used, then a single
;	DEM curve is calculated based on three EUV pixel values.
;	
;	NOTE-- By default, units of column emission measure are
;	used (see TRACE_T_RESP).
;CALLING SEQUENCE:
;	trace_dem_setup, [index, dn, basis=basis, Tavg=Tavg, em=em, A=A, $
;		K_inv=K_inv, /volume]
;INPUTS (OPTIONAL):
;	dn - 3-component array containing pixel values in DN for 171, 195, 284.
;		The three values may be for a single co-aligned pixel, or averaged
;		over many pixels. Pedestal (approx. must be already subtracted!
;	index - TRACE index structure containing, at minimum, the wavelength
;		(index.wave_len) and exposure time (index.sht_mdur) associated
;		with the three values contained in the dn parameter. If not
;		supplied, then it is assumed that the order of wavelengths is
;		171, 195, 284 and that the pixel values are normalized to
;		1.0 s exposure time.
;OPTIONAL INPUT KEYWORDS:
;	basis - A string containing the name of a function that will supply
;		the three basis elements necessary for constructing DEM curves.
;		The basis functions need not be normalized-- but be warned that
;		the program will normalize them so that they integrate to unity,
;		and all of the results will be expressed in terms of renormalized
;		basis functions. Calculated renormalization coefficients are
;		returned through the A keyword. If the basis keyword is not
;		used, then the trace_sbasis function (below) will supply a default
;		basis. The prototype of any supplied basis function should be
;		similar to that of trace_sbasis (see below).
;	volume - if set, use units of volume emission measure instead of
;		column emission measure. See TRACE_T_RESP.
;OPTIONAL OUTPUT KEYWORDS:
;	A - Contains normalization coefficients for the three basis functions.
;	K_inv - Inverse of the K-matrix, which specifies the linear mapping
;		from the basis elements to the pixel values (normalized to DN/s).
;		The inverse of K may then be used to find the coefficients for
;		the three normalized basis functions, and thereby reconstruct a
;		DEM curve from TRACE pixel values in three EUV wavelengths.
;	em - If the dn input is supplied, em contains coefficients for the three
;		normalized basis functions. Because of automatic renormalization, the
;		elements of em ALWAYS have units of emission measure. The total DEM 
;		function is then:
;
;			dem = em(0)*A(0)*trace_sbasis(0,T) + $
;				 em(1)*A(1)*trace_sbasis(1,T) + $
;				 em(2)*A(2)*trace_sbasis(2,T)
;
;		If the default basis is not used, then trace_sbasis in the above
;		expression would be replaced by a function returning the user's 
;		chosen basis.
;	Tavg - If the dn input is supplied, Tavg contains the resulting emission
;		measure averaged temperature.
;	tau - 3-vector giving characteristic temperatures of the DEM basis functions. 
;MODIFICATION HISTORY:
;	C. Kankelborg, 25-Sept-1998
;
;===== Default basis set ======
;NAME:
;	trace_sbasis
;PURPOSE:
;	Implement three DEM basis functions, indexed 0, 1, 2.
;	This basis has been tested extensively in comparison with
;	CHIANTI emission measure curves. It works well for
;	quiet sun and active regions.
;		The functions are based on sinusoidal bell curves
;	with maxima of 1.0 at log(T) = 5.95, 6.13, and 6.31. These
;	are close to the peaks in the response curves for the
;	TRACE 171, 195, and 284 angstrom passbands. Basis
;	element 0 evaluates to 1.0 for temperatures between
;	log(T) = 5.0 and 5.95, allowing better fits to DEM
;	curves derived from spectroscopic data.
;	NOTE: This basis is NOT normalized to unit total EM.
;CALLING SEQUENCE:
;	result = trace_sbasis(index, T)
;INPUTS:
;	index - denotes which basis function is to be evaluated
;		at temperature T. Valid values are 0, 1, or 2.
;	T - A temperature or array of temperatures at which to
;		evaluate the basis function.
;-

function sigmoid, theta
;a simple, sinusoidal bell curve
pi=3.14159265359
result = (sin(theta))^2
ss = where((theta gt pi) or (theta lt 0.0))
if ss(0) ne -1 then result(ss)=0.0
return,result
end



function trace_sbasis, index, T
logT = alog10(T)
pi = 3.14159265359

logTa = 5.95 ;peak of basis element 0 sigmoid
logTb = 6.31  ;peak of basis element 2 sigmoid
logTc = 5.0  ;low-T cutoff for basis element 0

;locate T0, the peak of the sigmoidal curve
case index of 
	0: logT0 = logTa
	1: logT0 = (logTa+logTb)/2.0
	2: logT0 = logTb
	else: message,'Index should be 0, 1, or 2.'
endcase

omega = 0.5*pi/(logTb-logTa)
theta = omega*(logT-logT0)+0.5*pi

;Compute sigmoidal curve
result = sigmoid(theta)

;Flatline basis element 0 below Ta. Makes it possible to produce
;more realistic DEM functions with this basis.
if index eq 0 then begin
	ss = where(logT lt logTa)
	if ss(0) ne -1 then result(ss)=1.0
	ss = where(logT lt logTc)
	if ss(0) ne -1 then result(ss)=0.0
endif
return, result

end





;============= Main Program =============
pro trace_dem_setup, index, dn, basis=basis, em=em, A=A, Tavg=Tavg, K_inv=K_inv, $
	volume=voume, tau=tau, debug=debug

;If no index passed, assume 171,195,284 order and exposures of 1s.
if n_elements(index) eq 0 then begin
	index = {wave_len:'171', sht_mdur:1.0}
	index = replicate(index,3)
	index(1).wave_len = '195'
	index(2).wave_len = '284'
endif

;If no basis specified, use default
if n_elements(basis) ne 1 then basis = 'trace_sbasis'

;Calculate normalization coefficients and characteristic
;temperatures for the basis elements
T = 10^(5 + findgen(1000)/500)  ;temperature abcissa, 1e5-1e7 Kelvin
K = fltarr(3,3)
A = fltarr(3)
tau = fltarr(3)
for j = 0,2 do begin
	integrand = call_function(basis,j,T)
	A(j) = 1.0/int_tabulated(T, integrand)
	integrand = A(j) * call_function(basis,j,T) * T
	tau(j) = int_tabulated(T, integrand)
endfor

;Calculate K matrix and its inverse
for i = 0, 2 do begin  ;cycle through telescopes
	for j = 0, 2 do begin  ;cycle through basis set
		integrand = trace_t_resp(index(i).wave_len+'ao', T, volume=voume) * A(j) * call_function(basis,j,T)
		K(i, j) = int_tabulated(T, integrand)
	endfor
endfor
K_inv = invert(K)

if n_elements(dn) eq 3 then begin
	;Convert pixel values to DN/s
	p = (dn/index.sht_mdur)

	;Calculate emission measure per basis element and EM averaged temperature
	em = K_inv # p
	Tavg = total(em*tau)/total(em)
endif

end





