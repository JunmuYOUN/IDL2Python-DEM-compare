;+
;NAME:
;	TRACE_DEM
;PURPOSE:
;	Calculate differential emission measure and EM weighted average
;	temperature for a set of TRACE images in all three EUV wavelengths
;	(171, 195, 284). By default, units of column emission measure are
;	used (see TRACE_T_RESP).
;CALLING SEQUENCE:
;	trace_dem, index, data[, basis=basis, Tavg=Tavg, emtot=emtot, emc=emc, $
;		A=A, K_inv=K_inv, /volume]
;INPUTS:
;	data - m x n x 3 array containing pixel values in DN for the 171, 195, 
;		and 284 angstrom EUV passbands of TRACE
;		(but not necessarily in that order). Readout pedestal (approx.
;		83 DN for full resolution images, 95 DN for 2x2 summed images)
;		must be subtracted prior to sending the data to TRACE_DEM.
;	index - TRACE index structure containing, at minimum, the wavelength
;		(index.wave_len) and exposure time (index.sht_mdur) associated
;		with the three images contained in the data parameter. If not
;		supplied, then it is assumed that the order of wavelengths is
;		171, 195, 284 and that the pixel values are normalized to
;		1.0 s exposure time (i.e., units of DN/s).
;OPTIONAL INPUT KEYWORDS:
;	basis - A string containing the name of a function that will supply
;		the three basis elements necessary for constructing DEM curves.
;		See TRACE_DEM_SETUP for more information.
;		-----------------------------------------------------------
;		NOTE: The default basis is NOT normalized to unit total EM.
;			Use keyword A to recover the DEM curve from the 
;			original basis functions. (example below under "em")
;		-----------------------------------------------------------
;	volume - if set, use units of volume emission measure instead of
;		column emission measure. See TRACE_T_RESP.
;OPTIONAL OUTPUT KEYWORDS:
;	A - Contains renormalization coefficients for the three basis functions.
;	K_inv - Inverse of the K-matrix, which specifies the linear mapping
;		from the basis elements to the pixel values (normalized to DN/s).
;		The inverse of K may then be used to find the coefficients for
;		the three normalized basis functions, and thereby reconstruct a
;		DEM curve from TRACE pixel values in three EUV wavelengths.
;	emc - contains emission measure coefficients for the three
;		normalized basis functions. Because of automatic renormalization,
;		the elements of em ALWAYS have units of emission measure.
;			In effect, this m x n x 3 data cube represents a 3-parameter 
;		differential emission measure curve at every pixel. The total DEM
;		function is then:
;
;			dem = emc(0)*A(0)*sbasis(0,T) + $
;				 emc(1)*A(1)*sbasis(1,T) + $
;				 emc(2)*A(2)*sbasis(2,T)
;
;		If the default basis is not used, then sbasis in the above
;		expression would be replaced by a function returning the user's 
;		chosen basis.
;	Tavg - emission measure averaged temperature, i.e. an m x n temperature map.
;	tau - 3-vector giving characteristic temperatures of the DEM basis functions. 
;	emtot - Map of total emission measure along the line of sight, 
;		integrated over the range of temperatures where the basis elements
;		are nonzero. This is just the sum of the three components of em:
;		emtot = total(emc, 3).
;MODIFICATION HISTORY:
;	C. Kankelborg, 25-Sept-1998
;-
pro trace_dem, index, data, basis=basis, emc=emc, A=A, Tavg=Tavg, K_inv=K_inv, $
	emtot=emtot, volume=voume, tau=tau, debug=debug

datasize = size(data)
if datasize(0) ne 3 or datasize(3) ne 3 then message, $
	'First parameter must be a data cube containing three 2D images.'
nx = datasize(1)
ny = datasize(2)
npix = nx*ny

;First, calculate K_inv and A for the chosen (or default) basis.
trace_dem_setup, index,[0.0,0.0,0.0], basis=basis, A=A, K_inv=K_inv, volume=voume, tau=tau

;Convert pixel values to DN/s
p = fltarr(nx,ny,3)
for i=0, 2 do p(*,*,i) = (data(*,*,i)/index(i).sht_mdur)

;Calculate emission measure per basis element and EM averaged temperature
emc = fltarr(nx,ny,3)
Tavg = fltarr(nx,ny)
for x = 0L, nx-1 do begin
for y = 0L, ny-1 do begin
	emc(x,y,*) = K_inv # reform(p(x,y,*),3)
	Tavg(x,y) = total(emc(x,y,*)*tau)/total(emc(x,y,*))
endfor
endfor

emtot = total(emc, 3)

end




