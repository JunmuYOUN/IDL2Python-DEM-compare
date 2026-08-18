;+
; function simple_reg_dem_aia_wrapper, dn, exptimes, logt, chi2=chi2, tr_struct=tr_struct, tresp=tresp
;
; A simple AIA wrapper for simple_reg_dem.pro. Inputs:
;	dn: A stack of AIA images, dimensions nx by ny by n_channels. For simplicity, assumes the
;		channels are in the same order as in output of aia_get_response,
;		and 304 is omitted. 
;	exptimes: Array of exposure times corresponding to dn. Dimensions n_channels.
;
; Returns:
;	Array of output DEMs, dimensions nx by ny by n_temperatures. See simple_reg_dem or paper
;	for details.
;
; Outputs:
;	logt:	Temperature array of output DEM. Same temperatures as tr_struct, but limited
;			to logtmin...logtmax
;
; Optional keywords:
;	chi2: 	Reduced chi squared array, dimension nx by ny, as returned by simple_reg_dem.
;	tr_struct:	temperature response structure. Should be formatted as returned by
;				aia_get_response(/temp,/dn). Will call aia_get_response with that
;				syntax if it's undefined on input.
;	tresp:	The temperature response array - same as in tr_struct, but with 304 removed	
;			and temperatures limited as in logt.
;	timedepend_date: Date of observations to take into account time dependence of AIA temperature
;			response functions.
;	quiet: Don't print DEM computation time benchmark.
;	logtmin: Minimum of logt for DEM computation (default: 5.5)
;	logtmax: Maximum of logt for DEM computation (default: 7.5)
;
; Author: Joseph Plowman
; 11/21/2017
;-
function simple_reg_dem_aiawrapper, dn, exptimes, logt, chi2=chi2, tr_struct=tr_struct, tresp=tresp, $
		timedepend_date=timedepend_date, quiet=quiet, logtmin=logtmin, logtmax=logtmax

	nx = n_elements(dn[*,0,0])
	ny = n_elements(dn[0,*,0])

	if(n_elements(logtmin) eq 0) then logtmin = 5.5 ; Min and Max (log10) temperature
	if(n_elements(logtmax) eq 0) then logtmax = 7.5

	; Get temperature response, take out 304 channel, and limit temperatures to logtmin...logtmax:
	if(n_elements(tr_struct) eq 0) then tr_struct = aia_get_response(timedepend_date=timedepend_date,/temp,/dn)
	ichan = where(tr_struct.channels ne "A304")
	ilogt = where(tr_struct.logte ge logtmin and tr_struct.logte le logtmax)
	nt0 = n_elements(ilogt)
	nchan = n_elements(ichan)
	channels = tr_struct.channels(ichan)
	tresp = tr_struct.all(ilogt#(1+lonarr(nchan)),(1+lonarr(nt0))#ichan)
	logt=tr_struct.logte[ilogt]

	; Populate the error array for the inversion using the DN:
	errs = dblarr(nx,ny,nchan)
	for i=0,nchan-1 do errs[*,*,i] = AIA_BP_ESTIMATE_ERROR(dn[*,*,i] > 0,replicate(strmid(channels[i],1),nx,ny))

	; Compute DEM:
	tstart = systime(1)
	dems_out = simple_reg_dem(dn, errs, exptimes, logt, tresp, chi2)
	if(keyword_set(quiet) eq 0) then begin
		print,'Computed '+strtrim(string(nx*long(ny)),2)+' DEMs in '+strtrim(string(systime(1)-tstart),2)+' seconds'
	endif

	return, dems_out
	
end