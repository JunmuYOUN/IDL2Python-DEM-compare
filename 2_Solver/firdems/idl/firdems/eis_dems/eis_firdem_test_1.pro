;+
; Name: eis_dem_test_1
;
; Performs EIS test inversions for a DEM consisting of a log-normal distribution with central
; temperature set by temp_in.
;
; Inputs:
;	-temp_in: The central temperature of the lognormal
;
; Outputs:
;	-temps_in: Abscissae in log t of the input DEM returned in dem_in
;	-dem_in: The input DEM
;	-temps_out: Abscissae in log t of the output DEM returned in dem_out
;	-dem_out: The output DEM for one of the noise realizations
;	-channels: The names of the lines in data
;	-data: The data values used for the DEM inversions
;	-errs: The errors corresponding to data
;-
pro eis_firdem_test_1, temp_in, temps_in, dem_in, temps_out, dem_out, channels=channels, $
		data=data, errs=errs, sigma_logn=sigma_logn, niter=niter, exptimes=exptimes, $
		etot=etot, other_noise=other_noise

	common aia_dem_iterative_test_1_static_variables, seed

;	channels = ["fe_09_188.497", "fe_10_184.537", "fe_12_195.119", "fe_14_264.79", "fe_15_284.163", "fe_16_262.976", "si_07_275.361", "si_10_258.371"]
;	channels = ["fe_09_188.497", "fe_10_184.537", "fe_12_195.119", "fe_15_284.163", "fe_16_262.976"]
	channels = ["fe_09_188.497", "fe_10_184.537", "fe_12_195.119", "fe_15_284.163", "fe_16_262.976","mg_05_276.579","mg_06_270.394","mg_07_280.737","fe_09_197.862","fe_11_188.216","fe_11_180.401","s_10_264.233","fe_12_192.394","fe_13_202.044","fe_13_203.826","fe_14_270.519","s_13_256.686","ca_14_193.874","ca_15_200.972","ca_16_208.604","ca_17_192.858","si_07_275.361","si_10_258.371","fe_14_264.79"]

	nchan = n_elements(channels)
	if(n_elements(sigma_logn) eq 0) then sigma_logn = 0.1
	if(n_elements(pix_scales) eq 0) then pix_scales = [2.0,1.0]
	if(n_elements(exptimes) eq 0) then exptimes = 30+dblarr(nchan)
	if(n_elements(etot) eq 0) then etot = 5.0e28

	logtmin=5.0
	logtmax=8.0
	tmin = 10.0^logtmin
	tmax = 10.0^logtmax

	nt = 10*(logtmax-logtmin)/sigma_logn
	logt = logtmin+(logtmax-logtmin)*findgen(nt)/(nt-1)
	logt_in = alog10(temp_in)
	temps = 10^logt
	nx = 50
	ny = 50
	
	dem_in = etot*exp(-(logt-logt_in)^2/(2*sigma_logn^2))/sqrt(4*acos(0)*sigma_logn^2)
	
	images = eis_firdem_synthesize_data(logt, dem_in, channels, exptimes, tr_struct=tr_struct, $
		nx=nx, ny=ny, errs=errs, data0=data0, pix_scales=pix_scales, other_noise=other_noise)	
	
	tstart=systime(1)

	eis_firdem_wrapper, channels, images, errs, dem_out, chi2s=chi2s, tr_struct=tr_struct, $
			logtlo=logtmin, logthi=logtmax, niter=niter

	telapsed = systime(1)-tstart

	firdem_get_pixel,dem_out,0,0,t2,dem2
	dem2=dem2>0

	chi2_95 = estimate_quantile(chi2s,0.95)
;	title = string(format='(%"%7.1e%s%4.3g")',telapsed/nx/ny, $
;			" Seconds per DEM; 95% Chi Squared=",chi2_95)
	title = string(format='(%"%s%4.3g")',"95% Chi Squared = ",chi2_95)

	yrange = [-0.0*max(dem_in), 1.3*([max(dem_in)])]
	xrange=alog10([tmin,tmax])
	plot, alog10(t2), dem2, linestyle=1, yrange = yrange, xrange=xrange, ystyle=1, thick=3, $
			xmargin=[14,3], xtitle="Log!D10!N(T) (Kelvin)", ytitle="DEM [cm^-5/Log10(T)]", title=title

	oplot,alog10(temps),dem_in,linestyle=0, thick=3

	legend,["Input","Output"], linestyle=[0,1], thick=3

	nxp = min([nx,3])
	nyp = min([ny,3])
	ii=1
	print,"Plotted chi squareds: "
	plotted_chi2s=dblarr(nxp*nyp)
	for i=0,nxp-1 do begin
		for j=0,nyp-1 do begin
			ii++
			firdem_get_pixel,dem_out,i,j,t2,dem2
			dem2=dem2>0
			oplot,alog10(t2),dem2,linestyle=1,color=128, thick=3
			plotted_chi2s(i*nyp+j) = chi2s(i,j)
		endfor
	endfor
	print,plotted_chi2s
	
;	images2 = dblarr(nchan)
;	tr2 = dblarr(n_elements(t2),nchan)
;	for i=0,nchan-1 do begin
;		tr02 = spl_init(temps,tr(*,i),/double)
;		tr2(*,i) = spl_interp(temps,tr(*,i),tr02,t2)
;		integrand = dem2*tr2(*,i)
;		images2(i) += integrate(alog10(t2),integrand)
;	endfor
	temps_out = t2
	dem_out = dem2
	temps_in = temps
	
	
end
