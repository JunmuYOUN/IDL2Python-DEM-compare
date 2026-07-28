;+
; Name: aia_firdem_test_4
;
; Performs AIA test inversions for a DEM consisting of npeaks summed log-normal distributions,
; with uniformly distributed EM, width, and central temperature.
;
; Inputs:
;	npeaks: The number of peaks
;	emrange[2]: The upper and lower limits on the randomly chosen total EM for each lognormal
;	sigrange[2]: The upper and lower limits on the widths of each lognormal
;	trange[2]: The upper and lower limits on the central log temperature of each lognormal
;
; Optional Inputs
;	niter: Maximum number of iterations for the DEM inversion (optional)
;	tr_struct: Temperature response, as returned by aia_get_response (optional)
;	seed:	The random number seed variable to be used. By default, a fixed value is used
;			each time, which causes the same random DEM to be generated each time (The AIA noise
;			realizations are different each time, however).
;
; Optional Outputs:
;	temps_in: Abscissae in log t of the input DEM returned in dem_in
;	dem_in: The input DEM
;	temps_out: Abscissae in log t of the output DEM returned in dem_out
;	dem_out: The output DEM for one of the noise realizations
; 	its: Array containing the number of iterations taken for each DEM inversion performed
;	chi2s: Array containing the final chi-squared for each DEM inversion performed
;	channels: The names of the channels in data, consistent with those in aia_get_response
;	data: The data values used for the DEM inversions
;	errs: The errors corresponding to data
;-
pro aia_firdem_test_4, npeaks, emrange, sigrange, trange, temps_in=temps_in, dem_in=dem_in, $
		temps_out=temps_out, dem_out=dem_out, niter=niter, tr_struct=tr_struct,its=its, $
		chi2s=chi2s, channels=channels, data=data, errs=errs,seed=seed

	if(n_elements(seed) eq 0) then seed=34985739857398

	if(n_elements(tr_struct) eq 0) then	tr_struct = aia_get_response(/temp,/dn)
	ichan = where(tr_struct.channels ne "A304")
	channels = tr_struct.channels(ichan)
	nchan = n_elements(channels)
	tr0 = tr_struct.all(*,ichan)
	logt0=tr_struct.logte
	t0 = 10^logt0

	nx=50
	ny=50
	exptimes = 1.0+dblarr(nchan)

	nt = 481
	tr = dblarr(nt,nchan)

	logtmin = 5.5
	logtmax = 7.5
	tmin = 10.0^logtmin
	tmax = 10.0^logtmin
	
	logt = logtmin+(logtmax-logtmin)*findgen(nt)/(nt-1)
	t=10^logt

	dem_in = dblarr(nt)
	for i=0,npeaks-1 do begin
		em = emrange(0)+(emrange(1)-emrange(0))*randomu(seed)
		sig = sigrange(0)+(sigrange(1)-sigrange(0))*randomu(seed)
		logt_in = trange(0)+(trange(1)-trange(0))*randomu(seed)
		dem_in += em*exp(-(logt-logt_in)^2/(2*(sig)^2))/sqrt(4*acos(0)*(sig)^2)
	endfor

	images = aia_firdem_synthesize_data(logt, dem_in, channels, exptimes, tr_struct=tr_struct, $
			nx=nx, ny=ny, errs=errs, data0=images0)		
	
	sigmas = dblarr(nchan)
	for i=0,nchan-1 do begin
		sigmas(i)=stddev(images(*,*,i))
	endfor
	
	tstart=systime(1)
	
	aia_firdem_wrapper, channels, exptimes, images, dem_out, niter=niter, tr_struct=tr_struct, its=its,chi2s=chi2s
	
	telapsed = systime(1)-tstart

	help,dem_out.basis

	firdem_get_pixel,dem_out,0,0,t2,dem2,/zeronegs

	qlvls = [0.25,0.5,0.75]
	nqlvls = n_elements(qlvls)

	ymin = 0; min([min(dem2),min(dem_in)])
	ymax = 1.3*max([max(dem_in)])
	yrange = [ymin,ymax]

	chi2_95 = estimate_quantile(chi2s,0.95)

	title = string(format='(%"%8.5g%s%5.3g")',1000.0*telapsed/nx/ny," ms per DEM; 95% Chi Squared=",chi2_95)

	plot,alog10(t),dem_in,linestyle=0, title=title, xmargin=[14,3], thick=3, yrange=yrange,xrange=[tmin,tmax],xtitle="Log!D10!N(T) (Kelvin)",ytitle="DEM [cm!E-5!N/Log!D10!N(T)]"
	nxp = min([nx,5])
	nyp = min([ny,5])
	ii=1
	print,"Plotted chi squareds: "
	plotted_chi2s=dblarr(nxp*nyp)
	for i=0,nxp-1 do begin
		for j=0,nyp-1 do begin
			ii++
			firdem_get_pixel,dem_out,i,j,t2,dem2,/zeronegs
			oplot,alog10(t2),dem2,linestyle=1,color=128, thick=3
			plotted_chi2s(i*nyp+j) = chi2s(i,j)
		endfor
	endfor
	print,plotted_chi2s
	legend,["Input","Output"], linestyle=[0,1],charsize=0.7, /right, thick=3
	oplot,alog10(t),dem_in,linestyle=0	
	
	negs = where(dem2 lt 0,count)
	if (count gt 0) then dem2(negs) = 0
	
	temps_out = t2
	dem_out = dem2
	temps_in = t
	
	
end
