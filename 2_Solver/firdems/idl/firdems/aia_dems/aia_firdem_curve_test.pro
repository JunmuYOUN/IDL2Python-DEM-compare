;+
; Name: aia_firdem_curve_test
;
; This procedure tests AIA DEM inversion of arbitrary curves provided by the calling
; procedure. The test is repeated 2500 times for benchmarking and to examine variability
; due to shot and read noise. 25 Of these results are plotted, and one is returned to the
; calling procedure.
;
; Inputs:
;	
;	temps_in0:	Temperature grid for the input DEM.
;	dem_in0:	The input DEM curve to be tested.
;
; Outputs:
;
;	temps_out:	Temperature grid for the output DEM.
;	dem_out:	The extracted output DEM for one of the noise realizations generated.
; 
; Optional Arguments:
;
;	niter:		The maximum number of iterations in the DEM inversion (otherwise the default
;				in aia_firdem_wrapper is used).
;	tr_struct:	The temperature response structure returned by aia_get_response.
;	its:		The number of iterations taken for each noise realization.
;	chi2s:		The chi squared for each noise realization.
;	data:		AIA channel intensities for each noise realization.
;	seed:		Optional input seed.
;	name:		String to add to plot title.
;
; Joseph Plowman (plowman@physics.montana.edu) 09-17-12
;-
pro aia_firdem_curve_test, temps_in0, dem_in0, temps_out, dem_out, niter=niter, tr_struct=tr_struct,its=its,chi2s=chi2s, data=data, seed=seed,name=name

	common aia_firdem_curve_test, seed2
	
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
	if(logtmin lt min(temps_in0)) then logtmin = min(temps_in0)
	logtmax = 7.5
	if(logtmax gt max(temps_in0)) then logtmax = max(temps_in0)
	tmin = 10.0^logtmin
	tmax = 10.0^logtmin
	
	logt = logtmin+(logtmax-logtmin)*findgen(nt)/(nt-1)
	t=10^logt

	dem_in = dblarr(nt)
	dem_in02 = spl_init(temps_in0,dem_in0,/double)
	dem_in = spl_interp(temps_in0,dem_in0,dem_in02,logt)

	images = aia_firdem_synthesize_data(logt, dem_in, channels, exptimes, tr_struct=tr_struct, $
			nx=nx, ny=ny, errs=errs, data0=images0)		
	
	sigmas = dblarr(nchan)
	for i=0,nchan-1 do begin
		sigmas(i)=stddev(images(*,*,i))
	endfor
	
	aia_firdem_wrapper, channels, exptimes, images, dem_out, niter=niter, tr_struct=tr_struct, $
			its=its,chi2s=chi2s, telapsed=telapsed

	firdem_get_pixel,dem_out,0,0,t2,dem2,/zeronegs

	qlvls = [0.25,0.5,0.75]
	nqlvls = n_elements(qlvls)

	ymin = min([min(dem2),min(dem_in)])
	ymax = 1.1*max([max(dem2),max(dem_in)])
	yrange = [ymin,ymax]
	chi2_95 = estimate_quantile(chi2s,0.95)

	title = string(format='(%"%7.1e%s%7.1e")',telapsed/nx/ny," Seconds per DEM; 95% Chi Squared=",chi2_95)
	if(n_elements(name) gt 0) then title = name+title
	plot,alog10(t),dem_in,linestyle=0, title=title, xmargin=[14,3], yrange=yrange,xrange=[tmin,tmax],xtitle="Log10(T)",ytitle="DEM [cm^-5/Log10(T)]"

	nxp = min([nx,5])
	nyp = min([ny,5])
	ii=1
	print,"Plotted chi squareds: "
	plotted_chi2s=dblarr(nxp*nyp)
	for i=0,nxp-1 do begin
		for j=0,nyp-1 do begin
			ii++
			firdem_get_pixel,dem_out,i,j,t2,dem2,/zeronegs
			oplot,alog10(t2),dem2,linestyle=1,color=128
			plotted_chi2s(i*nyp+j) = chi2s(i,j)
		endfor
	endfor
	print,plotted_chi2s
	legend,["Input","Output"], linestyle=[0,1]
	oplot,alog10(t),dem_in,linestyle=0	
	
	negs = where(dem2 lt 0,count)
	if (count gt 0) then dem2(negs) = 0
	
	temps_out = t2
	dem_out = dem2
	temps_in = t
	
end
