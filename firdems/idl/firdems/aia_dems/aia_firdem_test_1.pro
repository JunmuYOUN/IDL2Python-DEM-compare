;+
; Name: firdem_test_1
;
; Performs AIA test inversions for a DEM consisting of a log-normal distribution with central
; temperature set by temp_in.
;
; Inputs:
;	temp_in: The central temperature of the lognormal
;	niter: Maximum number of iterations for the DEM inversion (optional)
;	tr_struct: Temperature response, as returned by aia_get_response (optional)
;
; Outputs:
;	temps_in: Abscissae in log t of the input DEM returned in dem_in
;	dem_in: The input DEM
;	temps_out: Abscissae in log t of the output DEM returned in dem_out
;	dem_out: The output DEM for one of the noise realizations
; 	its: Array containing the number of iterations taken for each DEM inversion performed
;	chi2s: Array containing the final chi-squared for each DEM inversion performed
;	channels: The names of the channels in data, consistent with those in aia_get_response
;	data: The data values used for the DEM inversions
;	errs: The errors corresponding to data
;
; Joseph Plowman (plowman@physics.montana.edu) 09-17-12
;-
pro aia_firdem_test_1, temp_in, temps_in, dem_in, temps_out, dem_out, niter=niter, $
		tr_struct=tr_struct, its=its, chi2s=chi2s, channels=channels, data=data, errs=errs

	common aia_dem_iterative_test_1_static_variables, seed

	if(n_elements(tr_struct) eq 0) then	tr_struct = aia_get_response(/temp,/dn,/evenorm)
	ichan = where(tr_struct.channels ne "A304")
	channels = tr_struct.channels(ichan)
	tr0 = tr_struct.all(*,ichan)
	logt0=tr_struct.logte
	t0 = 10^logt0

	nx=100
	ny=100

	nt = 481
	nchan = n_elements(channels)
	tr = dblarr(nt,nchan)

	exptimes = 1.0+dblarr(nchan)


	logtmin = 5.5
	logtmax = 7.5
	tmin = 10.0^logtmin
	tmax = 10.0^logtmin
	
	logt = logtmin+(logtmax-logtmin)*findgen(nt)/(nt-1)
	t=10^logt
	logt_in = alog10(temp_in)

	dem_in = dblarr(nt)
	dem_in = 5.0e28*exp(-(logt-logt_in)^2/(2*(.1)^2))/sqrt(4*acos(0)*(.1)^2)
		
	images = aia_firdem_synthesize_data(logt, dem_in, channels, exptimes, tr_struct=tr_struct, $
			nx=nx, ny=ny, errs=errs, data0=images0)

	sigmas = dblarr(nchan)
	for i=0,nchan-1 do begin
		sigmas(i)=stddev(images(*,*,i))
	endfor	

	tstart=systime(1)
	
	aia_firdem_wrapper, channels, exptimes, images, dem_out, niter=niter, tr_struct=tr_struct, its=its,chi2s=chi2s,logtlo=logtmin,logthi=logtmax
	
	telapsed = systime(1)-tstart

	help,dem_out.basis

	firdem_get_pixel,dem_out,0,0,t2,dem2
	chi2_95 = estimate_quantile(chi2s,0.95)
	
;	title = string(format='(%"%7.1e%s%4.3g")',telapsed/nx/ny," Seconds per DEM; 95% Chi Squared=",chi2_95)

	title = string(format='(%"%s%4.3g")',"95% Chi Squared = ",chi2_95)

;	yrange = 1.3*[min([min(dem2),min(dem_in)]), max([max(dem2),max(dem_in)])]
	yrange = [-0.0*max(dem_in), 1.3*([max(dem_in)])]
	xrange=alog10([tmin,tmax])

	plot, alog10(t), dem_in, linestyle=0, yrange = yrange, xrange=xrange, ystyle=1, thick=3, $
			xmargin=[14,3], xtitle="Log!D10!N(T) (Kelvin)", ytitle="DEM [cm^-5/Log10(T)]", title=title
			
	nxp = min([nx,5])
	nyp = min([ny,5])
	ii=1
	print,"Plotted chi squareds: "
	plotted_chi2s=dblarr(nxp*nyp)
	for i=0,nxp-1 do begin
		for j=0,nyp-1 do begin
			ii++
			firdem_get_pixel,dem_out,i,j,t2,dem2
			dem2=dem2>0
			oplot,alog10(t2),dem2,linestyle=1,color=128,thick=3
			plotted_chi2s(i*nyp+j) = chi2s(i,j)
		endfor
	endfor
	print,plotted_chi2s
	legend,["Input","Outputs"], linestyle=[0,1],thick=3
	
	negs = where(dem2 lt 0,count)
	if (count gt 0) then dem2(negs) = 0
	
	images2 = dblarr(nchan)
	tr2 = dblarr(n_elements(t2),nchan)
	for i=0,nchan-1 do begin
		tr02 = spl_init(t,tr(*,i),/double)
		tr2(*,i) = spl_interp(t,tr(*,i),tr02,t2)
		integrand = dem2*tr2(*,i)
	endfor

	temps_out = t2
	dem_out = dem2
	temps_in = t
	
	
end
