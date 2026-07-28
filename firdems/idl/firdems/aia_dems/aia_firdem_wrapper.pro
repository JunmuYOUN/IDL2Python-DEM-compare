;+
; Name: aia_firdem_wrapper
;
; This is a wrapper for computing differential emission measures (DEMs) from
; AIA data using the firdem code. It operates on a set of AIA images, one for each channel.
;
; Inputs:
;
;	channels: 	The names of the channels corresponding to the set of images 
;				in datain. Must match the channel names returned by 
;				aia_get_response.
;
;	exptimes:	The exposure times for each image.
;
;	datain:		An 3-index [nx,ny,nchannels] array containing a set of AIA 
;				images. The DEM will be computed for each pixel in these 
;				images. Units should be DN.
;
; Outputs:
;
;	dem_out:	Structure containing the results of the DEM computation. Has
;				the following elements:
;					t(nt):	Array of temperatures at which the basis elements 
;							have been evaluated.
;					basis(nt,nb2): The amplitudes of the basis elements, in
;						units of differential emission measure (cm^-5 per unit
;						alog10(T))
;					coffs(nx,ny.nb2): The coefficients of each basis element,
;						which define the DEM. May have small residual negative
;						elements, which the user may safely zero.
;					t0a(n_bases): Array listing the centers of each basis element
;
; Optional Arguments:
;
;	tr_struct:	Temperature response structure, as returned by 
;				aia_get_response(/temp,/all,/dn). Will be computed
;				automatically if undefined. repeated calls with, e.g.,
;				tr_struct=tr_struct will only call aia_get_response
;				once, which saves some time.
;
;	niter:		The maximum number of iterations to take. Default is 2000.
;
;	its(nx,ny):		Optional output - The number of iterations taken at
;					each pixel location.
;
;	chi2s(nx,ny):	Optional output - The chi squared of the DEM (with any
;					residual negative coefficients zeroed) at each pixel
;					location.
;
;	logthi,logtlo:	Minimum and maximum alog10(T) (defaults: 7.5 and 5.5, set in
;					fastdem_firstpass)
;
;	errors: Array with same dimensions as datain which specifies the errors in the data.
;			If not supplied, shot noise proportional to the input data and read noise will
;			be assumed.
;
;	minconds: Minimum condition numbers for SVD used in iterative removel of negative EM. 
;			There should be one for each element of chi2_reg_probs, set below. Generally
;			speaking, a stronger regularization works better accompanied by a more stringent
;			minimum condition. Default: [0.01,0.05,0.1] for the default chi_2_reg_probs of
;			[0.9,0.5,0.01]. Note that if chi2_reg_probs is set but minconds is unassigned, it
;			will instead be set to a constant 0.02 for all input regularization strengths
;
;	chi2_reg_probs: The set of regularization probabilities to try, related to corresponding
;					chi squareds by calling chisqr_cvf(chi2_reg_probs(i),nchan). Trying more
;					probabilities will take longer. The default is [0.9,0.5,0.01], which tries
;					regularization strengths corresponding to typical errors, then one
;					corresponding to relatively large errors.
;
;	chi2_target: Target chi squared for output DEMs. If chisqr_cvf(chi2_reg_probs(0),nchan) is
;				greater than chi2_target, it will be changed so that it is smaller than
;				chi2_target. This may lead to the poorer performance if the resulting
;				regularization strength is too weak (small regularization chi squared). In other
;				words, don't try to overfit the data by asking for smaller chi squared than is
;				consistent with the majority of data values (i.e., reduced chi squared of
;				approximately one or less). Default: chisqr_cvf(0.05,nchan)/nchan.
;
;	telapsed:	The total elapsed time for computation of the DEMs.
;
;	--Joseph Plowman (plowman@physics.montana.edu) 09/17/12
;-

; Returns an array containing the conversion from AIA photons to DN for an input list of
; by channels, which are assumed to be an array of string labels like those returned by
; aia_get_response. Taken from Solar Physics 275, 41. 
function aia_firdem_dnperphot, channels

	refchannels=['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
	refg = [2.128,1.523,1.168,1.024,0.946,0.658,0.596]
	nchan=n_elements(channels)
	g=dblarr(nchan)
	for i=0,nchan-1 do begin
		g(i) = refg(where(refchannels eq channels(i)))
	endfor
	
	return,g
	
end

;+
; Name: aia_firdem_readnoise
;
; Returns an array containing standard deviation of the AIA read noise for an input list of
; channels, which are assumed to be an array of string labels like those returned by
; aia_get_response. Taken from Solar Physics 275, 41. 
;
;	--Joseph Plowman (plowman@physics.montana.edu) 09/17/12
;-
function aia_firdem_readnoise, channels

	refchannels=['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
	refn = [1.14,1.18,1.15,1.2,1.2,1.14,1.18]
	nchan=n_elements(channels)
	n=dblarr(nchan)
	for i=0,nchan-1 do begin
		n(i) = refn(where(refchannels eq channels(i)))
	endfor
	
	return,n
	
end

;+
; Name: aia_firdem_estimate_error
;
; Estimates the AIA error in a data array (units of DN), assuming shot and read noise specified
; by Solar Physics 275, 41. 
;
; Inputs:
;	
;	channels:	The channels labels for the data, which are assumed to be like those returned by
;				aia_get_response.
;	data:		The data array (dimensions [nx,ny,n_elements(channels)]) whose error is to be
;				estimated.
;
;	--Joseph Plowman (plowman@physics.montana.edu) 09/17/12
;-
function aia_firdem_estimate_error, channels, data

	nchan = n_elements(channels)
	dnpp = aia_firdem_dnperphot(channels) ; Photon to DN conversion.
	rdn = aia_firdem_readnoise(channels) ; Read noise.
	sigmas=data
	for i=0,nchan-1 do begin
		sigmas(*,*,i) = sqrt((data(*,*,i)*dnpp(i)) > 0.0*rdn(i) + 0.0 + rdn(i)^2)
	endfor
	
	return,sigmas
	
end

;+
; Name: aia_firdem_synthesize_data
;
; Returns an array of synthetic AIA data from an input DEM, with randomly generated read and
; shot noise. Array will have dimensions [nx,ny,n_elements(channels)]
;
; Inputs:
;	
;	logt:		The log_10(T) axis for the DEM.
;	dem_in:		The input DEM, assumed to have units 1/cm^5/log_10(T).
;	channels:	The channels to use (Note that the 304 Angstrom channel is not generally
;				optically thin, so its emission is not well represented by a simple
;				temperature response function and it is not considered appropriate for
;				DEM reconstruction). Should be a string array like that returned by $
;				aia_get_response.
;	exptimes:	Array of exposure times for each channel.
;
; Optional Inputs:
;	
;	tr_struct:	A temperature response structure computed by aia_get_response. Default:
;				aia_get_response(/temp,/dn,/evenorm).
;	nx:			x dimension of pixel array to generate. Default: 100
;	ny:			y dimension of pixel array to generate. Default: 100
;
; Optional Outputs:
;	
;	errs:		An array containing estimated errors for each element of the output data
;				array.
;	data0:		Same as the output data array, but without errors.
;
;	--Joseph Plowman (plowman@physics.montana.edu) 09/17/12
;-
function aia_firdem_synthesize_data, logt, dem_in, channels, exptimes, tr_struct=tr_struct, $
		nx=nx, ny=ny, errs=errs, data0=data0

	common aia_firdem_synthesize_data_static_variables, seed

	nt = n_elements(logt)
	nchan = n_elements(channels)
	if(n_elements(tr_struct) eq 0) then	tr_struct = aia_get_response(/temp,/dn,/evenorm)
	ichan = lonarr(nchan)
	for i=0,nchan-1 do begin
		ichan(i) = where(tr_struct.channels eq channels(i))
		if(ichan(i) eq -1) then message,"Error: Channel " + channels(i) + $
				" Not found in temperature response structure."
	endfor
	tr0 = tr_struct.all(*,ichan)
	logt0=tr_struct.logte
	t0 = 10^logt0

	if(n_elements(nx) eq 0) then nx=100
	if(n_elements(ny) eq 0) then ny=100

	tr = dblarr(nt,nchan)
	logtmin = min(logt)
	logtmax = max(logt)
	tmin = 10.0^logtmin
	tmax = 10.0^logtmin
	t=10^logt

	data0 = dblarr(nx,ny,nchan)
	data = dblarr(nx,ny,nchan)
	errs = dblarr(nx,ny,nchan)
	
	for i=0,nchan-1 do begin
		tr02 = spl_init(logt0,tr0(*,i),/double)
		tr(*,i) = spl_interp(logt0,tr0(*,i),tr02,logt)
		integrand = dem_in*tr(*,i)
		data0(*,*,i) += int_tabulated(logt,integrand)*exptimes(i)
	endfor
	
	G=aia_firdem_dnperphot(channels)
	print,"Data (Noise-Free)", transpose(data0(0,0,*))
	for j=0,nx-1 do begin
		for k=0,ny-1 do begin
			for i=0,nchan-1 do begin
				G = aia_firdem_dnperphot(channels(i))
				errs(j,k,i)=sqrt(data0(j,k,i)/G + aia_firdem_readnoise(channels(i)))
				photons = randomu(seed,poisson=data0(j,k,i)/G)
				data(j,k,i) = photons*G+randomn(seed)*aia_firdem_readnoise(channels(i))
				errs(j,k,i)=sqrt(data(j,k,i)*G + aia_firdem_readnoise(channels(i))^2)
			endfor
		endfor
	endfor
	print, "Data with noise", transpose(data(0,0,*))
	
	return, data

end

pro aia_firdem_wrapper, channels, exptimes, datain, dem_out, tr_struct = tr_struct, $
		niter = niter, its=its, chi2s=chi2s, logthi=logthi, logtlo=logtlo, errors=errors, $
		chi2_reg_probs=chi2_reg_probs, chi2_target=chi2_target, telapsed = telapsed, $
		minconds = minconds

	; Get the temperature response if it's not defined by the calling procedure:
	if(keyword_set(tr_struct) eq 0) then begin
		tr_struct= aia_get_response(/temp,/dn,/evenorm)
	endif

	nchan = n_elements(channels)
	
	data0 = datain > 0

	; Estimate errors in each pixel:	
	if(n_elements(errors) eq 0) then errors = aia_firdem_estimate_error(channels, data0)
		
	if(n_elements(chi2_target) eq 0) then chi2_target = chisqr_cvf(0.05,nchan)
	
	if(n_elements(logtlo) eq 0) then logtlo = 5.5
	if(n_elements(logthi) eq 0) then logthi = 7.5
	
	nb2=round(20*(logthi-logtlo))
				
	if(n_elements(chi2_reg_probs) eq 0) then begin
		chi2_reg_probs = [0.9,0.5,0.01]
		minconds = [0.01,0.05,0.1]
	endif
			
	firdem, data0, errors, exptimes, channels, tr_struct.all, tr_struct.logte, $
		tr_struct.channels, dem_out, its=its, chi2s=chi2s, niter_max=niter, logtmax=logthi, logtmin=logtlo, chi2_target = chi2_target, chi2_reg_probs=chi2_reg_probs, a_struc=a_struc, nb2=nb2, nbad_iter_max=250, telapsed = telapsed, minconds=minconds
		
end

