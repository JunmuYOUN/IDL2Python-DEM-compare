;+
; Name: eis_firdem_wrapper
;
; This is a code for fast computation of differential emission measures from
; EIS data. It operates on a set of EIS images, one for each spectral line
; (e.g., from eis_auto_fit). Requires a file called 'eis_tresp.sav' which
; defines the emissivities for each spectral line. This file is created by
; repeated calls (one for each spectral line) to 'add_eis_line_to_dem.pro'.
;
; Inputs:
;
;	channels: 	The names of the spectral lines corresponding to the set of images 
;				in datain. Must match the channel names found in 
;				eis_tresp.sav.
;
;	datain:		An 3-index [nx,ny,nchannels] array containing a set of EIS 
;				line intensities. The DEM will be computed for each pixel in these 
;				images. Units should be erg/cm^2/s/sr.
;
;	errors:		[nx,ny,nchannels] estimates of the errors corresponding
;				to each element of datain.
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
;
; Optional Arguments:
;
;	tr_struct:	Temperature response structure, as computed by eis_firdem_add_line_to_dem.pro.
;				If not supplied, will attempt to load a temperature response structure
;				from 'eis_tresp_d90.sav' (assumes a density of 10^9.0 g/cm^3). This file
;				can be created by the script 'eis_firdem_tresp_setup_script.pro'.
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
;	minconds:	Minimum condition numbers for SVD used in iterative removel of negative EM. 
;				There should be one for each element of chi2_reg_probs, set below. Generally
;				speaking, a stronger regularization works better accompanied by a more stringent
;				minimum condition. Default: [0.01,0.1,0.2] for the default chi_2_reg_probs of
;				[0.8,0.6,0.01]. Note that if chi2_reg_probs is set but minconds is unassigned,
;				it will instead be set to a constant 0.02 for all input regularization 
;				strengths. If chi2_reg_probs is not set, minconds will be set to the default
;				(i.e., chi2_reg_probs must by assigned for this argument to have any effect).
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
;	--Joseph Plowman (plowman@physics.montana.edu) 09/17/12;
;-

; Reads files listing the EIS effective area and returns the effective area for a requested
; wavelength:
function eis_firdem_effarea, wavelength, filea=filea, fileb=fileb

	if(n_elements(filea) eq 0) then filea = '/ssw/hinode/eis/response/EIS_EffArea_A.004'
	if(n_elements(fileb) eq 0) then fileb = '/ssw/hinode/eis/response/EIS_EffArea_B.005'
	
	areaa=rd_tfile(filea,/nocom,/auto,/conv)
	areab=rd_tfile(fileb,/nocom,/auto,/conv)
	area = [[areab],[areaa]]

	return,interpol(area(1,*),area(0,*),wavelength)
	
end

pro eis_firdem_erg_to_dnphot_factors,waves_in,dnconv_out,phconv_out

	common eis_conversion_factors_fastdem,dnconv,phconv,dnconv2,phconv2
	
	wavmin = 100
	wavmax = 300
	nw = 500
	wavesarr=wavmin+(wavmax-wavmin)*findgen(nw)/(nw-1)
	
	if(n_elements(dnconv) eq 0) then begin
		dnconv = eis_ergs_to_dn(wavesarr,1.0,phot=phconv)
		dnconv2 = spl_init(wavesarr,dnconv)
		phconv2 = spl_init(wavesarr,phconv)
	endif
	
	dnconv_out = spl_interp(wavesarr,dnconv,dnconv2,waves_in)
	phconv_out = spl_interp(wavesarr,phconv,phconv2,waves_in)
	
end

; Computes the conversion from erg/cm^2/s/sr to photons observed by EIS at a given wavelength
; (in angstroms). Not used directly by the EIS DEM inversion, but is instead intended for use 
; by programs calling the DEM inversion to allow estimation of errors from the intensities. 
function eis_firdem_int_to_photons, int, wavelength, exptimes, pixscal, noarea=noarea

	nw = n_elements(wavelength)
	nd = n_elements(size(int,/dimensions))
		
	area = eis_firdem_effarea(wavelength)
	if(n_elements(noarea) ne 0) then area(*)=1.0

	pho = int
	eis_firdem_erg_to_dnphot_factors,wavelength,dnconv,phconv
	dnconv *= pixscal(0)*pixscal(1)
	phconv *= pixscal(0)*pixscal(1)
;	print,wavelength,area,pixscal(0),pixscal(1),photons_to_ergs
	if(nd eq 1) then begin
		if(nw eq 1) then begin
			pho = int*phconv(0)
		endif else begin
			for i=0,nw-1 do pho(i)=int(i)*phconv(i)*exptimes(i)
		endelse
	endif
	if(nd eq 2) then begin
		if(nw ne 1) then begin
			for i=0,nw-1 do pho(*,i) = int(*,i)*phconv(i)*exptimes(i)
		endif else begin
			pho = int*phconv(0)
		endelse			
	endif
	if(nd eq 3) then begin
		for i=0,nw-1 do begin
			pho(*,*,i)=int(*,*,i)*phconv(i)*exptimes(i)
		endfor
	endif
	return,pho
	
end

function eis_firdem_readnoise, wavelength, linewidth=linewidth, pixscal=pixscal, rdn_dn=rdn_dn

	if(n_elements(linewidth) eq 0) then linewidth = 10 ; Default line width in pixels
	
	if(n_elements(rdn_dn) eq 0) then rdn_dn = 2.4 ; Default size of read noise in DN
	
	if(n_elements(pixscal) eq 0) then pixscal=[1,1] ; Default pixel binning
	
	eis_firdem_erg_to_dnphot_factors,wavelength,dnconv,phconv
	
	return, sqrt(linewidth*pixscal(0)*pixscal(1))*rdn_dn/dnconv
	
end
	

; Reverse of eis_firdem_int_to_photons:
function eis_firdem_photons_to_int, pho, wavelength, exptimes, pixscal, noarea=noarea

	nw = n_elements(wavelength)
	
	nd = n_elements(size(pho,/dimensions))
	
	area = eis_firdem_effarea(wavelength)
	if(n_elements(noarea) ne 0) then area(*)=1.0

	int=pho
	eis_firdem_erg_to_dnphot_factors,wavelength,dnconv,phconv
	dnconv *= pixscal(0)*pixscal(1)
	phconv *= pixscal(0)*pixscal(1)
	if(nd eq 1) then begin
		if(nw eq 1) then begin
			int = pho/phconv
		endif else begin
			for i=0,nw-1 do int(i)=pho(i)/phconv(i)/exptimes(i)
		endelse
	endif
	if(nd eq 2) then begin
		if(nw ne 1) then begin
			for i=0,nw-1 do begin
				int(*,i)=pho(*,i)/phconv(i)/exptimes(i)
			endfor
		endif else begin
			int = pho/phconv(0)
		endelse			
	endif
	if(nd eq 3) then begin
		for i=0,nw-1 do begin
			int(*,*,i)=pho(*,*,i)/phconv(i)/exptimes(i)
		endfor
	endif
	
	return,int

end

;+
; Name: eis_firdem_estimate_error
;
; Estimates the EIS error in a data array. data and error units are assumed to be erg/cm^2/s/sr.
;
; Inputs:
;	
;	waves:		The wavelengths of each line contained in the data.
;	data:		The data array (dimensions [nx,ny,n_elements(waves)]) whose error is to be
;				estimated. Units should be erg/cm^2/s/sr.
;	exptimes:	The exposure times for each line in the data array
;
; Optional Inputs:
;
;	pix_scales: 2-element vector giving the x and y scales of the pixels. Default: [2.0,1.0]
;	other_noise:	An estimate of sources of error other than the read and shot noise,
;					Such as that due to the line fitting and uncertainty in emission models.
;					Default: 0.1 (Note that this is a fractional error in the intensities)
;
; Optional Outputs:
;
;	photons:		Estimated photon counts.
;	errors_shot:	Estimated shot noise.
;	errors_nonshot:	Estimated errors from sources other than shot noise.
;
;	--Joseph Plowman (plowman@physics.montana.edu) 09/17/12
;-
function eis_firdem_estimate_error, waves, data, exptimes, pix_scales=pix_scales, $
		other_noise=other_noise, photons=photons, errors_shot=errors_shot, $
		errors_nonshot=errors_nonshot

	if(n_elements(pix_scales) eq 0) then pix_scales = [2.0,1.0]
	if(n_elements(other_noise) eq 0) then other_noise = 0.1
	nchan = n_elements(waves)

	; Convert to photons to estimate shot noise:
	photons = eis_firdem_int_to_photons(data,waves,exptimes,pix_scales) > 0
	
	; Convert shot noise back from photons:
	errors_shot = eis_firdem_photons_to_int(sqrt(photons),waves,exptimes,pix_scales)
	
	; Compute read and other noise:
	rdn = eis_firdem_readnoise(waves,pixscal=pix_scales)/exptimes
	errors_nonshot = other_noise*data
	for i=0,nchan-1 do begin
		errors_nonshot(*,*,i) = sqrt((errors_nonshot(*,*,i))^2+(rdn(i))^2)
	endfor

	; Combine errors in quadrature:
	return, sqrt((errors_shot)^2+(errors_nonshot)^2)	
	
end

;+
; Name: eis_firdem_synthesize_data
;
; Returns an array of synthetic EIS data from an input DEM, with randomly generated read and
; shot noise. Array will have dimensions [nx,ny,n_elements(channels)]
;
; Inputs:
;	
;	logt:		The log_10(T) axis for the DEM.
;	dem_in:		The input DEM, assumed to have units 1/cm^5/log_10(T).
;	channels:	The line names to use. Line names should be in the same format as in
;				tr_struct.eis_windows, and each line must be included in tr_struct.
;	exptimes:	Array of exposure times for each channel.
;	tr_struct:	Temperature response structure, as computed by add_eis_line_to_dem.pro.
;
; Optional Inputs:
;	
;	nx:			x dimension of pixel array to generate. Default: 100
;	ny:			y dimension of pixel array to generate. Default: 100
;	pix_scales: 2-element vector giving the x and y scales of the pixels. Default: [2.0,1.0]
;	other_noise:	An estimate of sources of error other than the read and shot noise,
;					Such as that due to the line fitting and uncertainty in emission models.
;					Default: 0.1 (Note that this is a fractional error in the intensities)
;
; Optional Outputs:
;	
;	errs:		An array containing estimated errors for each element of the output data
;				array.
;	data0:		Same as the output data array, but without errors.
;
;	--Joseph Plowman (plowman@physics.montana.edu) 09/20/12
;-
function eis_firdem_synthesize_data, logt, dem_in, channels, exptimes, tr_struct=tr_struct, $
		nx=nx, ny=ny, errs=errs, data0=data0, pix_scales=pix_scales, other_noise=other_noise

	; Get the temperature response if it's not defined by the calling procedure:
	if(n_elements(tr_struct) eq 0) then begin
		trfile = 'eis_tresp_d90.sav'
		print,"Loading temperature response from ", trfile
		restore,trfile
	endif

	channels_all = tr_struct.eis_windows
	nchan_all = n_elements(channels_all)
	nchan = n_elements(channels)
	tr_logt=tr_struct.logte
	tr0 = tr_struct.all
	waves = dblarr(nchan)
	nt = n_elements(logt)
	tr = dblarr(nt,nchan)
	data0 = dblarr(nx,ny,nchan)

	; Compute noise-free intensities from the DEM:
	for i=0,nchan_all-1 do begin
		for j=0,nchan-1 do begin
			if channels_all(i) eq channels(j) then begin
				tr02 = spl_init(tr_logt,tr0(*,i),/double)
				tr(*,j) = spl_interp(tr_logt,tr0(*,i),tr02,logt)
				tr = tr > 0
				integrand = dem_in*tr(*,j)
				data0(*,*,j) += int_tabulated(logt,integrand)
				waves(j) = tr_struct.line_centers(i)
			endif
		endfor
	endfor

;	if(n_elements(pix_scales) eq 0) then pix_scales = [2.0,1.0]
;	if(n_elements(other_noise) eq 0) then other_noise = 0.01
;	photons0 = int_to_eis_photons(exptimes*data0,waves,pix_scales)
;	photons = photons0
;	errs = eis_photons_to_int(sqrt(photons0),waves,pix_scales)/exptimes
;	rdn = get_eis_read_noise(waves,pixscal=pix_scales)/exptimes
;	errors_nonshot = other_noise*data0
;	for i=0,nchan-1 do begin
;		errors_nonshot(*,*,i) = sqrt((errors_nonshot(*,*,i))^2+(rdn(i))^2)
;	endfor
;	errs = sqrt((errs)^2+(errors_nonshot)^2)

	; Compute errors and photon counts:
	errs = eis_firdem_estimate_error(waves, data0, exptimes, pix_scales=pix_scales, $
			other_noise=other_noise, photons=photons, errors_shot=errors_shot, $
			errors_nonshot=errors_nonshot)
	
	; Generate Poisson distributed photon counts:
	for j=0,nx-1 do begin
		for k=0,ny-1 do begin
			for i=0,nchan-1 do begin
				if(photons(j,k,i) gt 0) then begin
					photons(j,k,i) = randomu(seed,poisson=photons(j,k,i))
				endif
			endfor
		endfor
	endfor
	
	data = eis_firdem_photons_to_int(photons, waves, exptimes, pix_scales) + $
			errors_nonshot*randomn(seed,nx,ny,nchan)
	
	; Estimate the errors using the noisy intensities (since in reality we'll only
	; know the noisy ones):
	errs = eis_firdem_estimate_error(waves, data, exptimes, pix_scales=pix_scales, $
			other_noise=other_noise, photons=photons, errors_shot=errors_shot, $
			errors_nonshot=errors_nonshot)
			
	return, data
	
end


pro eis_firdem_wrapper, channels, datain, errors, dem_out, tr_struct = tr_struct, $
		niter = niter, its=its, chi2s=chi2s, minconds=minconds, exptimes=exptimes, $
		logthi=logthi, logtlo=logtlo, chi2_reg_probs=chi2_reg_probs, $
		chi2_target=chi2_target, nb2=nb2, telapsed=telapsed

	; Get the temperature response if it's not defined by the calling procedure:
	if(n_elements(tr_struct) eq 0) then begin
		trfile = 'eis_tresp_d90.sav'
		print,"Loading temperature response from ", trfile
		restore,trfile
	endif

	nchan = n_elements(channels)
	
	data0 = datain > 0
	
	if(n_elements(nb2) eq 0) then nb2 = 80

	exptimes = 1.0+dblarr(nchan) ; Exposure times have been divided out with these units...
		
	if(n_elements(chi2_target) eq 0) then chi2_target = chisqr_cvf(0.05,nchan)

	if(n_elements(chi2_reg_probs) eq 0) then begin
		chi2_reg_probs = [0.9,0.5,0.01]
		minconds = [0.01,0.05,0.1]
	endif
			
	if(n_elements(nbad_iter_max) eq 0) then nbad_iter_max = 200
		
	firdem, data0, errors, exptimes, channels, tr_struct.all, tr_struct.logte, $
		tr_struct.eis_windows, dem_out, its=its, chi2s=chi2s, niter_max=niter, logtmax=logthi, $
		logtmin=logtlo, chi2_target = chi2_target, chi2_reg_probs=chi2_reg_probs, $
		a_struc=a_struc, nbad_iter_max = nbad_iter_max, nb2 = nb2, minconds=minconds, $
		telapsed=telapsed
		
end


