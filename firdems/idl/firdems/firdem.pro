;+
; Name: firdem
;
; Fast, iterative, regularized routine for computation of differential emission measures (DEMs)
; using density insensitive spectral data from optically thin plasmas.
;
; Inputs:
;
;	datain: Array (dimensions [nx,ny,nchan]) containing intensities in each of the spectral
;			channels/lines from which the DEM will be computed. Note that negative data
;			elements are not zeroed by this routine.
;	errsin: The uncertainties corresponding to datain, with the same array dimensions.
;	exptimes: Exposure times for each channel.
;	data_channeltags: Names of each channel contained in datain.
;	tresps: Array (dimensions [n_tresps,n_temps]) containing a set of temperature response
;			functions which must include all those found in the input data.
;	tresp_logt: Array (dimensions [n_temps]) giving the logt axis for the tresps array.
;	tresp_channeltags: Array (dimensions [n_tresps]) labeling each channel contained
;						in tresps. Labels must be a superset of data_channeltags.
;
; Outputs:
;
;	dem_out: Structure containing the output DEMs for each pixel in datain. Has the following
;			 elements:
;					t(nt):	Array of temperatures at which the basis elements 
;						have been evaluated.
;					basis(nt,n_bases): The amplitudes of the basis elements, in
;						units of differential emission measure (cm^-5 per unit
;						alog10(T)). The number of basis elements used is set by
;						the optional input 'nb2'.
;					coffs(nx,ny,n_bases): Sets of basis element coefficients, which define the
;						DEM for each pixel in the input data. May have small residual negative
;						elements, which the user may safely zero.
;					t0a(n_bases): Array listing the centers of each basis element
;
; Optional Inputs:
;	
;	niter_max: The maximum number of iterations to attempt when iterating away negative
;				emission. Default (set in the subroutine 'firdem_iterate'): 2000
;
;	logtmin: Minimum of temperature range for DEM and response fuctions. Default (set in
;			firdem_firstpass): 5.5
;
;	logtmax: Maximum of temperature range for DEM and response fuctions. Default (set in
;			firdem_firstpass): 7.5
;
;	nb2: Number of higher resolution basis elements. Default: 47
;	
;	nt: Number of temperature grid points. Default: (nb2+1)*10+1
;
;	mincond1: Minimum condition number for SVD used in regularized first-pass DEM inversion.
;           Limited by machine precision. Default (set in firdem_firstpass): 1.0e-12
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
;	nbad_iter_max: Stop iterating (at this chi2_reg_probs level) if this many iteration steps
;			have been taken without improving chi squared. Default (set in firdem_iterate): 200.
;
;	gaussbasis, squarebasis:	If set, use Gaussian or box functions, respectively, for the
;								higher-resolution basis instead of the default triangle
;								functions. Note: squarebasis supersedes gaussbasis
;
;
; Optional Outputs:
;
;	its(nx,ny): Total number of iterations taken for each pixel;
;
;	chi2s(nx,ny): Reduced chi squared residuals for each pixel's final DEM.
;
;	a_struc:	Contains information regarding the mapping from DEMs (expressed in terms of
;				coefficients of the basis functions) to intensities, and vice versa. See
;				comments for firdem_firstpass.
;
;	telapsed:	The total elapsed time for computation of the DEMs.
;
; Joseph Plowman	(plowman@physics.montana.edu) 09-27-12
;
;-

; This function returns a triangle with base b and height height centered at x0, evaluated
; at the points contained in x:
function firdem_triangle,x,x0,b,h

	sep = abs(x-x0)

	y = h*(0.5*b-sep)/(0.5*b)

	negs = where(y lt 0, count)
	if(count gt 0) then y(negs) = 0
	
	return,y
	
end

; This procedure computes a basis of nb triangle functions on the grid defined by logt.
; Results are returned in basis(nt,nb), and the optional argument t0a(nb) contains the
; locations of the peaks.
pro firdem_triangle_basis,nb,logt,basis,t0a=t0a

	logtlo = min(logt)
	logthi = max(logt)
	logdt = (logthi-logtlo)/(nb+1)
	
	nt = n_elements(logt)
	t0a = logdt+logtlo+logdt*findgen(nb)
	
	logtm = logt-logdt
	logtp = logt+logdt
	basis = dblarr(nt,nb)
	
	for j=0,nb-1 do begin
		basis(*,j) = firdem_triangle(logt,t0a(j),2.0*logdt,1.0)
	endfor
		
end
		
; This procedure computes a basis of nb gaussians functions on the grid defined by logt.
; Results are returned in basis(nt,nb), and the optional argument t0a(nb) contains the
; locations of the peaks. Standard deviations are equal to half the grid spacing.
pro firdem_gauss_basis,nb,logt,basis,t0a=t0a

	logtlo = min(logt)
	logthi = max(logt)
	logdt = (logthi-logtlo)/(nb+1)
	sig = 0.5*logdt
	pi = 2*acos(0)	
	nt = n_elements(logt)
	t0a = logdt+logtlo+logdt*findgen(nb)
	
	logtm = logt-logdt
	logtp = logt+logdt
	basis = dblarr(nt,nb)
	
	for j=0,nb-1 do begin
		basis(*,j)=logdt*exp(-0.5*((logt-t0a(j))/(sig))^2)/sqrt(2*pi*sig^2)
	endfor

end

; This procedure computes a basis of nb box functions on the grid defined by logt.
; Results are returned in basis(nt,nb), and the optional argument t0a(nb) contains the
; locations of the box centers.
pro firdem_square_basis,nb,logt,basis,t0a=t0a

	logtlo = min(logt)
	logthi = max(logt)
	logdt = (logthi-logtlo)/(nb+1)
	
	nt = n_elements(logt)
	t0a = logdt+logtlo+logdt*findgen(nb)
	
	logtm = logt-logdt
	logtp = logt+logdt
	basis = dblarr(nt,nb)
	
	for j=0,nb-1 do begin
		tin = where((logt ge t0a(j)-0.5*logdt) and (logt lt t0a(j)+0.5*logdt))
		basis(tin,j)=1.0
	endfor

end

; Name: firdem_regularize_data:
;
; Performs the regularization of an input data vector.
;
; Inputs:
;	
;	datavec0:	The input data vector.
;	sigmas: 	The estimated standard deviations in datavec0.
;	a_struc: 	Structure containing information pertaining to the response functions.
;				inversion matrices, response function norms, etc.
;
; Optional inputs:
;
;	chi2_end:	Desired final chi squared relative to the inital input data. Should be
;				consistent with typical random variation in the data, e.g.,
;				chisqr_cvf(0.5,nchan).
;	chi2_tol:	Regularization is considered complete once the fractional difference between 
;				the current chi squared and chi2_end is less than chi2_tol. Typically 0.1 or
;				less - smaller values take longer, but result in a more consistent
;				regularization.
;	niter_max:	Maximum number of iterations to perform. More iterations take longer, but may
;				be required to reach a given chi2_end if chi20 or maxdelta are small. 50 is a
;				typical value.
;
; Joseph Plowman (plowman@physics.montana.edu) 07-12-11 
;
;+
function firdem_regularize_data, datavec0, sigmas, a_struc, chi2_end=chi2_end, $
		chi2_tol=chi2_tol, niter_max=niter_max

	nchan = n_elements(datavec0)
	; Check to see if the data values are so small that zero is an acceptable solution:
	if(total((datavec0/sigmas)^2) lt chi2_end) then return,dblarr(nchan)
	if(chi2_end/nchan lt 1.0e-4) then return,datavec0 ; Essentially no regularization required.
	data2vec = dblarr(nchan)
	deltadata = dblarr(nchan)
	
	; Defaults for optional inputs:
	if(n_elements(chi2_end) eq 0) then chi2_end = chisqr_cvf(0.5,nchan)
	if(n_elements(chi2_tol) eq 0) then chi2_tol = 0.05
	if(n_elements(niter_max) eq 0) then niter_max = 50

	tr_norm = a_struc.tr_norm*a_struc.exptimes
	a_inv = a_struc.a_inv	
	
	; a_array is computed with normalized temperature response functions so that it has unit
	; diagonals and off-diagonal elements less than one. The data vector must be divided by
	; the same normalization factor to match.
	sigsmallest = min(sigmas)
	datavec = reform((datavec0)/tr_norm)
	sigsmall = 0.1*sigmas
	a_inv_scaled = a_inv/(tr_norm#tr_norm)
	sigs2_diag = diag_matrix(reform(1.0/sigmas))^2
	lastalpha = 0.0
	lastchi2 = 0.0
	
	chi2_low = 0.0
	alpha_low = 0.0
	bisect_start = 0
	
	for k=0,niter_max-1 do begin
		datavec(*) = datavec0/tr_norm
		; Intermediate vector appearing multiple times in the calculation:
		data2vec(*) = (a_inv#(datavec))/tr_norm 
		; Multiplicative factor ensuring chi squared between the regularized data vectors
		; for this step vs. the previous step is chi20:
		if(k eq 0) then alpha = sqrt(9.0*chi2_end/(total(sigmas^2*data2vec^2)))
		; The difference between the regularized data vector for this step and the previous
		; step:
		b_inv = invert(alpha*a_inv_scaled+sigs2_diag,/double)
		datavec(*) = datavec-(alpha*b_inv#data2vec)/tr_norm
		; Finish if chi squared between datavec and datavec0 reaches chi2_end:
		chi2=total(((datavec0-datavec*tr_norm)/sigmas)^2)
		if(k eq 0) then chi2_high = chi2

		if(bisect_start eq 1) then begin
			if(chi2 lt chi2_end) then alpha_low = alpha
			if(chi2 gt chi2_end) then alpha_high=alpha
			alpha = alpha_low+0.5*(alpha_high-alpha_low)
		endif

		if(chi2 lt chi2_end and bisect_start eq 0) then alpha = 5*alpha

		if(chi2 gt chi2_end and bisect_start eq 0) then begin
			bisect_start = 1
			chi2_high = chi2
			alpha_high = alpha
			alpha = alpha_low+0.5*(alpha_high-alpha_low)
		endif


		if(abs(chi2-chi2_end)/chi2_end lt chi2_tol) then break
	endfor

	return,datavec*tr_norm
	
end

;+
; Name: firdem_firstpass
;
; This code computes the `first pass' regularized DEMs for firdem,
; along with matrices used for the applying and inverting the convolution that maps
; the DEM (with instrument response as a basis) to channel/spectral line intensities. It
; operates on a set of images, one for each channel.
;
; Inputs:
;
;	datain:	An 3-index [nx,ny,nchannels] array containing a set of AIA 
;			images. The DEM will be computed for each pixel in these 
;			images. Units should be DN.
;
;	errsin: Array with same dimensions as datain which specifies the errors in the data.
;
;	exptimes:	The exposure times for each image.
;
;	data_channeltags: Names of each channel contained in datain.
;
;	tresps: Array (dimensions [n_tresps,n_temps]) containing a set of temperature response
;			functions. Must include all those found in the input data.
;
;	tresp_logt: Array (dimensions [n_temps]) giving the logt axis for the tresps array.
;
;	tresp_channeltags: Array (dimensions [n_tresps]) labeling each channel contained
;						in tresps. Labels must be a superset of data_channeltags.
;
; Outputs:
;
;	dem_out:	Structure containing the results of the DEM computation. Has
;				the following elements:
;					t(nt):	Array of temperatures at which the basis elements 
;						have been evaluated.
;					basis(nt,n_bases): The amplitudes of the basis elements, in
;						units of differential emission measure (cm^-5 per unit
;						alog10(T)). There are as many basis elements as channels.
;					coffs(nx,ny,n_bases): The coefficients of each basis element,
;						which define the DEM. May have small residual negative
;						elements, which the user may safely zero.
;
; Optional Inputs:
;
;	logtmin: Minimum log temperature (default is 5.5).
;
;	logtmax: Maximum log temperature (default is 7.5).
;
;	nt: Number of temperature points in interpolation grid for response functions 
;	(default is 100).
;
;	a_struc:	Contains information regarding the mapping from DEMs (expressed in terms of
;				coefficients of the basis functions) to intensities, and vice versa. Has the
;				following elements:
;					t: 			Array of temperatures; equivalent to dem_out.t.
;					tr:			An array containing the temperature response functions,
;								interpolated at the points specified by t. Essentially
;								the input tresps array interpolated to a_struc.t.
;					tr_norm:	Normalization of response functions, chosen to be
;								tr_norm(i)=sqrt(integral(R_i^2 dlog_10(T))).
;					basis:		Array listing the basis functions. Chosen to be the 
;								response functions divided by the corresponding tr_norm.
;					a_array:	matrix which maps from DEMs coefficients to intensities. It is
;								formed by convolving the basis functions with the response
;								functions, i.e., a_ij = integral(R_i*B_j/tr_norm(i) dlog_10(T)).
;								By construction, its diagonal elements are unity, and off
;								diagonals are positive and less than or equal to one.
;					v,u,wp:		matrices resulting from the SVD of a_array. Defined so that the 
;								(pseudo)inverse of a_array is equal to v#wp#u. wp reflects the
;								minimum condition number given by mincond1
;					wvec:		The original vector of singular values returned by the SVD.
;								Allows the pseudoinverse to be recomputed with a different
;								choice of minimum condition number (i.e., a_inv = v#wp#u where
;								wp is formed by taking diag_matrix(wvec) and zeroing entries
;								whose corresponding to elements of wvec less than									;								mincond*max(wvec)).
;					a_inv: 		Inverse (or pseudo-inverse, depending on condition of a_array)
;								matrix of a_array corresponding to wp. Used for the regularized
;								first pass inversion.
;
;	mincond1: Minimum condition number for SVD used in regularized first-pass DEM inversion.
;           Limited by machine precision. Default: 1.0e-12
;
;	noreg_flag: Don't regularized data if this keyword is set.
;
;	chi2_reg_end: Regularize the data so that it has (approximately) this chi squared.
;
;-
pro firdem_firstpass, datain, errsin, exptimes, data_channeltags, tresps, tresp_logt, $
		tresp_channeltags, dem_out, logtmin=logtmin, logtmax=logtmax, nt=nt, $
		a_struc=a_struc, mincond1=mincond1, noreg_flag=noreg_flag, $
		chi2_reg_end=chi2_reg_end

	; Set up the a structure described above:
	if(n_elements(a_struc) eq 0) then begin
		
		; Default minimum and maximum temperatures for the DEM inversion.
		if(n_elements(logtmin) eq 0) then logtmin=5.5
		if(n_elements(logtmax) eq 0) then logtmax=7.5
		
		tmin = 10^logtmin
		tmax = 10^logtmax
		
		; Default minimum condition number for the regularized first pass inversion. Due
		; to the regularized inversion, mincond1 is limited primarily by the machine precision.
		if(n_elements(mincond1) eq 0) then mincond1 = 1.0e-12
			
		; Determine which channels are in the data and downselect the temperature response
		; functions accordingly:
		nchan = n_elements(data_channeltags)
		logt = tresp_logt
		nt0 = n_elements(logt)

		tr0 = dblarr(nt0,nchan)
		n_tresps = n_elements(tresp_channeltags)
		ichan = 0
		for i=0,nchan-1 do begin
			for j=0,n_tresps-1 do begin
				if(data_channeltags(i) eq tresp_channeltags(j)) then begin
					tr0(*,ichan) = tresps(*,j)
					ichan++
				endif
			endfor
		endfor

		; Interpolate the temperature response functions to nt points uniformly spaced
		; in Log(T):
		if(n_elements(nt) eq 0) then nt=100
		logt2 = alog10(tmin)+(alog10(tmax)-alog10(tmin))*findgen(nt)/(nt-1)
		tr = dblarr(nt,nchan)
		for i=0,nchan-1 do begin
			tr02 = spl_init(logt,tr0(*,i),/double)
			tr(*,i) = spl_interp(logt,tr0(*,i),tr02,logt2)
		endfor
		t=10^logt2
		tout=10^logt2
		
		basis=tr ; Use the temperature response functions as the basis.
		nb = n_elements(basis(0,*))
	
		; Array Setup:
		a_array = dblarr(nchan,nb)
		basis_norm = dblarr(nb)
		tr_norm = dblarr(nb)
		flat_resp_intens = dblarr(nb)
		
		; Normalize the basis elements so that their squared integral over
		; Log_10(T) is one:
		for i=0,nb-1 do begin
			basis_norm(i) = sqrt(int_tabulated(logt2,basis(*,i)*basis(*,i)))
			basis(*,i) /= basis_norm(i)
			flat_resp_intens(i) = int_tabulated(logt2,basis(*,i))
		endfor
		
		; Compute the norm of the temperature response functions (same as for basis):
		for i=0,nchan-1 do tr_norm(i) = sqrt(int_tabulated(logt2,tr(*,i)*tr(*,i)))
		
		; Compute the a matrix, which maps from the DEM (expressed as coefficients of the basis
		; functions) to the intensities:
		for i=0,nb-1 do begin
			for j=0,nchan-1 do begin
				integrand = basis(*,i)*tr(*,j)/tr_norm(j)
				a_array(j,i) = int_tabulated(logt2,integrand)
			endfor
		endfor
		
		; Compute the SVD of the a matrix:
		svdc,a_array,wvec,u,v,/double,/column
		
		; Threshold the diagonal w matrices, returned by the SVD, to enforce the minimum 
		; condition numbers:
		smallw = where(wvec lt mincond1*max(wvec),count)
		wpvec=1.0/wvec
		if(count ne 0) then begin
			print,count," degrees of freedom discarded from inversion due to poor condition"
			wpvec(smallw)=0
		endif
		wp = diag_matrix(wpvec)
		
		flat_coffs=v#wp#transpose(u)#flat_resp_intens
	
		; Compute the inverses of the a matrix corresponding to wp:
		a_inv = v#wp#transpose(u)
		
		; Create the a structure:
		a_struc = {v:v,u:transpose(u), wp:wp, t:tout, basis:basis, tr_norm:tr_norm, a_array:a_array, nchan:nchan, nb:nb, a_inv:a_inv, tr:tr, flat_coffs:flat_coffs, exptimes:exptimes, wvec:wvec}
		
	endif

	; Make sure that the input data has dimensions (nx,ny,nchan):
	nchan = n_elements(data_channeltags)
	if(n_elements(datain) eq nchan) then begin
		data = reform(datain,1,1,nchan)
		errors = reform(errsin,1,1,nchan)
	endif
	if(n_elements(datain(0,*,*)) eq nchan) then begin
		data = reform(datain,n_elements(datain(*,0,0)),1,nchan)
		errors = reform(errsin,n_elements(errsin(*,0,0)),1,nchan)
	endif
	if(n_elements(datain(0,0,*)) eq nchan) then begin
		data = datain
		errors = errsin
	endif
	
	; Zero negative elements (unless no regularization has been requested).
;	if (keyword_set(noreg_flag) eq 0) then data = data > 0
	
	nx = n_elements(data(*,0,0))
	ny = n_elements(data(0,*,0))
	dems = dblarr(nx,ny,nb)
	
	; Loop over each individual pixel (slightly inefficient):
	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			datavec0 = reform(data(i,j,*))
			datavec = datavec0
			sigmas = reform(errors[i,j,*])
			if(keyword_set(noreg_flag) eq 0) then begin
				datavec=firdem_regularize_data(datavec0, sigmas, a_struc, chi2_end=chi2_reg_end)
			endif
			dems(i,j,*) = a_inv#(datavec/tr_norm)
		endfor
	endfor

	; Create the DEM output structure:
	dem_out = {t:a_struc.t, basis:a_struc.basis, coffs:dems, flat_coffs:a_struc.flat_coffs}

end

; +
; Name: firdem_iterate
;
; Iterative removal of negative emission measure for a single pixel.
;
; Inputs:
;
;	dem_in: The initial DEM coefficients in the higher resolution basis.
;	sigmas: The uncertainties in the data corresponding to the DEM.
;	a: Structure returned by firdem_firstpass containing information about the a matrix and 
;		its inverse.
;	basis22: matrix which maps from instrument response basis to the higher resolution basis.
;	a2_array: matrix which maps from the higher resolution basis to the intensities.
;
; Outputs:
;
;	dem_out: High resolution basis coefficients of final DEM. May contain some residual negative
;			 emission measure, which the user may set to zero. If the iteration was successful,
;			 the resulting all-positive DEM will match the data to within reduced chi-squared
;			 of chi2thold.
;
; Optional inputs:
;
;	datain: The data values which the DEM will be iterated towards, removing zeroes in the
;			process. If not supplied, it will be computed from the input DEM, and
;			the iteration will only try to iterate away any negative emission in the input
;			DEM while staying consistent with that data.
;	niter: Absolute maximum number of iterations to take. Default: 2000.
;	chi2thold: Threshold chi squared, relative to input data, for finishing iteration.
;			Default: chisqr_cvf(0.2,nchan)/nchan.
;	nbad_max: Stop if this many iteration steps have been taken without improving chi squared.
;			Default: 200.
;	extrap_fac: Extrapolation factor (fraction of current iteration number). 
;			Default: 0.025.
;	extrap_start: Don't start extrapolating until after this many regular steps. 
;			Default: 1/extrap_fac (no point in extrapolating by less than one step).
;
;	mincond:	The minimum condition number to impose on the inverse matrix used in the
;				iteration. Default: 0.02
;
; Optional outputs:
;
;	datapos: array of intensities corresponding to the positive DEM.
;	dempos_out: High resolution basis coefficients of final DEM, with negatives zeroed
;	chi2_final: Final chi squared of datapos (relative to input data).
;	itcount: number of iterations taken.
;
; Joseph Plowman (plowman@physics.montana.edu) 09-17-12
;
;+
pro firdem_iterate, dem_in, sigmas, a, basis22, a2_array, dem_out, datapos=datapos, dempos_out=dempos_out, chi2_final=chi2_final, itcount=itcount, niter=niter, chi2thold=chi2thold, nbad_max=nbad_max, extrap_start=extrap_start, extrap_fac=extrap_fac, datain=datain, mincond=mincond
	
	nb2 = n_elements(a2_array(0,*))
	nchan = n_elements(a2_array(*,0))

	; Defaults for optional inputs:
	if(n_elements(niter) eq 0) then niter=1000
	if(n_elements(chi2thold) eq 0) then chi2thold = chisqr_cvf(0.2,nchan)/nchan
	if(n_elements(nbad_max) eq 0) then nbad_max = 200
	if(n_elements(extrap_fac) eq 0) then extrap_fac = 0.02
	if(n_elements(extrap_start) eq 0) then extrap_start = ceil(1.0/extrap_fac)
	
	; Set up for iteration: 
	dem = dem_in ; Compute basis2 coefficients
	dem_out = dem
	dempos = dem > 0
	dempos_out = dempos
	if(n_elements(datain) eq 0) then datain = reform(a2_array#dem) ; Data values for DEM
	data=datain
	datapos = reform(a2_array#dempos) ; Data values for DEM with negatives zeroed
	deltadata = datapos-data

	if(n_elements(mincond) eq 0) then mincond = 0.02
	
	wpvec = 1.0/a.wvec
	wsmall = where(a.wvec lt mincond*max(a.wvec), count)
	if(count gt 0) then wpvec(wsmall) = 0.0
	wp = diag_matrix(wpvec)
	ainv_arr = a.v#wp#a.u

	normfac = a.tr_norm*a.exptimes
	coffs = reform(ainv_arr#(data/normfac))	
	deltadem = dblarr(nb2)
	
	nbad_recent=0 ; Counter for steps since last chi squared improvement.

	; When extrapolating, DEMs are initially assigned to these. If they improve chi squared,
	; they are copied into the main DEM and data arrays: 
	dem_test = dblarr(nb2)
	dempos_test = dblarr(nb2)
	datapos_extrap = dblarr(nchan)

	chi2totarr = 1.0/nchan+dblarr(nchan)
	deltatot = dblarr(1,nchan)
	
	; Compute inital chi squared:
	deltatot(0,*) = ((datapos)/sigmas)^2
	chi2min = deltatot^2#chi2totarr
	chi2last=chi2min

	; Begin iteration:
	for itcount=0,niter-1 do begin
		; Stop if we've gone over nbad_recent points since last chi squared improvement:
		if(nbad_recent gt nbad_max) then break

		; Compute DEM correction coefficients in the instrument response basis.
		; These will restore the DEM to agreement with the data, but generally 
		; introduce some negative emission.
		coffs(*) = ainv_arr#(deltadata/normfac)
		
		; Convert these coefficients to basis2:
		deltadem(*) = (basis22#coffs)
		dem(*) = dempos - deltadem; + 0.01*dem_sigmas*randomn(seed,nb2)
		dempos(*) = dem > 0 ; Zero negative basis2 coefficients.
		
		; Compute difference between data for positive DEM and initial data:
		datapos(*) = a2_array#dempos
		deltadata(*) = datapos-data

		; Compute chi squared for this step:
		deltatot(0,*) = ((deltadata)/sigmas)^2
		chi2 = deltatot#chi2totarr
		
		dem_noextrap=dem
		chi2_noextrap=chi2
		
		if(itcount ge extrap_start+1) then begin
			
			; If current chi squared is better than the previous one, extrapolate. Otherwise,
			; do nothing:
			if(chi2 lt lastchi2) then begin
				; Simple linear extrapolation scheme: Change is computed between
				; the most recent DEM and the previous DEM (before extrapolation).
				; This is multiplied by extrap_fac*itcount and added to the most
				; recent DEM:
				deltadem(*) = (dem - lastdem)*extrap_fac*itcount
				dem_test(*)=dem+deltadem
				dempos_test(*) = dem_test > 0
				
				; Compute chi squared resulting from the extrapolation:
				datapos_extrap(*) = a2_array#dempos_test
				deltatot(0,*) = ((datapos_extrap-data)/sigmas)^2
				chi2_test = deltatot#chi2totarr
								
				; If the resulting chi squared is better than the previous one,
				; update the DEM with the extrapolation:
				if(chi2_test lt chi2) then begin
					chi2=chi2_test
					datapos(*) = datapos_extrap
					deltadata(*) = datapos-data
					dem(*)=dem_test
					dempos(*)=dempos_test
				endif
			endif

		endif 
		; If chi squared got better, update output chi squared and dems:
		if(chi2 lt chi2min) then begin
			nbad_recent=0 
			chi2min=chi2
			dem_out=dem 
			dempos_out=dempos
			if(chi2 lt chi2thold) then break ; Finish if chi squared is good enough.
		endif
		
		lastdem=dem_noextrap
		lastchi2=chi2_noextrap
		nbad_recent++
	endfor

	; Finished iterating. dem contains the DEM corresponding to the best chi squared
	; (with negatives zeroed) found. Although chi squared is calculated with negatives
	; zeroed, the returned DEM retains them for good measure (the user will generally
	; want to zero the negatives before using the DEMs; a trivial exercise):
	datapos = a2_array#dempos_out
	deltatot(0,*) = ((datapos-data)/sigmas)^2
	chi2_final = deltatot#chi2totarr
	
end


pro firdem, datain, errsin, exptimes, data_channeltags, tresps, tresp_logt, $
		tresp_channeltags, dem_out, niter_max = niter_max, its=its, chi2s=chi2s, $
		logtmax=logtmax, logtmin=logtmin, nb2=nb2, nt=nt, mincond1=mincond1, $
		minconds=minconds, chi2_reg_probs=chi2_reg_probs, chi2_target=chi2_target, $
		nbad_iter_max=nbad_iter_max, a_struc=a_struc, telapsed=telapsed, $
		gaussbasis=gaussbasis, squarebasis=squarebasis

	tstart=systime(1) 	; Record the start time for benchmarking purposes...

	nchan = n_elements(data_channeltags)
	nx = n_elements(datain(*,0,0))
	ny = n_elements(datain(0,*,0))

	if(n_elements(nb2) eq 0) then nb2=47 ; Number of basis functions.
	if(n_elements(nt) eq 0) then nt=(nb2+1)*10+1; Number of temperature grid points...
	
	; Setup for computing 'first pass' DEMs. Here, firdem_firstpass is only used
	; to create the 'a_struc' structure:
	firdem_firstpass, datain, errsin, exptimes, data_channeltags, tresps, tresp_logt, $
			tresp_channeltags, dems, logtmin=logtmin, logtmax=logtmax, nt=nt, $
			a_struc=a_struc, mincond1=mincond1, /noreg_flag

	normfac = a_struc.tr_norm*a_struc.exptimes
	logt = alog10(a_struc.t)
	tr = a_struc.tr

	its = dblarr(nx,ny) ; Number of iterations for each pixel;
	chi2s = dblarr(nx,ny) ; chi squared for each pixel.
	
	nb = n_elements(dems.basis(0,*)) ; Number of instrument response basis elements.

	
	; Higher-resolution basis of gaussian functions:
	if(n_elements(gaussbasis) ne 0) then firdem_gauss_basis,nb2,logt,basis2,t0a=t0a
	; Higher-resolution basis of square functions:
	if(n_elements(squarebasis) ne 0) then firdem_square_basis,nb2,logt,basis2,t0a=t0a
	; Higher-resolution basis of triangle functions:
	if(n_elements(basis2) eq 0) then firdem_triangle_basis,nb2,logt,basis2,t0a=t0a

	dem_out = {t:10^logt, basis:basis2, coffs:fltarr(nx,ny,nb2), t0a:t0a}

	a2_array = dblarr(nchan,nb2) ; This array maps basis2 coefficients to channel intensities.
	for i=0,nchan-1 do begin
		for j=0,nb2-1 do a2_array(i,j) = int_tabulated(logt,tr(*,i)*basis2(*,j)*exptimes(i))
	endfor
	
	; This array maps DEMs in the instrument response basis ('dems.basis') to the high
	; resolution basis ('basis2'):
	basis22 = dblarr(nb2,nb)
	for i=0,nb-1 do begin
		basis0_derivs = spl_init(logt,dems.basis(*,i),/double)
		basis22(*,i) = spl_interp(logt,dems.basis(*,i),basis0_derivs,t0a)
	endfor
	
	; The regularization chi squareds which will be tried, in order, stopping once acceptable
	; chi squared has been reached:
	if(n_elements(chi2_reg_probs) eq 0) then begin
		chi2_reg_probs = [0.9,0.5,0.01]
		minconds = [0.01,0.05,0.1]
	endif
	nconfs = n_elements(chi2_reg_probs)
	if(n_elements(minconds) eq 0) then minconds = 0.02+dblarr(nconfs)
	chi2_reg_ends = dblarr(nconfs)
	for i=0,nconfs-1 do chi2_reg_ends(i) = chisqr_cvf(chi2_reg_probs(i),nchan)
	
	if(n_elements(chi2_target) eq 0) then chi2_target = chisqr_cvf(0.05,nchan)
	
	; If the initial regularization chi squared is greater than the target chi squared
	; (rendering it unachievable), decrease it so that it is less than the target chi squared:
	if(chi2_reg_ends(0) gt 0.7*chi2_target) then chi2_reg_ends(0) = 0.7*chi2_target
	; Stopping chi squared for the iterations, relative to regularized data:
	chi2_iter_thold = (chi2_target - chi2_reg_ends(0))/nchan
	
	chi2thold = 1.1*max(chi2_reg_ends)/nchan
	dempos_test = dblarr(nb2)
	datapos_test = dblarr(nb)
	
	; Loop over each pixel:
	for i=0,nx-1 do begin
		for j=0,ny-1 do begin

			itcount=0
			; Loop over regularization chi squareds:
			for k=0,nconfs-1 do begin
				chi2_reg_current = chi2_reg_ends(k)
				
				chi2_iter_thold = max([	chi2_target - chi2_reg_ends(0), $
										chi2_target - chi2_reg_ends(k)])/nchan
				
				; Regularize the data so that its chi squared (relative to the original data)
				; is equal to chi2_reg_current:
				data_out=reform(firdem_regularize_data(datain(i,j,*), errsin(i,j,*), a_struc, chi2_end=chi2_reg_current, niter_max=niter_reg_max));>0
				
				; Compute the first pass DEM corresponding to the regularized data:
				dem_initial = basis22#reform(a_struc.a_inv#(data_out/normfac))
				
				; Check to see if the first pass DEM needs iteration to remove negative
				; emission:
				dempos_test(*) = dem_initial > 0
				datapos_test(*) = a2_array#dempos_test
				deltatot = reform(((datapos_test-datain(i,j,*))/errsin(i,j,*))^2)
				chi2_current = total(deltatot/nchan)
				if(chi2_current le chi2thold) then begin
					chi2=chi2_current
					dem_out.coffs(i,j,*)=dempos_test
					break
				endif else begin	
					; Attempt to iterate away negative emission in the first pass DEM:
					firdem_iterate, dem_initial, errsin(i,j,*), a_struc, basis22, a2_array, dem, datapos=datapos, itcount=iters, chi2thold=chi2_iter_thold, datain=data_out, niter=niter_max, nbad_max=nbad_iter_max, mincond=minconds(k)

					itcount+=iters
				
					; Compute chi squared relative to the original data:
					deltatot = reform(((datapos-datain(i,j,*))/errsin(i,j,*))^2)
					chi2_current = total(deltatot/nchan)
					if(k eq 0) then chi2=chi2_current
					; Update the DEM if the new chi squared is better than the previous best:
					if(k eq 0 or chi2_current lt chi2) then begin
						chi2=chi2_current
						dem_out.coffs(i,j,*)=dem
					endif
					if(chi2 le chi2thold) then break ; Stop if chi squared is good enough
				endelse
			endfor
			
			chi2s(i,j) = chi2
			its(i,j) = itcount
			
		endfor				
	endfor

	print,"95th %ile chi squared:",estimate_quantile(chi2s,.95),", ", mean(its)," iterations average"

	telapsed = systime(1)-tstart
	print,nx*ny," DEMs computed, ", telapsed, " seconds elapsed."
	

end

