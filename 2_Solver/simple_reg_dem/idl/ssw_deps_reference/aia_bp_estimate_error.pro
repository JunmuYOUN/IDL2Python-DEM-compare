function AIA_BP_ESTIMATE_ERROR, counts, channel, $
	nocompress = nocompress, num_images = num_images, n_sample = n_sample, $
	darknoise = darknoise, quantnoise = quantnoise, readnoise = readnoise, $
	compressnoise = compressnoise, calibnoise = calibnoise, $
	chiantinoise = chiantinoise, $
	cal = cal, evenorm = evenorm, temperature = temperature, $
	loud = loud, pstop = pstop, noread = noread, force_read = force_read, $
	table = table

;+
;
;	$Log: aia_bp_estimate_error.pro,v $
;	Revision 1.7  2015/02/04 23:33:40  boerner
;	Improve dark estimate, compression for low pixel values, handling of summing
;	
;	Revision 1.6  2014/10/24 22:57:27  boerner
;	Use long loop counters
;	
;	Revision 1.5  2014/04/07 19:50:16  boerner
;	Added comments, fixed calibnoise keyword
;	
;	Revision 1.4  2012/02/10 23:56:29  boerner
;	Added keywords for specifying each component of error
;
;	Revision 1.3  2012/02/10 01:40:05  boerner
;	Update documentation, error handling, output table
;
;
;   Purpose: given an observed number of counts (DN/pixel), returns an estimate of the 
;		uncertainty (sigma) in the measurement (in DN/pixel)
;
;   Input Parameters:
;		counts		-	the measured counts (in DN/pixel); can be a scalar or array
;		channel		-	an integer specifying the wavelength of the AIA channel used to record
;						the observations; should have the same dimensions as counts
;
;	Keyword Parameters:
;		/loud		-	if set, a table indicating the size of the various error terms is printed
;		/pstop		-	set to halt inside the routine prior to returning (for debugging)
;		/noread		-	if set, then a default set of hard-coded error parameter estimates
;						are used (otherwise, the error parameters are loaded from a file
;						distributed with the SSWDB)
;		/force_read	-	if set, then the error values are re-read from disk even if they have been
;						previously read. By default, they are only read once in each IDL session,
;						and stored in a common block thereafter
;
;		n_sample	-	an integer specifying how many measurements (adjacent pixels, or consecutive
;						images) were averaged to produce the measured counts. Defaults to 1, 
;						i.e. no averaging. Correlated errors are reduced by sqrt(n_sample)
;		num_images	-	redundant with n_sample
;		/nocompress	-	if set, then the error due to on-board compression is not included
;
;		It is possible to include systematic errors that do not technically apply to the 
;		observations, but rather to the response functions that are generally used to predict
;		observations from a physical model. The following keywords pertain to these uncertainties:
;		/evenorm	-	if set, then the systematic error due to normalizing the effective
;						area using cross-calibration with SDO/EVE is included
;		/cal		-	if set, then the systematic error due to uncertainty in the preflight
;						photometric calibration is included (ignored if /evenorm is set; these
;						errors are larger than the errors associated with EVE normalization)
;		/temperature-	if set, then the systematic error due to uncertainty in the CHIANTI
;						data used to generate the temperature response function is included
;
; 	Output Parameters:
;		darknoise, quantnoise, readnoise, compressnoise, calibnoise, chiantinoise
;		There are all arrays of n_elements(counts) specifying the contribution
;		to the resulting error from various sources.
;		table 		-	A structure containing all of the above
;
;   Returns:
;		sigmadn		-	the uncertainty in the counts (in DN/pixel) (same dimensions as counts)
;
;-

common aia_bp_error_common, common_errtable

if N_ELEMENTS(n_sample) eq 0 then begin
	if N_ELEMENTS(num_images) eq 0 then begin
		n_sample = 1 
	endif else begin
		n_sample = num_images
	endelse
endif
sumfactor = SQRT(n_sample) / n_sample
numcounts = N_ELEMENTS(counts)
nullcount = FLTARR(numcounts)

okwaves = [94, 131, 171, 193, 211, 304, 335, 1600, 1700, 4500]
wave_ids = nullcount
if N_ELEMENTS(channel) eq 0 then channel = INTARR(numcounts)
for i = 0l, numcounts - 1 do wave_ids[i] = WHERE(okwaves eq channel[i])
if MIN(wave_ids) lt 0 then begin
	BOX_MESSAGE, ['AIA_BP_ESTIMATE_ERROR : Bad (or no) wavelength specified', $
		'Please set the second argument to one or more of:', STRJOIN(STRTRIM(okwaves, 2), '  '), $
		'Using a default set of error parameters...']
	noread = 1
endif

; Read in the error table
chanid = nullcount
if not KEYWORD_SET(noread) then begin
	if KEYWORD_SET(force_read) or N_ELEMENTS(common_errtable) eq 0 then begin
		errtable = AIA_BP_READ_ERROR_TABLE()
		common_errtable = errtable
	endif else begin
		errtable = common_errtable
	endelse
	for i = 0l, numcounts - 1 do chanid[i] = WHERE(errtable.wavelnth eq channel[i])
endif else begin	;	Create a dummy default error table right here
	errtable = {dnperpht : 1.122, calerr : 0.25, chianti : 0.25, compress : 14.4, eveerr : 0.02}
endelse

; Determine the typical number of photons per DN for the observed channel based on wavelength
; and camera gain
if N_ELEMENTS(dnperpht) eq 0 then dnperpht = errtable[chanid].dnperpht

if N_ELEMENTS(shotnoise) eq 0 then begin
	numphot = counts / dnperpht
	stdevphot = SQRT(numphot)
	shotnoise = (stdevphot * dnperpht)
	shotnoise = nullcount + shotnoise / SQRT(n_sample)		;	DN uncertainty due to shot noise
endif

;	DN uncertainty due to dark subtraction 
if N_ELEMENTS(darknoise) eq 0 then darknoise = nullcount + 0.18

;	DN uncertainty due to read noise
if N_ELEMENTS(readnoise) eq 0 then readnoise = nullcount + 1.15 * sumfactor		

;	DN uncertainty due to quantization
if N_ELEMENTS(quantnoise) eq 0 then quantnoise = nullcount + 0.288819 * sumfactor

;	DN uncertainty due to onboard compression
if N_ELEMENTS(compressnoise) eq 0 then begin
	if not KEYWORD_SET(nocompress) then begin
		compressratio = errtable[chanid].compress
		compressnoise = (shotnoise / compressratio) > 0.288819
		lowcounts = WHERE(counts lt 25, numlow)	;	Linear portion of lookup table
		if numlow gt 0 then compressnoise[lowcounts] = 0.
		compressnoise = compressnoise * sumfactor
	endif else begin
		compressnoise = nullcount
	endelse
endif

;	DN uncertainty due to errors in photometric calibration
case 1 of						;	% uncertainty due to errors in photometric calibration
	KEYWORD_SET(evenorm)	:	calerr = errtable[chanid].eveerr		
	KEYWORD_SET(cal)		:	calerr = errtable[chanid].calerr	
	else					:	calerr = nullcount
endcase
if N_ELEMENTS(calibnoise) eq 0 then calibnoise = calerr * counts	

;	DN uncertainty due to errors in CHIANTI model
if KEYWORD_SET(temperature) then begin
	chiantierr = errtable[chanid].chianti	;	% uncertainty due to errors in CHIANTI model
endif else begin
	chiantierr = nullcount
endelse
if N_ELEMENTS(chiantinoise) eq 0 then chiantinoise = chiantierr * counts			

sigmadn = SQRT(shotnoise^2. + darknoise^2. + readnoise^2. + quantnoise^2. + compressnoise^2. + $
	chiantinoise^2. + calibnoise^2.)
snr = counts / sigmadn

if KEYWORD_SET(loud) then begin
	PRINT
	PRINT, 'Counts [DN]', 'RMS Error', 'SNR', 'Shot', 'Dark', 'Read', 'Quant', 'Compress', 'Chianti', $
		'Calibrat', format = '(2a12, 8a9)'
	for i = 0l, N_ELEMENTS(counts)-1 do PRINT, counts[i], sigmadn[i], snr[i], shotnoise[i], darknoise[i], $
		readnoise[i], quantnoise[i], compressnoise[i], chiantinoise[i], calibnoise[i], $
		format = '(2f12.2, 8f9.2)'
	PRINT
endif

if KEYWORD_SET(pstop) then STOP

table = CREATE_STRUCT('darknoise', darknoise, 'quantnoise', quantnoise, $
	'readnoise', readnoise, 'compressnoise', compressnoise, 'calibnoise', calibnoise, $
	'chiantinoise', chiantinoise)

RETURN, sigmadn

end
