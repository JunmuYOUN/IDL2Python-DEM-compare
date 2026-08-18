FUNCTION PREP_VDEM_SUMER,index1,data1,WAVE=wave,BACKGROUND=background, $
		VERBOSE=verbose,wbins=wbins,ybins=ybins,$
		num_bck=num_bck,percent_bck=percent_bck,$
		arraynum=arraynum,specs=specs,exptime=exptime

;
; NAME:
;
;	prep_vdem_sumer.pro
;
; PURPOSE:
;
;	Part of a package of routines associated with the calculation
;	of the Velocity Differential Emission Measure (VDEM).
;	This routine extracts data from SUMER data array, totals certain pixels,
;	calibrates those pixels, subtracts the background, and converts
;	specific intensity to flux. The output of this routine becomes the
;	input in the VDEM.pro routine.
;
; CATEGORY:
;
;       Data analysis
;
; CALLING SEQUENCE:
;
;	flux_parameters=PREP_VDEM_SUMER(index,data,[WAVE=wave,$
;			BACKGROUND=background,wbin=wbin,ybins=ybins,$
;			num_bck=num_bck,percent_bck=percent_bck,$
;			arraynum=arraynum,specs=specs,exptime=exptime,$
;			/verbose])
;
; INPUTS:
;
;	index=index array associated with SUMER data set (from rd_sumer.pro)
;	data= data array associated with SUMER data set (from rd_sumer.pro)
;
; OPTIONAL INPUTS
;
;	wave= wavelength array calculated from reference spectra or some
;		other method (if not defined, sgt_wave.pro will be used)
;	background=constant or array with length as bins to be analyzed
;		for background subtraction.
;       specs=array containing [minimum spectrum number to be analyzed,
;               maximum spectrum number to be analyzed]
;               All spectra from the minimum number to the maximum
;               number will be analyzed and returned.  Default is all
;               spectra in the data array.
;	exptime= exposure time of the spectra being analyzed.  (In general,
;		the exposure time can be found from the index structure, and
;		this input in not required.  Occassionally, however,
;		the index structure contains an exposure time of 0 or an
;		incorrect exposure time.)
;
;	(If the following optional inputs are not provided in the call
;	sequence, the user will be prompted for the information.)
;
;	arraynum= an array containing [data array number]. If more than
;		one data array is in the SUMER data structure, all arrays
;		will be printed and user will be asked for the array
;               number.  (Note: array number of first array is
;               0, second array is 1 and so on.)
;	wbins= array containing [minimum horizontal pixel number,
;		maximum horizontal number] to be analyzed.
;	ybins= array containing [minimum vertical pixel number,
;		maximum vertical pixel number] to be analyzed.
;	num_bck= number of pixels to be used in background subtraction
;		(suggest 5).  This number of pixels on either
;		side of the line will be fitted with a first degree
;		degree polynomial.  [If a background is supplied (see
;		background optional input above), this input will not be used.]
;	percent_bck= percentage of background to subtract (suggest .7).
;		[If a background is supplied (see background
;		optional input above), this input will not be used.]
;
; KEYWORDS:
;
;	/verbose= if set, user will be told the process as the function
;		calculates the flux parameters.
;
; OUTPUTS:
;
;	flux_parameters= structure containing the flux, error in flux,
;		and wavelength array to be used as input into vdem.pro.
;
; OPTIONAL OUTPUTS:
;
;	None.
;
; COMMON BLOCKS:
;
;	None.
;
; SIDE EFFECTS:
;
;	None.
;
; RESTRICTIONS:
;
;	Original SUMER index and data array from rd_sumer.pro must be
;	the inputs.  (Data array must be in counts and should not have been
;       changed by the count2countrate or sgt_calibration (/apply keyword)
;       programs.  If the data array HAS been changed by these programs,
;       you will receive an error message.)
;
;	If the index array incorrectly contains an exposure time of 0, an
; 	error message will be printed.
;
;	Note that routine converts SUMER's specific intensity to flux.
;
; PROCEDURE:
;
;	Find index and data structures using rd_sumer routine.  Suggest
;	also calculating wavelength array using reference lines.  Locate
;	an event and determine the:
;		a) data array number
;		b) minimum and maximum vertical pixels
;		c) minimum and maximum horizontal pixels
;		d) number of pixels to use to determine bachground -or-
;			constant background to subtract
;		e) percentage of background to subtract
;	for entry into the program when prompted or as optional inputs.
;
; EXAMPLE:
;
;	Read in some data using rd_sumer:
;
; IDL> file='/data1/soho/sumer/sumer_fits/amy/sum_961102_132633.genx'
; IDL> rd_sumer,file,index,data
;
;	After finding all necessary information (see PROCEDURE), run
;	prep_vdem_sumer.
;
; IDL> flux_parameters=prep_vdem_sumer(index,data,/verbose)
;
; There is only one data array in data structure.
; Choosing that array.
;
; Setting wavelength array using sgt_wave or supplied wavelength array.
;
;
; Calculating wavelength array using sgt_wave.
;
;
; Getting calibration array
;
;
; Finding date from SUMER filename and
; changing data from specific intensity to flux.
;
; Input the minimum and maximum vertical pixel numbers to be totaled
; Minumum:52
; Maximum:58
;
; Input the minumum and maximum horizontal pixel numbers.
; Make the window large enough to capture the entire event,
; but do not include excess on either side of the spectral line.
;
; Minimum horizontal pixel number:195
; Maximum horizontal pixel number:245
;
; Input number of pixels to use on either side of line to
; calculate background (suggest 5).
;
; Number of pixels:5
;
; Input percentage of background to subtract (suggest 0.7)
;
; Percentage of background: 0.7
;
; IDL> help,flux_parameters,/str
; ** Structure <422650>, 3 tags, length=4692, refs=1:
;    FLUX            FLOAT     Array(51, 11)
;    EFLUX           FLOAT     Array(51, 11)
;    WAVE            FLOAT     Array(51)
; IDL>
;
;	To plot the flux of nth spectrum as a function of wavelength
;	with errors:
;
; IDL> ploterr,flux_parameters.wave,flux_parameters.flux(*,n),$
; IDL> flux_parameters.eflux(*,n)
;
; MODIFICATIONS:
;
;	A. Winebarger 3/98- Original Code
;
index=index1
data=data1


au=1.495979e13          ;astronomical unit (cm)
rr=6.9626e10		;solar radius (cm)
nmax=N_tags(data)

;;******************************************************************
;;
;;      Get data array number if not provided by user or 0
;;
;;******************************************************************

IF nmax EQ 1 THEN BEGIN
arraynum=0
	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,''
		PRINT,'There is only one data array in data structure.'
		PRINT,'Choosing that array.'
		PRINT,''
	ENDIF
ENDIF ELSE IF KEYWORD_SET(arraynum) THEN arraynum=arraynum(0) ELSE BEGIN
help,data,/str
arraynum=0
PRINT,''
PRINT,'Input data array to be analyzed.'
PRINT,'Note that first data array is number 0, second is 1, and so on.'
READ,arraynum,PROMPT="Data array number:"
ENDELSE
;;******************************************************************
;;
;;      Make sure data has not been changed by sgt_calibration or
;;      counts2countrate
;;
;;******************************************************************

datacheck_c=strpos(index.(arraynum).spectrum.history(0),'counts2countrate')
datacheck_s=strpos(index.(arraynum).spectrum.history(0),$
        'radiometric calibration')
FOR i=0, n_elements(datacheck_c)-1 DO BEGIN
IF (datacheck_c(i) NE -1) OR (datacheck_s(i) NE -1) THEN BEGIN
PRINT,'%% ERROR:  DATA must not be altered by counts2countrate or '
PRINT,'%% radiometric calibration.  No flux returned.'
        GOTO, DONE
ENDIF
ENDFOR

;;******************************************************************
;;
;;      Get wave array from user or sgt_wave
;;
;;******************************************************************



IF KEYWORD_SET(WAVE) THEN BEGIN
	wave=wave
	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,''
		PRINT,'Setting wavelength array to user-defined values.'
		PRINT,''
	ENDIF
ENDIF ELSE BEGIN
wave=sgt_wave(index)
wave=wave(*,arraynum)
	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,''
		PRINT,'Calculating wavelength array using sgt_wave.'
		PRINT,''
	ENDIF
ENDELSE

;;******************************************************************
;;
;;      Get specs number if specs keyword is set.  Otherwise
;;      finding flux for all specs in data set.
;;
;;******************************************************************

IF keyword_set(specs) THEN BEGIN
        specmin=specs(0)
        specmax=specs(1)
ENDIF ELSE BEGIN
        specmin=0
        specmax=n_elements(data.(arraynum)(0,0,*))-1
ENDELSE

;;******************************************************************
;;
;;      Get calibration array from sgt_calibration.  Default
;;      units are photons, unless is keyword /ergs is set.
;;
;;******************************************************************


	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,''
	PRINT,'Getting calibration array'
		PRINT,''
	ENDIF
radcal=sgt_calibration(index,data,/photons)
radcal=radcal(*,arraynum)

;;******************************************************************
;;
;;	Get exposure time and check to make sure it is not 0.
;;
;;******************************************************************

IF KEYWORD_SET(exptime) THEN $
        exptime=fltarr(n_elements(data.(arraynum)(0,0,*)))+exptime $
ELSE exptime=sgt_exptime(index.(arraynum))      ;to change counts to countrate
FOR i=0,n_elements(exptime)-1 DO IF exptime(i) EQ 0 THEN BEGIN
	print,'%% ERROR: Exposure time cannot be 0.'
	print,'%% Possible error in index structure,'
	print,'%% suggest using EXPTIME optional input.'
	print,'%% No flux returned.'
	ENDIF

;;******************************************************************
;;
;;      Find data from filename, use to find radius in arcsecs.
;;
;;******************************************************************


	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,''
		PRINT,'Finding date from SUMER filename and'
		PRINT,'changing data from specific intensity to flux.'
		PRINT,''
	ENDIF
year=strmid(index.gen.filename,4,2)
month=strmid(index.gen.filename,6,2)
	IF month EQ '01' THEN month='jan'
	IF month EQ '02' THEN month='feb'
	IF month EQ '03' THEN month='mar'
	IF month EQ '04' THEN month='apr'
	IF month EQ '05' THEN month='may'
	IF month EQ '06' THEN month='jun'
	IF month EQ '07' THEN month='jul'
	IF month EQ '08' THEN month='aug'
	IF month EQ '09' THEN month='sep'
	IF month EQ '10' THEN month='oct'
	IF month EQ '11' THEN month='nov'
	IF month EQ '12' THEN month='dec'
day=strmid(index.gen.filename,8,2)

date=string(day,'-',month,'-',year)

radius_arcsecs=get_rb0p(date,/radius)

;;******************************************************************
;;
;;      Get vertical bin number from user
;;
;;******************************************************************

IF KEYWORD_SET(ybins) THEN BEGIN
	ymin=ybins(0)
	ymax=ybins(1)
ENDIF ELSE BEGIN
ymin=0
ymax=0
PRINT,''
PRINT,'Input the minimum and maximum vertical pixel numbers to be averaged.'
READ,ymin,PROMPT="Minumum:"
READ,ymax,PROMPT="Maximum:"
ENDELSE

;;******************************************************************
;;
;;      Find flux [photons/cm^2/s/A] from specific intensity
;;      [photons/cm^2/s/A/str] by multiplying by the area of
;;      the (2-dimnesional) sun in one pixel/sun_sat_distance^2.
;;      Specific intensity profile is found by averaging the pixels
;;      along the slit (ymin:ymax), multiplying by radcal, and
;;      dividing by exptime.
;;      Error in counts is assumed to be Poisson [error=sqrt(counts)].
;;      Error in spec. intensity and flux is found by multiplying
;;      error in counts by same thing as counts.
;;
;;******************************************************************

nflux=n_elements(data.(arraynum)(*,0,0))
nspec=specmax-specmin+1
nspat=ymax-ymin+1

spec_intensity=dblarr(nflux,nspat,nspec)
espec_intensity=dblarr(nflux,nspat,nspec)

tempflux=dblarr(nflux,nspec)
tempeflux=dblarr(nflux,nspec)

FOR p=0,nspec-1 DO FOR i=0, nspat-1 DO BEGIN
        spec_intensity(*,i,p)=data.(arraynum)(*,ymin+i,specmin+p)*$
                        (radcal/exptime(specmin+p))
        espec_intensity(*,i,p)=sqrt(data.(arraynum)(*,ymin+i,specmin+p))*$
                        (radcal/exptime(specmin+p))
ENDFOR

slit_n=sgt_slit(index)
IF (slit_n EQ 1) THEN slit_dim=4.*300. ELSE $
IF (slit_n EQ 2) THEN slit_dim=1.*300. ELSE $
IF (slit_n EQ 4 OR slit_n EQ 3 OR slit_n EQ 5) THEN slit_dim=1.*120. ELSE $
IF (slit_n EQ 7 OR slit_n EQ 6 OR slit_n EQ 8) THEN slit_dim=.3*120. ELSE BEGIN
        Print,'%% ERROR:  Slit number is undetermined'
        PRINT,'%% No flux returned.'
        GOTO,DONE
ENDELSE

darea=(slit_dim/n_elements(data.(arraynum)(0,*,0)))*(rr/radius_arcsecs)^2


FOR p=0,nspec-1 DO BEGIN
        IF (nspat EQ 1) THEN BEGIN
tempflux(*,p)=spec_intensity(*,*,p)*darea/(.99*au)^2
                                                ;photons/s/cm^2/A
tempeflux(*,p)=espec_intensity(*,*,p)*darea/(.99*au)^2
                                                ;photons/s/cm^2/A
ENDIF ELSE BEGIN

tempflux(*,p)=average(spec_intensity(*,*,p)*darea/(.99*au)^2,2)
                                                ;photons/s/cm^2/A
tempeflux(*,p)=(darea/(.99*au)^2)*sqrt(total(espec_intensity(*,*,p)^2,2))/$
        float(ymax-ymin+1)
ENDELSE
                                                ;photons/s/cm^2/A
ENDFOR

;;******************************************************************
;;
;;      Get minimum and maximum bin number (along lambda to include
;;      line of interest)
;;
;;******************************************************************

IF KEYWORD_SET(wbins) THEN BEGIN
	wfirst=wbins(0)
	wlast=wbins(1)
ENDIF ELSE BEGIN
wfirst=0
wlast=0
PRINT,''
PRINT,'Input the minumum and maximum horizontal pixel numbers.'
PRINT,'Make the window large enough to capture the entire event,'
PRINT,'but do not include excess on either side of the spectral line.'
PRINT,''
READ,wfirst,PROMPT="Minimum horizontal pixel number:"
READ,wlast,PROMPT="Maximum horizontal pixel number:"
ENDELSE

;;******************************************************************
;;
;;      Calculate and subtract background
;;
;;******************************************************************

flux=dblarr(wlast-wfirst+1,nspec)
eflux=dblarr(wlast-wfirst+1,nspec)

IF KEYWORD_SET(BACKGROUND) THEN BEGIN
	IF KEYWORD_SET(VERBOSE) THEN BEGIN
		PRINT,''
		PRINT,'Subtracting user-supplied background'
		PRINT,''
	ENDIF
bck=fltarr(nflux)
bck(*)=background
for i=0,nspec-1 do flux(*,i)=flux(*,i)-bck
ENDIF ELSE BEGIN

IF KEYWORD_SET(num_bck) THEN num_bck=num_bck ELSE BEGIN
PRINT,''
PRINT,'Input number of pixels to use on either side of line to '
PRINT,'calculate background (suggest 5).'
PRINT,''
num_bck=5
READ,num_bck,PROMPT="Number of pixels:"
ENDELSE

IF KEYWORD_SET(percent_bck) THEN percent_bck=percent_bck ELSE BEGIN
PRINT,''
PRINT,'Input percentage of background to subtract (suggest 0.7).'
PRINT,''
percent_bck=0.7
READ,percent_bck,PROMPT="Percentage of background:"
ENDELSE

x=intarr(2*num_bck)
FOR j=0,num_bck-1 DO x(j)=wfirst-num_bck+j
FOR j=num_bck,2*num_bck-1 DO x(j)=wlast+(j-num_bck+1)

FOR i=0,nspec-1 DO BEGIN

	y=tempflux(x,i)
	l=poly_fit(x,y,1)
	bck=(l(0)+findgen(nflux)*l(1))*percent_bck
	tempflux(*,i)=tempflux(*,i)-bck
	FOR p=0, nflux-1 DO IF tempflux(p,i) LE 0. THEN tempflux(p,i)=1.
	FOR p=0, nflux-1 DO IF tempeflux(p,i) LE 0. THEN tempeflux(p,i)=1.

ENDFOR
ENDELSE

flux=tempflux(wfirst:wlast,*)
eflux=tempeflux(wfirst:wlast,*)

flux_parameters={flux:flux,eflux:eflux, wave:wave(wfirst:wlast)}

RETURN,flux_parameters

DONE:

END
