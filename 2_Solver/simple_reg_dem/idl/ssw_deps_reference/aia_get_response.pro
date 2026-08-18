function AIA_GET_RESPONSE, $
	effective_area=effective_area, area=area, temperature=temperature, emissivity = emissivity, $
	full=full, all=all, uv=uv, dn=dn, phot=phot, noblend=noblend, $
	use_photospheric = use_photospheric, $
	evenorm=evenorm, timedepend_date=timedepend_date, chiantifix=chiantifix, $
	version=version, emversion = emversion, respversion = respversion, $
	silent=silent, loud=loud
   
;+
;
;	Name: aia_get_response
;
;	Purpose: return SDO/AIA instrument response data structures
;
;	Input Parameters:
;
;	Output:
;		function returns AIA wavelength response, spectral model, or temperature response
;		Descriptions of the structures returned by this routine are available via: 
;			http://sohowww.nascom.nasa.gov/solarsoft/sdo/aia/response/README.txt
;
;	Calling Examples:
;		IDL> effarea=aia_get_response(/area, /dn)		;	return per wavelength effective area 
;		IDL> effarea=aia_get_response(/area,/full,/dn)	;	same with component details (filter/ccd...)
;		IDL> emiss=aia_get_response(/emiss, /full)		;	return default CHIANTI model with line list
;		IDL> tresp=aia_get_response(/temp,/dn,/evenorm)	;	return temperature response including
;														; 	constant scale factors to give good
;														;	overall agreement with SDO/EVE
;
;	Keyword Parameters:
;
;		SPECIFY ONE OF THE FOLLOWING TO DETERMINE THE TYPE OF RESPONSE STRUCTURE:
;
;		area - returns AIA instrument response/effective area data
;		emissivity - returns the default CHIANTI emissivity model
;		temperature - returns AIA temperature response (= effective area -> CHIANTI model)
;			NOTE that the temperature response is generated on-the-fly,
;			while the area and emissivity are loaded from dbase files
;		effective_area - synonym for area 
;		(If none of these is specified, an effective area structure is used)
;
; 		DEPENDING ON WHICH TYPE OF RESPONSE STRUCTURE IS SELECTED, some or all
;		of the following keyword switches can be set:
;
;		dn	- (AREA and TEMP ONLY) if set, the response function includes the 
;			photon-to-DN units conversion 
;		phot - (AREA and TEMP ONLY) if set, the response function does not include the
;			photon-to-DN units conversion
;		(NOTE that either /DN or /PHOT should be set if you are requesting an AREA or TEMP)
;		full - if set, return the detailed version of the data.
;			for /area, this includes component-level efficiency measurements
;			for /emiss, this includes the full CHIANTI line list (~100s of MB)
;			for /temp, this includes the grid of wavelength and temperature response
;		all - if set, return the response for all EUV channels (including thin and thick
;			focal plane filters). By default, only return the response for the thin filters
;		noblend - (AREA and TEMP ONLY) by default (if this keyword is NOT set), the 
;			telescope 1 and telescope 4 mirror reflectivities are combined to account 
;			for possible cross-contamination in the channels on those ATAs. This 
;			affects the 94, 131 and (especially) 304 and 335 channels. Set this keyword
;			to consider each channel individually.
;		uv	- (AREA ONLY) if set, return response for the UV channels (1600/1700/4500).
;			By default, return the EUV channels (94/131/171/193/211/304/335). If /uv is
;			set, then /all and /noblend are ignored (/full and /dn are still available 
;			as options)
;		use_photospheric - (EMISS and TEMP ONLY) if set, then photospheric abundances
;			are used; by default, coronal abundances are used.
;
;		VARIOUS CORRECTIONS TO THE RESPONSE FUNCTION CAN BE APPLIED BASED ON ANALYSIS OF THE
;		INSTRUMENT CALIBRATION (DERIVED FROM CROSS-CALIBRATION WITH OTHER INSTRUMENTS)
;		Note that these corrections are included in the response structure, whether or not they
;		are applied, along with a flag indicating whether they were applied.
;			evenorm	-	(AREA and TEMP ONLY) if set, then the response functions are normalized
;				to give good agreement with full-disk EVE observations
;			timedepend_date	-	(AREA and TEMP ONLY) can be set as a switch (/timedepend_date) 
;				or set to a string specifying a date in any format recognized by ANYTIM.PRO.
;				If set, then the response function is calculated for the specified date.
;				(If set as a switch, the date at which the program is being called is used)
;				By default, the response at the start of the SDO mission (1-May-2010) is returned
;				(i.e. by default there is no correction for degradation on orbit)
;				Note that, if you want time dependent correction but do not want EVE normalization,
;				You should explicitly set evenorm = 0
;			chiantifix	-	(TEMP ONLY) if set, then an empirical correction is applied to the 
;				temperature response of the 94 and 131 A channels in order to account for emission
;				not in the CHIANTI database. For more detail on this correction, see
;				http://sohowww.nascom.nasa.gov/solarsoft/sdo/aia/response/chiantifix_notes.txt
;	
;		silent	-	if set, do not print any diagnostic information while reading
;			or generating the response functions
;		loud	-	if set, print, some extra diagnostic information
;		version	-	set to the integer version number of the response functions to use. 
;			Defaults to the	latest version (currently 6)
;		emversion	-	(TEMP ONLY) If requesting a temperature response, this keyword gives you the 
;			option of using a different version number for the emissivity and the instrumental
;			response characteristics (i.e. you can set version=2 and emversion=1 to see
;			what the T response looks like with the revised instrument response but pre-revision
;			emissivity)
;		respversion	-	The response table that tracks the time-dependent degradation by cross-
;			calibration with EVE is updated periodically without incrementing the version number
;			of the response function. By default, the latest response table is retrieved; set this
;			keyword to a string to force the routine to use a deprecated response table for calculating
;			time-dependent corrections
;
;	History:
;		15-Jun-2010 - Calibration by Paul Boerner boerner@lmsal.com; trivial SSW wrapper by S.L.Freeland	 
;		21-Jun-2010	- Updated with preflight calibration, and options for /uv, /dn, /all
;		30-Jan-2012	- Updated with V2 calibration, including corrections for time dependence, 
;						EVE normalization, missing CHIANTI lines, etc.
;		28-Sep-2012	-	Updated with V3 calibration
;		7-Jan-2013	-	Updated with V4 calibration
;		18-Mar-2014	-	Updated with V6 calibration
;
;		2017/11/30	- 	1) updated version handling, follow the revised aia_bp_get_corrections.pro
;						now include new V8 of the response table, while the other files remain the same as V6
;		2017/12/01, Fri: 1) take 7 out of allowedversions, it's just an intermediate version made before Paul left, but never released to SSW
;		2020/07/13  -   Updated with V9 calibration
;		2020/10/22  -   Updated with V10 calibration
;
;	$Id: aia_get_response.pro,v 1.28 2015/06/03 21:41:32 boerner Exp $
;
;	Contact:
;		Questions, problems or suggestions, email 
;			boerner ~at~ lmsal ~dot~ com
;
;-

aiaresp=GET_LOGENV('AIA_RESPONSE_DATA')
if aiaresp eq '' then aiaresp=CONCAT_DIR('$SSW_AIA','response')
aiaresp = (FILE_SEARCH(aiaresp))[0]

area=KEYWORD_SET(area) or KEYWORD_SET(effective_area)
temp=KEYWORD_SET(temperature)
emiss=KEYWORD_SET(emissivity)
full=KEYWORD_SET(full)
uv=KEYWORD_SET(uv)
all=KEYWORD_SET(all)
dn=KEYWORD_SET(dn)
phot=KEYWORD_SET(phot)
noblend=KEYWORD_SET(noblend)
silent = KEYWORD_SET(silent)
loud = KEYWORD_SET(loud)
if N_ELEMENTS(evenorm) gt 0 then user_evenorm = 1-evenorm else user_evenorm = 0 ;	Did the user explicitly set it to evenorm=0?
evenorm = KEYWORD_SET(evenorm)
timedepend = KEYWORD_SET(timedepend_date)
chiantifix = KEYWORD_SET(chiantifix)

;++++++++++++++++++++++++++
; Detect and warn of 
; incompatible keyword 
; settings
;--------------------------
if (temp + emiss + area) gt 1 then begin
	BOX_MESSAGE, ['Please specify only one of the following:', $
		'  /temperature  /emissivity  /area', $
		'Defaulting to /area...']
	temp = 0
	emiss = 0
	area = 1
endif	
if (temp + emiss + area) eq 0 then area = 1	;	Silently switch on /area
if uv then begin
	area = 0			;	Henceforth, (area=1) implies EUV area
	all = 1				; 	Always return all UV channels
endif
;if temp and uv then begin
;	BOX_MESSAGE, ['Cannot generate temperature response for UV channels', $
;		'Returning UV channel effective area']
;	temp = 0
;endif
;if emiss and uv then begin
;	BOX_MESSAGE, ['No emissivity defined for UV channels', $
;		'Returning UV channel effective area']
;	emiss = 0
;endif
if emiss and (dn or noblend or phot or timedepend or evenorm) then begin
	BOX_MESSAGE, ['/DN, /PHOT, /TIME, /EVENORM and /NOBLEND keywords are not defined for emissivity structure']
	dn = 0
	noblend = 0
	phot = 0
	timedepend = 0
	evenorm = 0
endif
if ((dn + phot) ne 1) and (not emiss) then begin	;	Units not explicitly specified - check with user
	BOX_MESSAGE, ['AIA_GET_RESPONSE...', $
				'Do you want to include the DN/phot factor in the response function?', $
				'Specify either /dn or /phot keyword to avoid this message']
	reply = ''
	READ, 'Include DN/phot correction? (Default is yes) [y/n] ', reply
	if STRUPCASE(reply) eq 'N' then begin
		BOX_MESSAGE, 'AIA_GET_RESPONSE: NOT including DN/phot factor'
		dn = 0
		phot = 1
	endif else begin
		BOX_MESSAGE, 'AIA_GET_RESPONSE: including DN/phot factor'
		dn = 1
		phot = 0
	endelse
endif
if timedepend and not evenorm and not user_evenorm then begin	;	Timedepend but not evenorm - check if user is sure
	BOX_MESSAGE, ['AIA_GET_RESPONSE...', $
				'You asked for time dependent correction but not EVE normalization;', $
				'Are you sure you do not want to include EVE normalization?']
	reply = ''
	READ, 'Include EVE normalization as well? (Default is no) [y/n] ', reply
	if STRUPCASE(reply) eq 'Y' then begin
		BOX_MESSAGE, 'AIA_GET_RESPONSE: Including EVE normalization'
		evenorm = 1
	endif else begin
		BOX_MESSAGE, 'AIA_GET_RESPONSE: Not including EVE normalization'
	endelse
endif

;++++++++++++++++++++++++++
; Determine which version
; of the SSWDB files should
; be used
;--------------------------
	;--- 2017/11/30, replaced below, follow revised aia_bp_get_corrections.pro ---------
	;defaultversion = 6
	;allowedversions = [1,2,3,4,6]
allowedversions = [1,2,3,4,6,8,9,10]		;[1,2,3,4,6,7,8]
defaultversion = max(allowedversions)

if N_ELEMENTS(version) eq 0 then version = defaultversion
okversion = WHERE(FIX(version) eq allowedversions, isok)
if not isok then  begin
	print, 'Warning: user-specified version ' + trim(version) + ' NOT allowed, now force version to be default of '  + trim(defaultversion)
	version = defaultversion	;	silently set the version to the default if a funny value is used
endif

if version eq 1 then begin
	vstring = 'aia_preflight_'
endif else begin
	vstring = 'aia_V' + trim(version) + '_'	
endelse
; Note that, for v10, only the response table file is different
tablevstring = vstring
if version eq 10 then vstring = 'aia_V9_'

if KEYWORD_SET(use_photospheric) then abundtag = 'photo_' else abundtag = ''
if KEYWORD_SET(uv) then abundtag = 'uv_' + abundtag

;++++++++++++++++++++++++++
; More keyword consistency
; checking (look at chianti
; fix keyword)
;--------------------------
if chiantifix then begin
	if temp then begin						;	 chiantifix only applies if /temp is set
		if not evenorm then begin
			BOX_MESSAGE, ['/chiantifix requires /evenorm - using EVE-normalized response']
			evenorm = 1
		endif
		if version lt 2 then begin
			BOX_MESSAGE, ['/chiantifix requires version > 1', 'Disabling chiantifix']
			chiantifix = 0
		endif
		if chiantifix then begin
			RESTGEN, file = aiaresp + '/' + vstring + 'chiantifix.genx', str = ch_str
		endif
	endif else begin
		BOX_MESSAGE, ['/chiantifix only applies for temperature response']
	endelse
endif


;++++++++++++++++++++++++++
; Check possible cases and 
; return the appropriate 
; structure
;--------------------------

; Emissivity
if emiss then begin
	filename = CONCAT_DIR(aiaresp, vstring) + abundtag + 'fullemiss.genx'
	if FILE_EXIST(filename) then begin
		RESTGEN,file=filename,struct = data
	endif else begin
		BOX_MESSAGE, ['Cannot find AIA emissivity file', filename]
		RETURN, filename
	endelse
	if full then begin
		output = data
	endif else begin
		output = CREATE_STRUCT(data.total, 'CH_Info', data.general)
	endelse
	output = CREATE_STRUCT(output, 'Version', FIX(version))
	RETURN, output
endif else begin
	;++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	; Everything except emissivity may be time dependent, 
	; so read in the data needed to perform the calculation
	; of the time dependence in the instrument response
	;--------------------------------------------------------
	tablefilebase = CONCAT_DIR(aiaresp, tablevstring)
	respfiles = FILE_SEARCH(tablefilebase + '*_response_table.txt')
	for i = 0, N_ELEMENTS(respfiles)-1 do $
		respfiles[i] = STRMID(respfiles[i], STRLEN(tablefilebase), 15)
	respfile_tais = ANYTIM2TAI(FILE2TIME(respfiles))
	if not KEYWORD_SET(respversion) then respversion = respfiles[WHERE(respfile_tais eq MAX(respfile_tais))]
	if not FILE_EXIST(tablefilebase + respversion + '*') then begin
		BOX_MESSAGE, ['File not found: ', tablefilebase + respversion + '*']
		respversion = respfiles[WHERE(respfile_tais eq MAX(respfile_tais))]
	endif
	tablefilebase = tablefilebase + respversion + '_*'
	ver_date = respversion
	tablefile = tablefilebase + 'response_table.txt'
	tabledat = AIA_BP_READ_RESPONSE_TABLE(tablefile, silent = silent)
	sampledate = RELTIME(/now)
	if timedepend then begin
		if timedepend_date ne '1' then sampledate = timedepend_date
	endif
	sampleutc = ANYTIM2UTC(sampledate, /ccsds)
endelse


; UV effective area
if uv then begin
	filename = CONCAT_DIR(aiaresp, vstring) + 'fuv_fullinst.genx'
	if FILE_EXIST(filename) then begin
		RESTGEN,file=filename,struct = data
	endif else begin
		BOX_MESSAGE, ['Cannot find AIA UV area file', filename]
		RETURN, filename
	endelse
	full_resp = AIA_BP_PARSE_EFFAREA(data, short_resp, dn = dn, $
		version = version, sampleutc = sampleutc, resptable = tabledat, $
		evenorm = evenorm, timedepend = timedepend, /uv, ver_date = ver_date)
	if full then result = full_resp else result = short_resp
	RETURN, result
endif


; Temperature Response
if temp then begin
	; Check to see if the user has requested a different version of the emissivity than
	; what is being used for the instrument response
	if KEYWORD_SET(emversion) then begin
		case emversion of
		1	:	evstring = 'aia_preflight_'
		else:	evstring = 'aia_V' + trim(emversion) + '_'		;2017/11/30

		;2	:	evstring = 'aia_V2_'
		;3	:	evstring = 'aia_V3_'
		;4	:	evstring = 'aia_V4_'
		;5	:	evstring = 'aia_V5_'
		;6	:	evstring = 'aia_V6_'
		endcase
	endif else begin
		emversion = version
		evstring = vstring
	endelse
	
	afilename = CONCAT_DIR(aiaresp, vstring) + 'all_fullinst.genx'
	efilename = CONCAT_DIR(aiaresp, evstring) + abundtag + 'fullemiss.genx'
	
	filemissing	= not( FILE_EXIST(afilename) and FILE_EXIST(efilename) )
	if filemissing then begin
		BOX_MESSAGE, ['Cannot find area or emissivity file:', afilename, efilename]
		RETURN, afilename + ',' + efilename
	endif else begin
		if not silent then BOX_MESSAGE, ['Generating temperature response function from', $
			afilename, efilename]
	endelse

	RESTGEN, file = afilename, str = areadat
	if not noblend then areadat = AIA_BP_BLEND_CHANNELS(areadat)
	if not all then chanlist = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
	full_area = AIA_BP_PARSE_EFFAREA(areadat, short_area, dn = dn, channels = chanlist, $
		version = version, sampleutc = sampleutc, resptable = tabledat, $
		evenorm = evenorm, timedepend = timedepend, ver_date = ver_date)	
	RESTGEN, file = efilename, str = emiss_str

	t_short_resp = AIA_BP_AREA2TRESP(short_area, emiss_str, t_full_resp, emversion)
	full_resp = AIA_BP_PARSE_TRESP(t_full_resp, short_resp, chiantifix = ch_str)
	if full then result = full_resp else result = short_resp
	RETURN, result
endif


; If you've gotten this far, you must want the EUV effective area
if area then begin
	filename = CONCAT_DIR(aiaresp, vstring) + 'all_fullinst.genx'
	if FILE_EXIST(filename) then begin
		RESTGEN, file=filename, struct = data
	endif else begin
		BOX_MESSAGE, ['Cannot find AIA response file', filename]
		RETURN, filename
	endelse
	if not noblend then data = AIA_BP_BLEND_CHANNELS(data)
	if not all then chanlist = ['A94', 'A131', 'A171', 'A193', 'A211', 'A304', 'A335']
	full_resp = AIA_BP_PARSE_EFFAREA(data, short_resp, dn = dn, channels = chanlist, $
		version = version, sampleutc = sampleutc, resptable = tabledat, $
		evenorm = evenorm, timedepend = timedepend, ver_date = ver_date)
	if full then result = full_resp else result = short_resp
	RETURN, result
endif


end	


