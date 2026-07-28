function AIA_BP_READ_ERROR_TABLE, tablefile, $
	silent = silent

;+
;
; $Log: aia_bp_read_error_table.pro,v $
; Revision 1.1  2012/02/10 01:39:22  boerner
; Initial version
;
; Revision 1.1  2011/10/03 23:21:06  boerner
; Initial version
;
; 
; INPUTS:
;	tablefile	--	string giving the name of the _error_table.txt file to be read in
;
; RETURNS:
;	tablestr	--	a structure array containing the data from the error table file
;
;-

silent = KEYWORD_SET(silent)
if N_ELEMENTS(tablefile) eq 0 then begin
	tablefile = GET_LOGENV('AIA_RESPONSE_DATA') + '/aia_V2_error_table.txt'
endif

hastable = FILE_EXIST(tablefile)
if hastable then begin
	tabledat = RD_TFILE(tablefile, 11)
endif else begin
	if not silent then BOX_MESSAGE, 'AIA_BP_READ_ERROR_TABLE: Missing file ' + tablefile
	RETURN, 0
endelse

table_template = {Date : '', t_start : '', t_stop : '', ver_num : 0, $
	wave_str : '', wavelnth : 0, dnperpht : 0d, compress : 0d, calerr : 0d, $
	chianti : 0d, eveerr : 0d}

tablesize = SIZE(tabledat)
for i = 1, tablesize[2] - 1 do begin
	thisline = table_template
	thisline.date 		= tabledat[0,i]
	thisline.t_start	= tabledat[1,i]
	thisline.t_stop		= tabledat[2,i]
	thisline.ver_num	= tabledat[3,i]
	thisline.wave_str	= tabledat[4,i]
	thisline.wavelnth	= tabledat[5,i]
	thisline.dnperpht	= tabledat[6,i]		;	number of DN per photon
	thisline.compress	= tabledat[7,i]		;	Ratio of shot noise to error from onboard (LUT) compression
	thisline.calerr		= tabledat[8,i]		;	Systematic error in photometric calibration (%)
	thisline.chianti	= tabledat[9,i]		;	Systematic error in CHIANTI portion of Tresponse (%)
	thisline.eveerr		= tabledat[10,i]	;	Error in effarea if /evenorm is set (rmse from resp table) (%)
	if i eq 1 then tablestr = thisline else tablestr = [tablestr, thisline]
endfor

RETURN, tablestr

end