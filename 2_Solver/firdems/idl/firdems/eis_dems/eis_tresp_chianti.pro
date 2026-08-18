pro eis_tresp_chianti,tresp

;	line_id_ref = [	'CA XVII 192.470', 'CA XIV 194.100', 'CA XV 201.050', 'CA XVI 208.500', $
;					'FE XI 180.400', 'CA XV 182.100', 'FE X 184.330', 'FE VIII 185.210', $
;					'FE XII 186.750', 'FE XI 188.400', 'FE XII 195.120', 'FE IX 197.860', $ 
;					'FE XIII 202.040', 'MG VI 269.000', 'FE XIII 203.830', 'FE XXIV 255.000', $
;					'FE XVI 262.980', 'FE XIV 270.520', 'FE XV 284.160', 'MG V 276.300', $ 
;					'MG VII 278.400', 'MG VII 280.740', 'S XIII 256.480', 'S X 264.300', $ 
;					'SI X 258.370', 'SI VII 275.400' ]
					
;	ion_names = [	'ca_17', 'ca_14', 'ca_15', 'ca_16', $
;					'fe_11', 'ca_15', 'fe_10', 'fe_08', $
;					'fe_12', 'fe_11', 'fe_12', 'fe_09', $
;					'fe_13', 'mg_06', 'fe_13', 'fe_24', $
;					'fe_16', 'fe_14', 'fe_15', 'mg_05', $
;					'mg_07', 'mg_07', 's_13', 's_10', $
;					'si_10', 'si_07' ]
					
;	eis_centers = [	192.470, 194.100, 201.050, 208.500, $
;					180.400, 182.100, 184.330, 185.210, $
;					186.750, 188.400, 195.120, 197.860, $
;					202.040, 269.000, 203.830, 255.000, $
;					262.980, 270.520, 284.160, 276.300, $
;					278.400, 280.740, 256.480, 264.300, $
;					258.370, 275.400 ]


	line_id_ref = [	'CA XVII 192.470', 'CA XV 201.050', 'FE XI 180.400', 'FE X 184.330', $ 
					'FE XII 186.750', 'FE XII 195.120', 'FE IX 197.860', 'FE XIII 202.040', $
					'FE XIII 203.830', 'FE XVI 262.980', 'FE XIV 270.520', 'FE XV 284.160', $
					'S X 264.300', 'SI X 258.370', 'CA XVII 192.850' ]

	ion_names = [	'fe_12', 'fe_13', 'fe_11', 'fe_10', 'fe_12', 'fe_12', 'fe_09', 'fe_13', $
					'fe_13', 'fe_16', 'fe_14', 'fe_15', 'fe_14', 'si_10', 'ca_17' ]

	line_centers = [192.3940, 201.1280, 180.4080, 184.5370, 186.8870, 195.1190, 197.8620, $
					202.0440, 203.8280, 262.9760, 270.522, 284.1630, 264.7900, 258.3710, $
					192.8532]

					
	wmins = line_centers - 0.001
	wmaxs = line_centers + 0.001
	
	abund_file = '/ssw/packages/chianti/dbase/abundance/sun_coronal.abund'
	ioneq_file = '/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'
	
	nlines = n_elements(line_id_ref)
	
	nt = 100
	
	tmin = 10^5.0
	tmax = 10^7.5
	print,tmin,tmax
	logt = alog10(tmin) + dindgen(nt)*(alog10(tmax)-alog10(tmin))/(nt-1)
	t = 10^logt
		
	tr_array = dblarr(nt,nlines)
	
	for i=0,nlines-1 do begin
		print,line_id_ref(i)
		gofnt,ion_names(i),wmins(i),wmaxs(i),t2,g,abund_name=abund_file,ioneq_name=ioneq_file
		g2 = spl_init(t2,g,/double)
		tmin2 = min(t2)
		tmax2 = max(t2)
		itmin2 = min(where(t ge tmin2))
		itmax2 = max(where(t le tmax2))
;		print,tmin,tmax
;		print,tmin2,tmax2,itmin2,itmax2
		if itmin2 ge 0 and itmax2 ge 0 then begin
			tr_array(itmin2:itmax2,i) = spl_interp(t2,g,g2,t(itmin2:itmax2))
		endif
		;print,logt

	endfor
	
	tresp = {eis_windows:line_id_ref, logte:logt, all:tr_array, line_centers:line_centers, $
			ion_names:ion_names}
			
end
