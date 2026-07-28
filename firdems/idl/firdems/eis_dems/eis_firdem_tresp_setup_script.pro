.compile eis_firdem_add_line_to_dem.pro
den = 10^(9.0)
trfile=string(format='(%"eis_tresp_d%g.sav")',long(round(alog10(den)*10)))
print, "Creating EIS temperature response file ", trfile
eis_firdem_add_line_to_dem,'fe_09',188.497,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		/make_new, density=den
eis_firdem_add_line_to_dem,'fe_10',184.537,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_12',195.119,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_15',284.163,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_16',262.976,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'mg_05',276.579,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'mg_06',270.394,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'mg_07',280.737,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_09',197.862,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_11',188.216,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_11',180.401,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'s_10',264.233,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_12',192.394,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_13',202.044,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_13',203.826,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_14',270.519,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'s_13',256.686,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'ca_14',193.874,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'ca_15',200.972,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'ca_16',208.604,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'ca_17',192.858,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'Al-thick',0.0,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		/xrt_emiss, density=den

eis_firdem_add_line_to_dem,'si_07',275.361,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'si_10',258.371,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
eis_firdem_add_line_to_dem,'fe_14',264.790,dat=dat,fit_struc=fit_struc,save_file=trfile, $
		density=den
