;+
; Name: eis_dem_warren_comparison
;
; This procedure computes the DEM using EIS data from Warren et al 2011 (ApJ 734, 90) and plots
; it. The data are taken from their Table 1.
;
; Outputs:
;	temps_out: Temperature grid for the output DEM
;	dem_out: The output DEM
;-
pro eis_firdem_warren_comparison, temps_out, dem_out

	common aia_dem_iterative_test_1_static_variables, seed
	restore,'eis_tresp_d95.sav'

	channels_all = tr_struct.eis_windows
	channels = ["fe_09_188.497", "fe_10_184.537", "fe_12_195.119", "fe_15_284.163", "fe_16_262.976"];, "ca_14_193.874","ca_15_200.972","s_13_256.686"]
;	channels_full = ["fe_09_188.497", "fe_10_184.537", "fe_12_195.119", "fe_15_284.163", "fe_16_262.976","mg_05_276.579","mg_06_270.394","mg_07_280.737","fe_09_197.862","fe_11_188.216","fe_11_180.401","s_10_264.233","fe_12_192.394","fe_13_202.044","fe_13_203.826","fe_14_270.519","s_13_256.686","ca_14_193.874","ca_15_200.972","ca_16_208.604","ca_17_192.858"]
	channels_full = ["fe_09_188.497", "fe_10_184.537", "fe_12_195.119", "fe_15_284.163", "fe_16_262.976","mg_05_276.579","mg_06_270.394","mg_07_280.737","fe_09_197.862","fe_11_188.216","fe_11_180.401","s_10_264.233","fe_12_192.394","fe_13_202.044","fe_13_203.826","fe_14_270.519","s_13_256.686","ca_14_193.874","ca_15_200.972","ca_16_208.604","ca_17_192.858","si_07_275.361","si_10_258.371","fe_14_264.79","Al-thick_0"]

	images0 = [72.4,	280.4,	1475.4,	10334.0,	1157.6];,	311.9,	238.9,	854.7]
	errors0 = [16.0,	61.8,	324.6,	2273.6,		254.7];,	68.6,	52.6,	188.1]
	
;	images0_full = [72.4,	280.4,	1475.4,	10334.0,	1157.6,	16.5,	35.9,	32.7,	40.0,	578.2,	926.1,	71.8,	437.8,	1248.3,	2533.9,	515.0,	854.7,	311.9,	238.9,	122.0,	146.5]
;	errors0_full = [16.0,	61.8,	324.6,	2273.6,		254.7,	3.8,	7.9,	7.4,	8.8,	127.2,	204.2,	73.9,	478.5,	665.9,	1331.4,	523.0,	853.2,	276.3,	195.7,	108.7,	128.5]
	
	images0_full = [72.4,	280.4,	1475.4,	10334.0,	1157.6,	16.5,	35.9,	32.7,	40.0,	578.2,	926.1,	71.8,	437.8,	1248.3,	2533.9,	515.0,	854.7,	311.9,	238.9,	122.0,	146.5,	47.0,	294.0,	1026.9, 4.5]
	errors0_full = [16.0,	61.8,	324.6,	2273.6,		254.7,	3.8,	7.9,	7.4,	8.8,	127.2,	204.2,	73.9,	478.5,	665.9,	1331.4,	523.0,	853.2,	276.3,	195.7,	108.7,	128.5,	10.4,	64.7,	226.0, 0.9]
	
;	nonmg_chans = where(channels_full ne "mg_05_276.579") and  where(channels_full ne "mg_06_270.394") and where(channels_full ne "mg_07_280.737")
;	channels_full = channels_full(nonmg_chans)
;	images0_full = images0_full(nonmg_chans)
;	errors0_full = errors0_full(nonmg_chans)

	nchan_all = n_elements(channels_all)
	nchan = n_elements(channels)
	nchan_full = n_elements(channels_full)
	
	tr_struct.all(*,where(channels_all eq "mg_05_276.579")) *= 1.7
	tr_struct.all(*,where(channels_all eq "mg_06_270.394")) *= 1.7
	tr_struct.all(*,where(channels_all eq "mg_07_280.737")) *= 1.7
	
	tr_struct.all(*,where(channels_all eq "Al-thick_0")) *= 1.5

	nx=25
	ny=40

	images = dblarr(nx,ny,nchan)
	errors = dblarr(nx,ny,nchan)
	images_full = dblarr(nx,ny,nchan_full)
	errors_full = dblarr(nx,ny,nchan_full)
	
	
	for i=0,nchan-1 do begin
		images(*,*,i) = images0(i)
		errors(*,*,i) = errors0(i)*1.0
	endfor

	for i=0,nchan_full-1 do begin
		images_full(*,*,i) = images0_full(i)
		errors_full(*,*,i) = errors0_full(i)*1.0
	endfor

	chi2_target=chisqr_cvf(0.5,nchan)

	eis_firdem_wrapper, channels, images, errors, dem_out, tr_struct=tr_struct

	chi2_target=chisqr_cvf(0.5,nchan_full)
	tstart=systime(1)
	eis_firdem_wrapper, channels_full, images_full, errors_full, dem_out_full, $
			tr_struct=tr_struct, chi2s=chi2s, chi2_target=chi2_target, $
			logtlo=5.5, logthi=7.5

	telapsed = systime(1)-tstart

	chi2_95 = estimate_quantile(chi2s,0.95)
	title = string(format='(%"%8.5g%s%5.3g")',1000.0*telapsed/nx/ny," ms per DEM; 95% Chi Squared=",chi2_95)
	
	pi = 2*acos(0)

	firdem_get_pixel,dem_out,0,0,t2,dem2
	firdem_get_pixel,dem_out_full,0,0,t2_full,dem2_full
	plot, alog10(t2_full), dem2_full*0.1, linestyle=0, yrange=[1e25,1e30], xrange=[5.25,7.35], $
			xtitle = "Log!D10!N(T) (Kelvin)", ytitle="DEM (cm!E-5!N [0.1 Log!D10!N(T)]!E-1!N)", $
			title=title, xstyle=1, /ylog, thick=3
	for i=0,nchan_full-1 do begin
		tr = tr_struct.all(*,where(tr_struct.eis_windows eq channels_full(i)))
		emloc = 1.0*images0_full(i)/tr
		plotted_trpoints = where(tr gt 5.0e-3*max(tr))
;		emloc(where(tr lt 5.0e-3*max(tr))) = 1.0e30
		oplot,tr_struct.logte(plotted_trpoints),emloc(plotted_trpoints),linestyle=2,thick=3
	endfor
	legend,['Extracted DEM','EM Locii'],linestyle=[0,2],/right,/bottom,thick=3
	print, "Total EMs:",int_tabulated(alog10(t2),dem2),int_tabulated(alog10(t2_full),dem2_full)
	temps_out = t2
	dem_out = dem2
	
	
end

pro plot_basis,dem_out
	plot, dem_out.t, dem_out.basis(*,0), /xlog, yrange=[0,max(dem_out.basis)], xtitle='Temperature(Kelvin)',title='Normalized EIS Emissivities/Basis Functions'
	for i=0,n_elements(dem_out.basis(0,*))-1 do begin
		oplot,dem_out.t,dem_out.basis(*,i)
	endfor
end
