;+
; Name: eis_firdem_add_line_to_dem
;
; Updates the EIS 'temperature response' structure with a new spectral line
; using emissivities computed by the CHIANTI subroutine gofnt. Saved
; to 'eis_tresp.sav'
;
; Inputs:
;
;	ion_name: 		The name of the ion to be added (e.g., 'fe_12').
;
;	line_center:	The line center in Angstroms.
;	
; Optional Inputs:
;
;	wdfilehead:		An EIS window data filename, minus the '.fits' suffix,
;					to fit the line at the same time at the same time
;					as it is added to the emissivity set.
;
;	dat:			The window data.
;
;	fit_struc: 		The line fit structure returned by eis_auto_fit. also saved
;					to wdfilehead+'_'+line_id_ref+'.fit_struc'.
;
;	make_new: 		Start over from scratch with 'eis_tresp.sav'. The file
;					will then have only the current line.
;
;	abund_file:		The name of the abundance file. Default is
;					'/ssw/packages/chianti/dbase/abundance/sun_coronal_ext.abund'
;
; 	ioneq_file:		The ionization equilibrium file. Default is
;					'/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'
;
;	density:		The density to assume. Default is 10^9.5.
;
;	xrt_emiss:		Add an effective emissivity for an XRT filter instead of for an
;					EIS line. In this case, ion_name should be the name of the filter as
;					returned by make_xrt_wave_resp and make_xrt_temp_resp.
;					line_center is added to the window name and output structure variable as
;					usual, but is otherwise ignored. Effective emissivity is computed using
;					make_xrt_wave_resp and make_xrt_temp_resp, with DN to erg conversion
;					based on Solar Physics 269, 169-236. So far, this functionality has only
;					been tested with the Al-thick filter.
;
;	contam_thick:	Contamination thickness vector of the type used by make_xrt_wave_resp. 
;					Ignored unless xrt_emiss is set. Zero thickness is assumed by default.
;
;	save_file:		Emissivities will be saved to this IDL save file. They are in the form of
;					a structure called tr_struct, which has a similar format to the result of
;					aia_get_response(/temp):
;						- eis_windows: String identifying each emissivity in the file.
;						- logte: array defining the temperature axis for all emissivities
;						- all: array (n_temps,n_windows) containing the emissivity values
;						- line_centers: array containing wavelength centers for each line
;						- ion_names: array listing ion names (or XRT filters).
;					Default save filename is 'eis_tresp.sav'
;
; Joseph Plowman (plowman@physics.montana.edu) 
; 07-11-12
;	
;-

;+
; Name: get_xrt_emiss
;
; Computes an effective emissivity for an XRT filter. uses
; make_xrt_wave_resp and make_xrt_temp_resp, with DN to erg conversion
; based on Solar Physics 269, 169-236. So far, this functionality has only
; been tested with the Al-thick filter.
;
; Inputs:
;	
;	filter: String identifying desired filter. should match filter names used by
;	make_xrt_wave_resp.
;
; Optional Inputs:
;	
;	contam_thick: contamination thickness array of the type used by make_xrt_wave_resp
;
; Outputs:
;	
;	t: Temperature array for the output emissivity
;	g: output emissivity array
;
; Joseph Plowman (plowman@physics.montana.edu)
; 07-11-12
;
;-
pro get_xrt_emiss,t,g,filter,contam_thick=contam_thick

	if(n_elements(contam_thick) eq 0) then contam_thick=[0.0,0.0]
	xrtwr=make_xrt_wave_resp(contam_thick=contam_thick)
	xrttr=make_xrt_temp_resp(make_xrt_wave_resp(contam_thick=[0.0,0.0]),/apec_default)
	
	filter_index=value_locate(xrtwr.name,filter)
	
	ccd_energy=3.65
	ev_conv=1.602e-12
	gain=57.5
	l_focal=270.8
	sr_conv=(13.5e-4)^2
	conversion_fac=ccd_energy*ev_conv*gain*l_focal^2/sr_conv
	
	t=reform(xrttr(filter_index).temp)
	it_valid = where(t gt 0)
	t=t(it_valid)
	g=conversion_fac*reform(xrttr(filter_index).temp_resp)
	g=g(it_valid)
	
end

pro eis_firdem_add_line_to_dem, ion_name, line_center, wdfilehead=wdfilehead, dat=dat, $
		fit_struc=fit_struc, make_new=make_new, abund_file=abund_file, ioneq_file=ioneq_file, $
		density=density,xrt_emiss=xrt_emiss, contam_thick=contam_thick, save_file=save_file

	line_id_ref = string(format='(%"%s_%g")',ion_name,line_center)
	if(n_elements(save_file) eq 0) then save_file = 'eis_tresp.sav'
	
	wmin = line_center - 0.01
	wmax = line_center + 0.01

	if(n_elements(abund_file) eq 0) then begin
		abund_file = '/ssw/packages/chianti/dbase/abundance/sun_coronal_ext.abund'
		print,'Abundance: using ',abund_file
	endif
	if(n_elements(ioneq_file) eq 0) then begin
;		ioneq_file = '/ssw/packages/chianti/dbase/ioneq/bryans_etal_09.ioneq'
		ioneq_file = '/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'
		print,'Ion equilibrium: using ',ioneq_file
	endif
	if(n_elements(density) eq 0) then begin
		density = 10^(9.0)
		print,'Density: using ',density
	endif

	nt = 100
	
	tmin = 10^5.0
	tmax = 10^8.0
	print,tmin,tmax
	logt = alog10(tmin) + dindgen(nt)*(alog10(tmax)-alog10(tmin))/(nt-1)
	t = 10^logt
		
	tr_array = dblarr(nt)

	print,'Please select the following spectral line: ',line_id_ref
	if(n_elements(xrt_emiss) gt 0) then begin
		get_xrt_emiss,t2,g,ion_name,contam_thick=contam_thick
	endif else begin
		gofnt,ion_name,wmin,wmax,t2,g,abund_name=abund_file,ioneq_name=ioneq_file, density=density
	endelse
	g2 = spl_init(t2,g,/double)
	tmin2 = min(t2)
	tmax2 = max(t2)
	itmin2 = min(where(t ge tmin2))
	itmax2 = max(where(t le tmax2))
	if itmin2 ge 0 and itmax2 ge 0 then begin
		tr_array(itmin2:itmax2) = spl_interp(t2,g,g2,t(itmin2:itmax2))
	endif

	if(keyword_set(make_new)) then begin
		nl = 1
		line_id_refs = replicate(line_id_ref,nl)
		line_centers = replicate(line_center,nl)
		ion_names = replicate(ion_name,nl)
		tr_arrays = dblarr(nt,nl)
		tr_arrays(*,nl-1) = tr_array
		tr_struct = {eis_windows:line_id_refs, logte:logt, all:tr_arrays, $
				line_centers:line_centers, ion_names:ion_names}
	endif else begin
		restore,save_file
	
		nl = n_elements(tr_struct.eis_windows)
		line_found = 0
		for i=0,nl-1 do begin
			if tr_struct.eis_windows(i) eq line_id_ref then begin
				tr_struct.all(*,i) = tr_array
				tr_struct.line_centers(i) = line_center
				tr_struct.ion_names(i) = ion_name
				line_found = 1
				break
			endif
		endfor
		
		if line_found eq 0 then begin
			nl = nl+1
			print,'nl =',nl
			line_id_refs = replicate(line_id_ref,nl)
			line_centers = replicate(line_center,nl)
			ion_names = replicate(ion_name,nl)
			tr_arrays = dblarr(nt,nl)
			tr_arrays(*,nl-1) = tr_array
			for i=0,nl-2 do begin
				line_id_refs(i) = tr_struct.eis_windows(i)
				line_centers(i) = tr_struct.line_centers(i)
				ion_names(i) = tr_struct.ion_names(i)
				tr_arrays(*,i) = tr_struct.all(*,i)
			endfor
			tr_struct = {eis_windows:line_id_refs, logte:logt, all:tr_arrays, $
					line_centers:line_centers, ion_names:ion_names}
		endif
	endelse
	
	save, tr_struct, filename = save_file
	
	
	if(n_elements(wdfilehead) ne 0) then begin
		wdfile = findfile(wdfilehead+'.fits')
		dat=eis_getwindata(wdfile,line_center)
		print,'Please select the relevant spectral region for ',line_id_ref
		eis_wvl_select, dat, wvl_select
		eis_auto_fit, dat, fit, wvl_select=wvl_select
		
		fit_int = reform(fit.int,fit.nx,fit.ny)
		ibadpx = where(fit.bad_pix gt 0, count)
		fit_int(ibadpx) = 0.0
		cen = reform(fit.aa(1,*,*),fit.nx,fit.ny) + fit.refwvl(0)
		fit_struc = create_struct('fit_int',fit_int,'cen',cen,'line_id',line_id_ref,fit)
		save,fit_struc,filename = wdfilehead+'_'+line_id_ref+'.fit_struc'
	endif

end
