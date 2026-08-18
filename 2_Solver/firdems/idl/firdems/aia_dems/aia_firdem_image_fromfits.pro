;+
; Name: aia_firdem_demimage_fromfits
;
; Computes DEMs for a set of AIA fits files and uses it to make a jpeg image showing
; the emission weighted median temperature (EMWMT) as hue and square root total emission measure
; as intensity.
;
; Inputs:
;	filenames: the names of the AIA fits files to use
;	tmin:		Minimum of (base 10 log) temperature range for the image.
;	tmax:		Maximum of (base 10 log) temperature range for the image.
;	channels:	The AIA channel names corresponding to each file. Must match names used
;				by aia_get_response.
;	outfile:	The name of the jpeg output image.
;
; Optional inputs:
;
;	outsize:	Pixel dimensions of output. Handed to read_sdo as 'dim'.
;	tr:			Temperature response function returned by aia_get_response(/temp,/dn).
;	xscl:		Image scaling factor in x direction.
;	yscl:		Image scaling factor in y direction.
;	pmax:		Maximum total emission measure for image intensity. If unassigned, it is
;				set based on the quantiles of the image. Default:
;				pm_fac*estimate_quantile(pm_qnt)
;	pm_fac:		Quantile multiple to use for pmax. Does nothing if pmax is already assigned. 
;				Default:1.25
;	pm_qnt:	Quantile to use for pmax. Does nothing if pmax is already assigned. 
;				Default:0.995
;	corner:		Image corner for cropping used by read_sdo.
;	subfield:	Image subfield for cropping used by read_sdo.
;	noborder: 	Don't draw title and axis labels
;	tiff:		Plot to tiff instead of jpg
;	sbtiff:		Plot to 16-bit tiff instead of jpg. Forces /noborder, supercedes /tiff.
;	virt_tresp: Compute the colors with a set of 'virtual' temperature response functions,
;				instead of the EMWMT and total emission measure. Defaults to a set of log-normal
;				functions; can also use user-supplied function.
;	vtr_cens:	The centers of the virtual temperature response functions mentioned above. 
;				Implies virt_tresp. Default: [5.85,6.3,6.75]
;	vtr_sigs:	The sigmas of the virtual log-normal temperature response functions mentioned
;				above. Note that these functions peak at unity, regardless of width. Implies
;				virt_tresp. Default: [0.24,0.08,0.24]
;	vtr_logt:	Log10(T) axis for virtual response functions. Implies virt_tresp. Default:
;				200 points uniform in logt ranging from 5.5 to 7.5.
;	vtrs:		Array (dimensions [3,n_elements(vtr_logt)]) containing virtual temperature
;				response functions. Must correspond to temperatures in vtr_logt, which must
;				also be supplied. Implies virt_tresp. Default: log-normal distributions with
;				centers and sigmas specified by vtr_cens and vtr_sigs.
;
; Optional Output:
;
;	dems: The DEM structure computed by aia_dem_iterative_extrap.
;
; Joseph Plowman (plowman@physics.montana.edu) 09-19-12
;
;-
pro aia_firdem_demimage_fromfits, filenames, channels, outfile, tmin=tmin, tmax=tmax, dim=dim, $ 
		tr=tr, xscl=xscl, yscl=yscl, pmax=pmax, pm_fac=pm_fac, pm_qnt=pm_qnt, corner=corner, $
		subfield=subfield, tiff=tiff, sbtiff=sbtiff, virt_tresp=virt_tresp, vtr_cens=vtr_cens, $
		vtr_sigs=vtr_sigs, vtr_logt=vtr_logt, vtrs=vtrs, dems=dems

	tstart=systime(1)
	nc = n_elements(channels)
	print,"Calculating DEM from",filenames
	
	read_sdo,filenames,indices0,/nodata
	index = indices0(0)
	nx = min(indices0.naxis1)-1
	ny = min(indices0.naxis2)-1
	if((n_elements(corner) gt 0) and (n_elements(subfield) gt 0)) then begin
		nx = subfield(0)
		ny = subfield(1)
	endif
	if(n_elements(xscl) gt 0) then nx = round(xscl*nx)
	if(n_elements(yscl) gt 0) then ny = round(yscl*ny)	
	dim=[nx,ny]
	
	data = dblarr(nx,ny,nc)
	for i=0,nc-1 do begin
		if((n_elements(corner) gt 0) and (n_elements(subfield) gt 0)) then begin
			read_sdo,filenames(i), indextemp, datatemp, corner(0), corner(1), subfield(0), %
					subfield(1), outsize=dim,/uncomp_delete
		endif else begin
			read_sdo,filenames(i),indextemp,datatemp,outsize=dim,/uncomp_delete
		endelse
		aia_prep,indextemp,datatemp,indextemp_prep,datatemp_prep,/use_ref,index_ref=refindex
		if(i eq 0) then begin
			indices = replicate(indextemp_prep,nc)
			indextemplate = indextemp_prep
			refindex = indextemplate
		endif else begin
			struct_assign,indextemp_prep,indextemplate
			indices(i) = indextemplate
		endelse
		data(*,*,i) = datatemp_prep
	endfor
	
	data /= (xscl*yscl)
	
	print,"nx=",nx," ny=",ny
		
	aia_firdem_wrapper, channels, indices.exptime, data, dems, tr=tr

	if(keyword_set(virt_tresp)) then begin
		firdem_vchan_image, dems, image, vtr_cens=vtr_cens, vtr_sigs=vtr_sigs, $
				vtr_logt=vtr_logt, vtrs=vtrs
		if(n_elements(vtr_cens) eq 0) then begin
			titlestr = indices(0).t_obs	
		endif else begin
			tminstr=string(format='(%"%g")',vtr_cens(0))
			tmidstr=string(format='(%"%g")',vtr_cens(1))
			tmaxstr=string(format='(%"%g")',vtr_cens(2))
			titlestr = indices(0).t_obs + ' log(tmin)=' + tminstr + ' log(tmid)=' + tmidstr + $
					' log(tmax)='+tmaxstr
		endelse
	endif else begin
		tminstr=string(format='(%"%g")',tmin)
		tmaxstr=string(format='(%"%g")',tmax)
		titlestr = indices(0).t_obs +' (EMWMT)' + ' log(tmin)=' + tminstr + ' log(tmax)=' + $
				tmaxstr
		firdem_emwtemp_image_color,dems,tmin,tmax,image
	endelse

	firdem_imagefile_plot, image, outfile, pmax=pmax, pm_fac=pm_fac, pm_qnt=pm_qnt, $
			title=titlestr, xlabel=xlabel, ylabel=ylabel, origin=origin, scale=scale, $
			noborder=noborder, tiff=tiff, sbtiff=sbtiff	

	print,"Wrote ",outfile
	print,"aia_firdem_demimage_fromfits finished. ",systime(1)-tstart," seconds elapsed."

end
