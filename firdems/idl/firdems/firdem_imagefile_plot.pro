;+
; Name: firdem_imagefile_plot
;
; Plots a truecolor image (such as one produced by fastdem_vchan_image or
; fastdem_emwtemp_image_color) to a jpeg (the default) or tiff file.
;
; Inputs:
;
;	image: The image to be plotted. Dimensions are assumed to be [nx,ny,nchan], with [0,0,0]
;	being the bottom right pixel of the red channel.
;	filename: The filename for the plot
;
; Optional arguments:
;
;	pmax:		Maximum total emission measure for image intensity. If unassigned, it is
;				set based on the quantiles of the image. Default:
;				pm_fac*estimate_quantile(pm_qnt)
;	pm_fac:		Quantile multiple to use for pmax. Does nothing if pmax is already assigned. 
;				Default:1.25
;	pm_qnt:		Quantile to use for pmax. Does nothing if pmax is already assigned. 
;				Default:0.995
;	title: 		Title to use for the plot.
;	xlabel:		x axis label.
;	ylabel:		y axis label.
;	origin:		Plot origin.
;	scale:		Plot scale.
;	noborder: 	Don't draw title and axis labels.
;	tiff:		Plot to tiff instead of jpg.
;	sbtiff:		Plot to 16-bit tiff instead of jpg. Forces /noborder, supercedes /tiff.
;
; Joseph Plowman (plowman@physics.montana.edu) 8-24-12
;-
pro firdem_imagefile_plot, image, filename, pmax=pmax, pm_fac=pm_fac, pm_qnt=pm_qnt, title=title, xlabel=xlabel, ylabel=ylabel, origin=origin, scale=scale, noborder=noborder, tiff=tiff, sbtiff=sbtiff

	nx=n_elements(image(*,0,0))
	ny=n_elements(image(0,*,0))
	
	if(n_elements(pmax) eq 0) then begin
		if(n_elements(pm_fac) eq 0) then pm_fac=1.25
		if(n_elements(pm_qnt) eq 0) then pm_qnt=0.995
		pmax = pm_fac*estimate_quantile(image,pm_qnt)
	endif

	if(keyword_set(sbtiff)) then noborder=1
	
	if(keyword_set(noborder) eq 0) then begin
		previous_device=!D.name
		set_plot,'z'
		charysize = max([15,15.0*ny/1024.0])
		charsize = round([0.6*charysize,charysize])
		device, set_pixel_depth=24, set_resolution=round([1.3*nx,1.3*ny]), set_character_size=charsize, z_buffering=0
		device,decomposed=1
		plot_image, image, origin=origin, scale=scale, min=0, max=pmax, title=title, xlabel=xlabel, ylabel=ylabel
		device,decomposed=0
		image_out=tvrd(true=3)
		set_plot,previous_device
	endif else begin
		image_out = image < pmax
		if(keyword_set(sbtiff)) then begin
			image_out=image_out*(2.0^16.0-1.0)/pmax
		endif else begin
			image_out=image_out*(255.0)/pmax
		endelse
	endelse

	if(keyword_set(tiff) or keyword_set(sbtiff)) then begin
		if(keyword_set(sbtiff)) then begin
			write_tiff,filename, image_out, planarconfig=2, compression=1, /short, orientation=0
		endif else begin
			write_tiff,filename, image_out, planarconfig=2, compression=1, orientation=0
		endelse
	endif else begin
		write_jpeg,filename,image_out,true=3,quality=99
	endelse
	
end

