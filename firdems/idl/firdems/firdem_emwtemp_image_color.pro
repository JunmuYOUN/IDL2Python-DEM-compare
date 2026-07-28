;+
; Name: firdem_emwtemp_image_color
;
; Makes a true color image from a DEM computed by fastdem_iterative (or its EIS or AIA
; wrappers), showing emission weighted median temperature (EMWMT) as hue and total emission
; measure as intensity.
;
; Inputs:
;
;	dem_struct: The DEM structure returned by aia_dem_iterative_extrap
;	logtmin:		Minimum of temperature range for color the image.
;	logtmax:		Maximum of temperature range for color the image.
;
; Outputs:
;
;	image:			True color image showing EMWMT (hue) and total emission measure (brightness)
;					for each pixel in the image.
;
; Optional Input:
;
;	keepnegs: 		Don't zero negative emission measure before computing median.
; 
; Optional Outputs:
;
;	emwmt:			Emission weighted median (base 10 log) temperature of each pixel.
;	emtot:			Total emission measure in each pixel.
;-
pro firdem_emwtemp_image_color, dem_struct, logtmin, logtmax,image, emwmt=emwmt, emtot=emtot,keepnegs=keepnegs

	t=dem_struct.t
	basis = dem_struct.basis
	coffs = dem_struct.coffs
	nb = n_elements(basis(0,*))

	nt = n_elements(t)
	logt = alog10(t)
	
	logtrange = logtmax-logtmin

	logtmid = logtmin+0.5*logtrange

	nc = 3.0	

	nx = n_elements(coffs(*,0,0))
	ny = n_elements(coffs(0,*,0))

	image = dblarr(nx,ny,nc)
	emwmt = dblarr(nx,ny)
	emtot = dblarr(nx,ny)

	coffsij = dblarr(nb)

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			coffsij(*) = coffs(i,j,*)
			dem = basis#coffsij
			if(keyword_set(keepnegs) eq 0) then dem = dem > 0
			logtavg = dist_median(logt,dem,tot=emtot_current)
			rgb = firdem_rgbmap_triangle(logtavg,logtmin,logtmax)
			for k=0,nc-1 do image(i,j,*) = sqrt(emtot_current)*rgb
			emwmt(i,j)=logtavg
			emtot(i,j)=emtot_current
		endfor
	endfor	
	
end
