;+
; Name: firdem_vchan_image
;
; Makes a true color image from a DEM computed by firdem using 3 'virtual' channels with
; specified temperature response functions. The response functions may either be specified
; directly using vtr_logt and vtrs, or they may be specified as log-normal distributions 
; by giving their centers and sigmas as vtr_cens, and vtr_sigs.
;
; Inputs:
;
;	dem_struc: The DEM structure returned by firdem or its EIS and AIA wrappers.
;
; Outputs:
;
;	image: True color image made from the DEM using the 3 supplied virtual response functions.
;
; Optional Outputs:
;
;	vtr_cens: 3-element array containing the centers for the 3 log-normal virtual
;				temperature response functions used if vtrs is not set. Default: [5.85,6.3,6.75]
;	vtr_sigs: 3-element array containing the sigmas for the 3 log-normal virtual
;				temperature response functions used if vtrs is not set. Default:
;				[0.24,0.08,0.24]
;	vtr_amps: 3-element array containing the peak amplitudes for the 3 log-normal virtual
;				temperature response functions used if vtrs is not set. Default:
;				[1.0,1.0,1.0]
;
;	vtr_logt: Array containing the logt values at which the response functions will be 
;				evaluated. Default: 200 points uniform from 5.5 to 7.5, inclusive.
;	vtrs: Array (Dimensions [3,ntemps]) containing the values of the 3 virtual temperature
;			response functions at the temperatures specified by vtr_logt. Does nothing
;			if vtr_logt is unassigned.
;
; Joseph Plowman (plowman@physics.montana.edu) 10-11-12
;-
pro firdem_vchan_image, dem_struc, image, vtr_cens=vtr_cens, vtr_sigs=vtr_sigs, $
		vtr_amps=vtr_amps, vtr_logt=vtr_logt, vtrs=vtrs

	t=dem_struc.t
	basis0 = dem_struc.basis
	nb = n_elements(basis0(0,*))

	if(n_elements(vtr_logt) eq 0) then begin
		vtr_logtmin=5.5
		vtr_logtmax=7.5
		nt=200
		vtr_logt = vtr_logtmin+(vtr_logtmax-vtr_logtmin)*findgen(nt)/(nt-1.0)
		vtrs=0
	endif
	nt = n_elements(vtr_logt)
	if(keyword_set(vtrs) eq 0) then begin
		vtrs = dblarr(3,nt)
		pi=2.0*acos(0)
		if(n_elements(vtr_cens) eq 0) then vtr_cens = [5.85,6.3,6.75]
		if(n_elements(vtr_sigs) eq 0) then vtr_sigs = [0.24,0.08,0.24]
		if(n_elements(vtr_amps) eq 0) then vtr_amps = [1.0,1.0,1.0]
		for i=0,2 do vtrs(i,*) = $
				vtr_amps(i)*exp(-(vtr_cens(i)-vtr_logt)^2/(2.0*vtr_sigs(i)^2))/sqrt(2*pi)
	endif

	nt_dem = n_elements(t)
	logt_dem = alog10(t)
		
	basis=dblarr(nt,nb)
	for i=0, nb-1 do begin
		basis02 = spl_init(logt_dem,basis0(*,i),/double)
		basis(*,i) = spl_interp(logt_dem,basis0(*,i),basis02,vtr_logt)
	endfor
	
	nresps = n_elements(vtrs(*,0))
	
	maparr = fltarr(nresps,nb)
	for i=0,nresps-1 do begin
		for j=0,nb-1 do begin
			maparr(i,j)=int_tabulated(vtr_logt,vtrs(i,*)*basis(*,j))
		endfor
	endfor

	nx = n_elements(dem_struc.coffs(*,0,0))
	ny = n_elements(dem_struc.coffs(0,*,0))

	image = fltarr(nx,ny,nresps)
	coffsij = fltarr(nb)

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			coffsij(*) = dem_struc.coffs(i,j,*) > 0
			image(i,j,*) = maparr#coffsij
		endfor
	endfor	

end
