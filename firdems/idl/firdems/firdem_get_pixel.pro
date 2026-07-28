;+
; Name: firdem_get_pixel
;
; Returns the DEM for a single pixel, given a DEM structure computed by firdem
; (or its EIS or AIA wrappers).
;
; Inputs:
;
;	dem_struct: The DEM structure created by firdem
;
;	ix,iy: The desired pixel indices.
;
; Optional Input:
;
;	zeronegs: Zero negative emission in the DEM
;
; Outputs:
;
;	t: The temperature grid for the output DEM
;
;	dem: The output DEM
;
; Joseph Plowman (plowman@physics.montana.edu) 07-31-12
;
;-
pro firdem_get_pixel, dem_struct, ix, iy, t, dem, zeronegs=zeronegs

	t=dem_struct.t
	basis = dem_struct.basis
	coffs = dem_struct.coffs(ix,iy,*)
	nb = n_elements(basis(0,*))
	
	dem = dblarr(n_elements(t))
		
	for i=0,nb-1 do begin
		dem += coffs(i)*basis(*,i)
	endfor

	if(n_elements(zeronegs) gt 0 ) then dem = dem > 0

end
