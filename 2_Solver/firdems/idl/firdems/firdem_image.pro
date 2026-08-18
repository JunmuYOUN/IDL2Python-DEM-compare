;+
; Name: firdem_image
;
; Creates an image at a specified temperature, from a dem structure returned by
; fastdem_iterative or its AIA and EIS wrappers.
;
; Inputs:
; 	temp: the desired temperature
; 	dem_struct: the DEM structure
;
; Returns:
;	image containing the DEM value at the specified temperature, for each pixel.
;
; Joseph Plowman (plowman@physics.montana.edu) 07-31-12
;
;-
function firdem_image,temp,dem_struct

	t=dem_struct.t
	basis = dem_struct.basis
	coffs = dem_struct.coffs
	nb = n_elements(basis(0,*))
	
	data = dblarr(n_elements(coffs(*,0,0)),n_elements(coffs(0,*,0)))
	
	iarr = where(abs(t-temp) eq min(abs(t-temp)))
		
	if(iarr eq n_elements(t)-1) then iarr--
	
	t1=t(iarr)
	t2=t(iarr+1)
	dt2=(temp-t1)/(t2-t1)
	dt1=(t2-temp)/(t2-t1)
	
	for i=0,nb-1 do begin
		a1 = basis(iarr,i)
		a2 = basis(iarr+1,i)
		data+=coffs(*,*,i)*(a1(0)*dt1(0)+a2(0)*dt2(0))
	endfor
	
	return,data
	
end
