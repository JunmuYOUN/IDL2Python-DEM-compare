FUNCTION GCVFUNCTION, trialq,flux, basekernel, A, H, data, kernel,num=num

; NAME:
;	gcvfunction.pro

; PURPOSE:
;	Part of the package of routines associated with computation of
;	Velocity Differential Emission Measure (VDEM) from spectral lines.
;	This function computes the Generalized Cross-Validation (GCV) function
;	which is minimized by the routine GCVCOMPUTE.PRO to find the
;	optimal smoothing parameter for the VDEM inversion.  Reference
;	Newton et al. 1996 on use.

; CATEGORY:
;	data analysis

; CALLING SEQUENCE:
;	gcv=gcvfunction(a1,flux,basekernel,A,H,tempdata,kernel,num=num)

; INPUTS:
;	trialq=Initial guess for optimal smoothing parameter
;	flux=Observed flux in spectral line
;	basekernel=Kernel function employed in VDEM inversion
;	A=Kernel^T * Kernel
;	H=Regularizaton matrix employed in VDEM inversion
;	data=Kernel^T*Normalized flux (normalized by errors)
;	kernel=Basekernel normalized by data errors
;	num=Number of velocity bins for which VDEM will be inverted

; OUTPUTS:
;	gcv=the GCV function that is minimized in GCVCOMPUTE.PRO

; KEYWORDS:
;	None

; COMMON BLOCKS:
;	None

; SIDE EFFECTS:
;	None

; RESTRICTIONS:
;	The inputs must have been computed in VDEM.PRO.

; ERRORS:
;	An error message is printed if the matrix inversion (performed by
;	Gaussian elimination) are suspect.

; PROCEDURE:
;	Run VDEM.PRO and utilize its default for computing the
;	optimal smoothing parameter for inversion with the GCV technique.

; MODIFICATION HISTORY:
;	E. Newton 1995: IDL implementation
;=========================================================================


;------Preparation----------------------------------------
au=1.495979e13	;Earth-Sun distance (cm)
sun_sat_dis=0.99*au	;Sun-satellite distance (cm)

nflux=n_elements(flux)

matrix0 = A + trialq*H
identity=fltarr(nflux,nflux) & for i=0,nflux-1 do identity(i,i)= 1.

;-------Compute GCV Numerator-----------------------------

if (fix(strmid(!version.release,0,1)) lt 5) then begin

   lu_matrix=transpose(matrix0)            ;nr_ludcmp requires column format
   lu_data=transpose(data)                ;nr_lubksb requires column format
   nr_ludcmp,lu_matrix,permut,/double     ;matrix is replaced by its LU
					  ;decomposition
   recovery=nr_lubksb(lu_matrix,permut,lu_data)

endif else begin

lu_matrix0=matrix0
lu_data = transpose(data)
recovery= v5_invert(lu_matrix0,lu_data)

endelse


gcvnum=total((flux/(4.*!pi*sun_sat_dis^2) - $
	(recovery # basekernel)/(4.*!pi*sun_sat_dis^2))^2)


;------Compute GCV Denominator----------------------------

inverse=invert(matrix0,status)	;finds inverse by using gaussian elimination
if (status ne 0) then $
print, 'Inverse of Matrix0 is suspect or invalid in computation of gcv'

temp1=transpose(kernel) # inverse
temp2=temp1 # kernel
temp3=identity - temp2

g_trace=temp3(0,0) & for i=0,n_elements(flux)-2 do $
	g_trace=g_trace + temp3(i+1,i+1)
gcvdenom=g_trace^2

;------Compute GCV Function-------------------------------
gcv=gcvnum/gcvdenom

return, gcv

end
