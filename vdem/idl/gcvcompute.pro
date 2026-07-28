FUNCTION GCVCOMPUTE,a1,b1,flux,basekernel,A,H,data,kernel,num=num

; NAME:
;	gcvcompute.pro

; PURPOSE:
;	Part of the package of routines associated with computation of
;	Velocity Differential Emission Measure (VDEM) from spectral lines.
;	This function computes the optimal smoothing parameter (optq)
;	for the VDEM inversion by minimizing the Generalized Cross-
;	Validation (GCV) function.  Reference Newton et al. 1996 on use.
;	Function uses Numerical Recipes' bracketing and brent
;	minimization routines for single dimension.


; CATEGORY:
;	data analysis

; CALLING SEQUENCE:
;	optq=gcvcompute(q1,q2,tempflux,basekernel,A,H,data,kernel,num=nvel)

; INPUTS:
;	q1,q2=Initial guesses to bracket optimal smoothing parameter
;	tempflux=Observed flux in spectral line (not normalized by data
;		uncertainty)
;	basekernel=Kernel function employed in VDEM inversion
;	A=Kernel^T * Kernel
;	H=Regularizaton matrix employed in VDEM inversion
;	data=Kernel^T*Normalized flux (normalized by observational errors)
;	kernel=Basekernel normalized by data errors
;	num=Number of velocity bins for which VDEM will be inverted

; OUTPUTS:
;	optq=the optimal smoothing parameter which should be employed for
;	the VDEM inversion.

; KEYWORDS:
;	None

; COMMON BLOCKS:
;	None

; SIDE EFFECTS:
;	None

; RESTRICTIONS:
;	The inputs must have been computed in VDEM.PRO.

; ERRORS:
;	If minimization routine is unsuccessful, prints error message,
;	'brent routine exceeded maximum iterations'

; PROCEDURE:
;	Run VDEM.PRO and utilize its default for computing the optimal
;	smoothing parameter for inversion with the GCV technique,
;	GCVCOMPUTE.PRO in turn calls GCVFUNCTION.PRO.

; MODIFICATION HISTORY:
;	R. Barrett 1994: Original coded in Fortran for use in Inversions
;	E. Newton 1995: IDL implementation
;=========================================================================



;------Find triplet that brackets the minimum-----------------------
glimit=100.			;max magnification for step
gold=1.618034			;golden ratio
tiny=1.e-20			;smallest fractional accuracy

tempdata=data
gcva=gcvfunction(a1,flux,basekernel,A,H,tempdata,kernel,num=num)
gcvb=gcvfunction(b1,flux,basekernel,A,H,tempdata,kernel,num=num)

;Orient so heading downhill from a1 to b1
if (gcvb gt gcva) then begin
	dummy=a1 & a1=b1 & b1=dummy & dummy=gcvb & gcvb=gcva & gcva=dummy
endif

c=b1 + gold*(b1-a1)		;first guess for third point
gcvc=gcvfunction(c,flux,basekernel,A,H,data,kernel,num=num)

;Enter loop if (gcvb ge gcvc).  Otherwise have already found the bracketing
;triplet
newtry:
if (gcvb ge gcvc) then begin
	r=(b1-a1)*(gcvb-gcvc)	;compute u by parabolic extrapolation from a,b,c
	q=(b1-c)*(gcvb-gcva)
	if (abs(q-r) ge tiny) then temp=q-r	;avoid division by zero
	if (abs(q-r) lt tiny) and (q-r ge 0.) then temp=tiny else temp=-tiny
	u=b1-((b1-c)*q -(b1-a1)*r)/(2.*temp)
	ulim=b1+glimit*(c-b1)	;limit to extrapolated u

	;test  3 possibilities
	;first possibility---u is between b&c
	if ((b1-u)*(u-c) gt 0.) then begin
		gcvu=gcvfunction(u,flux,basekernel,A,H,data,kernel,num=num)
		if (gcvu lt gcvc) then begin		;have min between b&c
			a1=b1 & gcva=gcvb & b1=u & gcvb=gcvu
			goto, bracket		;success in finding bracket
		endif else begin			;have min between a&u
			if (gcvu gt gcvb) then begin
				c=u & gcvc=gcvu
				goto, bracket	;success in finding bracket
			endif
		endelse
		goto,newu
	endif
	;second possibility---u is between c&limit
	if ((c-u)*(u-ulim) gt 0.) then begin
		gcvu=gcvfunction(u,flux,basekernel,A,H,data,kernel,num=num)
		if (gcvu lt gcvc) then begin	;still headed downhill
			b1=c & c=u & gcvb=gcvc & gcvc=gcvu
			goto,newu
		endif else goto, swap

	endif

	;third possibility---u must be limited to maximum allowable value/u is
 	;beyond ulim
	if ((u-ulim)*(ulim-c) ge 0.) then begin
	u=ulim & gcvu=gcvfunction(u,flux,basekernel,A,H,data,kernel,num=num)
		goto,swap
	endif

	newu:
	;parabolic fit is no good--use default magnification to compute new u
	u=c+gold*(c-b1)
	gcvu=gcvfunction(u,flux,basekernel,A,H,data,kernel,num=num)

	;eliminate oldest point and continue
	swap:
	a1=b1 & b1=c & c=u & gcva=gcvb & gcvb=gcvc & gcvc=gcvu
	goto, newtry

endif

bracket:		;now have the bracketing triplet (a1,b1,c)
			;and their gcv values (gcva,gcvb,gcvc)

;------Find minimum using Brent's method-------------------------------

tol=1.e-4
itermax=100                     ;maximum # of iterations
cgold=0.3819660                 ;golden ratio
zeps=1.e-10                     ;smallest fractional accuracy

;endpoints of bracket must be in ascending order
if (a1 lt c) then begin
	a2=a1 & b2=c
endif else begin
	a2=c & b2=a1
endelse

;initializing
v=b1 & w=v & x=v
e=0.				;distance moved step before last
gcvx=gcvfunction(x,flux,basekernel,A,H,data,kernel,num=num)
gcvv=gcvx & gcvw=gcvx

;main loop of program
for i=0,itermax do begin
	xm=0.5*(a2+b2)
	tol1=tol*abs(x)+zeps & tol2=2.*tol1
	if (abs(x-xm) le (tol2-0.5*(b2-a2))) then goto, minfound  ;optimal
							;found--exit loop

	if (abs(e) gt tol1) then begin		;construct trial parabolic fit
		r=(x-w)*(gcvx-gcvv) & q=(x-v)*(gcvx-gcvw)
		p=(x-v)*q - (x-w)*r & q=2.*(q-r)
		if (q gt 0.) then p=-p
		q=abs(q)& etemp=e & e=d

	;determine acceptability of parabolic fit--if not o.k. go to golden
	;section step
		if (abs(p) ge abs(0.5*q*etemp)) or (p le q*(a2-x)) or $
			(p ge q*(b2-x)) then goto,golden
		d=p/q		;parabolic fit was o.k., take parabolic step
		u=x+d
		if (u-a2 lt tol2) or (b2-u lt tol2) then begin
			if (xm-x ge 0.) then d=abs(tol1) else d=-abs(tol1)
		endif
		goto,skipgold	;skip over golden section step
	endif


	golden:		;take golden section step, as larger of 2 segments
	if (x ge xm) then e=a2-x else e=b2-x
	d=cgold*e


	skipgold:                    	;arrive here with d from parabolic
					;fit or from golden section
	if (abs(d) ge  tol1) then begin
		u=x+d
	endif else begin
		if (d ge 0.) then u=x+abs(tol1) else u=x-abs(tol1)
	endelse

	gcvu=gcvfunction(u,flux,basekernel,A,H,data,kernel,num=num)
	;housekeeping
	if (gcvu le gcvx) then begin
		if (u ge x) then a2=x else b2=x
		v=w & gcvv=gcvw & w=x & gcvw=gcvx & x=u & gcvx=gcvu
	endif else begin
		if (u lt x) then a2=u else b2=u
		if (gcvu le gcvw) or (w eq x) then begin
			v=w & gcvv=gcvw & w=u & gcvw=gcvu
		endif else begin
			if (gcvu le gcvv) or (v eq x) or (v eq w) then begin
				v=u & gcvv=gcvu
			endif
		endelse
	endelse
	;return for another iteration

endfor
print, 'brent routine exceeded maximum iterations'
minfound:
optq=x & brent=gcvx	;return best values for optimal q and gcv minimum

;.....................................................................

return, optq

end
