;+
; Name: dist_median
;
; Computes the median of a tabulated distribution function.
;
; Inputs:
;	
;	x:	Abscissa values for the function. Must be sorted in increasing order.
;	y:	The values of the function evaluated at the points in x. The function is assumed to be
;		weighted per the units of x, as well. Normalization is not required; the median returned
;		will be the point where the cumulative density is 50% of the integral over the entire
;		interval.
;
; Optional output:
;
;	tot: The integral of y over the range.
;
; Returns the median value of y integrated over x
;
; Joseph Plowman (plowman@physics.montana.edu) 09-13-12
;-
function dist_median,x,y,tot=tot

	; Evaluate the cumulative integral via a simple quadrature method:
	np = n_elements(x)
	dx0 = x(1:np-1)-x(0:np-2)
	yavg = 0.5*(y(1:np-1)+y(0:np-2))		
	integral0 = total(dx0*yavg,/cumulative)
	integrals = dblarr(np)
	integrals(1:np-1)=integral0
	
	tot=integrals(np-1)

	; Find the 50% value via interpolation and return it:
	return,interpol(x,integrals/integrals(np-1),0.5)
	
end
