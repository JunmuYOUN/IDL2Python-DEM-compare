;+
; Name: firdem_rgbmap_triangle
;
; A simple function which maps a number to a 3-dimensional vector (e.g., an RGB
; color) whose components sum to unity. The mapping is such that red is 1.0 for xin <= xlo
; (while green and blue are zero outside that range) and  decreases linearly from xlo to
; 0.5*(xlo+xhi), after which it is zero. Green is zero for xin < xlo or xin > xhi, and increases
; linearly to 1.0 as xin approaches 0.5*(xlo+xhi). Blue is zero for xin < 0.5*(xlo+xhi) and
; increases linearly to 1.0 as xin goes from that value to xhi, after which it remains at 1.0.
;
; Inputs:
;	xin: The input value at which the vector will be evaluated
;	xlo: The lower end of the range (below which red will be 1.0 and the other colors 0)
;	xhi: The upper end of the range (above which blue will be 1.0 and the other colors 0)
;
; Joseph Plowman (plowman@physics.montana.edu) 09-17-12
;-
function firdem_rgbmap_triangle,xin,xlo,xhi

	xrange = xhi-xlo

	xmid = xlo+0.5*xrange

	r=1.0-(xin-xlo)/(0.5*xrange)
	if(r lt 0.0) then r = 0.0
	if(r gt 1.0) then r = 1.0
	g=1.0-abs(xin-xmid)/(0.5*xrange)
	if(g lt 0.0) then g = 0.0
	b=1.0-abs(xhi-xin)/(0.5*xrange)
	if(b lt 0.0) then b = 0.0
	if(b gt 1.0) then b = 1.0

	return,[r,g,b]
	
end
