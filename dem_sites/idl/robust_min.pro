;+
;
; NAME:
; robust_min
;
; PURPOSE:
; Calculates the robust minimum of an input array
;
; CALLING SEQUENCE:
; min=robust_min(array[,percent,max=max])
;
; INPUTS:
;   array=input array, numerical
;
; OPTIONAL INPUT:
;   percent=the percentile minimum to use, default is 1%
;
; OUTPUTS:
;   function returns robust minimum
;
;
; OPTIONAL KEYWORD
;   max = returns robust maximum of the array at same percentile
;
; PROCEDURE:
;   Sorts array members into ascending order, then interpolates to the value
;   where number of array members is equal to the percentile
;
; USE & PERMISSIONS
;  If you reuse in your own code, please include acknowledgment to Huw Morgan (see below)
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 07/2019 - Huw Morgan hmorgan@aber.ac.uk
;
;
;-
function robust_min,y0,per,max=max;,nbin=nbin

if n_params() lt 2 then per=1

max=!values.f_nan
indok=where(finite(y0),n)
if n eq 0 then return,!values.f_nan
indsort=sort(y0[indok])
indmin=per*float(n-1)/100.
indmax=(100-per)*float(n-1)/100.
min=interpol(y0[indok[indsort]],findgen(n),indmin)
if arg_present(max) then $
  max=interpol(y0[indok[indsort]],findgen(n),indmax)

return,min

end