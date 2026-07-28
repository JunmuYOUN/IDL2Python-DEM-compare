; Author: Mark Cheung
; Purpose: Sample script for running sparse DEM inversion (see Cheung
; et al. 2015) for details.
; Revision history: 2015-10-20 First version
;                   2015-10-23 Added note about lgT axis
;                   2017-03-25 For initialization, use basis function
;                   sets [Dirac,sigma=0.1,sigma=0.2] instead of the
;                   default, which would also include sigma=0.6.

files = file_list('.','*genx')
aiadatafile = files[0]
; Restore some prepared AIA level 1.5 data (i.e. already been aia_preped)
restgen, file=aiadatafile, struct=s

; Initialize solver. 
; This step builds up the response functions (which are
; time-dependent) and the basis functions.
; If running inversions over a set of data spanning over a few hours
; or even days, it is not necessary to re-initialize

;IF (0) THEN BEGIN 
; As discussed in the Appendix of Cheung et al. (2015), the inversion
; can push some EM into temperature bins if the lgT axis goes below
; logT ~ 5.7 (see Fig. 18). It is suggested the user check the
; dependence of the solution by varying lgTmin.
lgTmin = 5.7   ; minimum for lgT axis for inversion 
dlgT   = 0.1   ; width of lgT bin
nlgT   = 21    ; number of lgT bins
aia_sparse_em_init, timedepend = s.oindex[0].date_obs, /evenorm, use_lgtaxis=findgen(nlgT)*dlgT+lgTmin, bases_sigmas=[0.0,0.1,0.2]
lgtaxis = aia_sparse_em_lgtaxis()
; ENDIF

; We will use the data in s.img to invert. s.img is a stack of level
; 1.5, exposure time normalized AIA pixel values. 
exptimestr = '['+strjoin(string(s.oindex.exptime,format='(F8.5)'),',')+']'

; This defines the tolerance function for the inversion. 
; y denotes the values in DN/s/pixel (i.e. level 1.5, exposure time
; normalized)

; aia_bp_estimate_error gives the estimated uncertainty for a given
; DN / pixel for a given AIA channel. Since y contains DN/pixel/s,
; we have to multiply by exposure time first to pass
; aia_bp_estimate_error. Then we will divide the output of
; aia_bp_estimate_error by the exposure time again.
; If one wants to include uncertainties in atomic data, suggest to add
; the /temp keyword to aia_bp_estimate_error()
tolfunc = 'aia_bp_estimate_error(y*'+exptimestr+', [94,131,171,193,211,335], num_images='+strtrim(string(s.binning^2,format='(I4)'),2)+')/'+exptimestr

; Do DEM solve. 
; Note!!! The solver assumes the third dimension is arranged according to [94,131,171,193,211,335]
aia_sparse_em_solve, s.img, tolfunc=tolfunc, tolfac=1.4, oem=emcube, status=status, coeff=coeff
; emcube contains the emission measure contained in each lgt bin
; status contains a mask indicating whether a solution was
; found. status = 0 means a solution was found within the given constriants.
; If you want to actual coefficients of the basis functions
; (i.e. x_i's of Kx = y), then add coeff=coeff to the call to aia_sparse_em_solve

; Now show results
dispimage= aia_sparse_em_display(em=emcube)
window, xs=1440, ys=640
tv, dispimage, /true
end
