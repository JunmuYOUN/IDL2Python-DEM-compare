;+
;
; NAME:
; dem_sites
;
; PURPOSE:
; Given (i) a set of input measurements, (ii) associated measurement uncertainties, 
; (iii) temperature response curves, (iv) response relative uncertainties and 
; (v) the temperature bin widths, this code implements the SITES inversion described in 
; https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract and returns the estimated differential emission measure (DEM). 
; This code processes multi-channel measurements for a single spatial pixel. For multiple pixels, 
; user should implement a loop that records the output DEM from this function at each pixel.
; See also Grid-SITES paper and code for efficient processing of large number of pixels.
;
; CALLING SEQUENCE:
; dem=dem_sites(obs_in,err_in,response,response_err,delta_temp[,convergence=convergence, $
;               obsmod=obsmod,demerr=demerr,irep=irep,ker=ker,res2=res2,totres2=totres2,gres=gres])
;
; INPUTS:
;   obs_in = multi-channel measurement. E.g. for 6-channel AIA/SDO measurement, this would be a 6-element vector in DN/s units.
;               The order of the channels is not important, as long as the same ordering is used
;               for all relevant inputs (i.e. response functions and errors).
;   err_in = the errors or uncertainties associated with the measurements, 
;               so vector same length (and in same channel order) as obs_in. 
;               These are absolute errors, not relative. Same units as obs_in (e.g. DN/s) 
;   response = the temperature response function for each channel, e.g.  as provided by the AIA/SDO Solarsoft routines.
;               E.g. for the 6-channel AIA/SDO example, this would be a [N-temperature,N-channel] array e.g. [43,6]
;   response_err = the RELATIVE uncertainties of the temperature response functions, one uncertainty for each channel. 
;               For the 6-channel AIA example, response_err will be a 6-element vector. Fractional
;               uncertainties (not percentage!). 
;               I use, for channels [94,131,171,193,211,304,335], uncertainties of [0.5,0.5,0.25,0.25,0.25,0.5,0.25]
;   delta_temp = the width of each temperature bin, in Kelvin. From trial and error, this code works best on 
;               fixed logarithmic temperature bin width, so that delta_temp in Kelvin increases with temperature
;               (bins are wider with increasing temperature). This is a vector with temperature
;               bins corresponding to the input response functions, e.g [43 bins]
;
; OUTPUTS:
;   This function returns the DEM at each temperature bin
;   
; OPTIONAL KEYWORD
;   convergence = sets the convergence threshold, default 1.e-2 (1%)
;
; OPTIONAL INPUTS
;   ker = user can provide their own smoothing kernel via this keyword. Default is Gaussian in logarithmic
;         temperature, with width set by the number of temperature bins as found by testing and described
;         in https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract. 
;   res2 = the response functions multiplied by the width of the temperature bins, required for the
;         inversion. If user is processing large number of pixels, this variable can be returned
;         from the first call to DEM_SITES at the first pixel, then passed back to DEM_SITES as a keyword
;         input for subsequent pixels, giving a small increased efficiency.
;   totres2 = as optional input/output res2, this variable can be returned
;         from the first call to DEM_SITES at the first pixel, then passed back to DEM_SITES as a keyword
;         input for subsequent pixels, giving a small increased efficiency.
;            
; OPTIONAL OUTPUTS
;   obsmod = returns the model measurement, so user can compare the input and model measurements
;   demerr = gives the (absolute) errors for the output DEM
;   irep = returns the number of iterations
;   res2 = the response functions multiplied by the width of the temperature bins, required for the
;         inversion. See res2 in optional inputs above.
;   totres2 = as optional input/output res2, see totres2 in optional inputs above
;   gres = returns the relative response functions weighted by the measurement and response errors 
;         (see equation 2 of https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract).
;         Is provided as an optional output just for user's interest, useful to see the relative
;         responses of each channel as a function of temperature (e.g. figure 1b of SITES paper).
;
; PROCEDURE:
; See https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract
;
; USE & PERMISSIONS
; If you use this code for DEM inversions, please cite https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract
; Any problems/queries, or suggestions for improvements, please email Huw Morgan, hmorgan@aber.ac.uk
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC Studentship to Aberystwyth University (Pickering)
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;  
; MODIFICATION HISTORY:
; Created at Aberystwyth University 07/2019 - Huw Morgan hmorgan@aber.ac.uk
; 
; 
;-

function dem_sites,obs_in,err_in,response,response_err,delta_temp,convergence=convergence, $
        obsmod=obsmod,demerr=demerr,irep=irep,ker=ker,res2=res2,totres2=totres2,gres=gres

if n_elements(convergence) eq 0 then convergence=1.e-2

sz=size(response)
nt=sz[1] & nwl=sz[2]

wt=1/sqrt((err_in/obs_in)^2+response_err^2)
wwt=rebin(reform(wt,1,nwl),nt,nwl)
gres=response*wwt/rebin(total(response*wwt,2),nt,nwl)

obs=obs_in

nker=(0.08*nt)>0.5
if ~keyword_set(ker) then ker=gaussian_function(nker,/norm)

irep=0l
obsprev=obs_in
demmain=fltarr(nt)
if ~keyword_set(res2) then res2=response*rebin(delta_temp,nt,nwl)
if ~keyword_set(totres2) then totres2=total(res2^2,1)
res3=res2*gres

if arg_present(demerr) then begin
  demerr=rebin(reform((err_in/obs_in)^2+response_err^2,1,nwl),nt,nwl)*gres
  demerr=sqrt(total(demerr,2))
endif

totwt=total(wt)
obswt=wt/obs_in
convobsprev=1
repeat begin
  dem=rebin(reform(obs/totres2,1,nwl),nt,nwl)*res3
  demmain=(temporary(demmain)+convol(total(dem,2),ker,edge_ZERO=1,edge_trunc=0))>0
  obs=obs_in-total(rebin(demmain,nt,nwl)*res2,1)
  convobs=total(abs(obs*obswt))/totwt
  irep++
endrep until irep gt 300 or (convobs lt convergence)

if arg_present(obsmod) then obsmod=total(rebin(demmain*delta_temp,nt,nwl)*response,1)
if arg_present(demerr) then demerr=demerr*demmain

return,demmain

end