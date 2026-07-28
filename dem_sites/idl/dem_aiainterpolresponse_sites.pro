;+
;
; NAME:
; dem_aiainterpolresponse_sites
;
; PURPOSE:
; Given user choice of channels (specified by wavelength), number of temperature bins, a logarithmic
; temperature range (min,max) and a temperature response structure returned by aia_get_response.pro, 
; this procedure takes the temperature response profiles within the structure and interpolates, for 
; each channel specified by user, to the required temperature range and number of bins. It also returns
; the uncertainty for each channel's response function. User can choose to have a regular logarithmic
; bin otherwise the bins are regular in non-logarithmic (linear) temperature 
;
; CALLING SEQUENCE:
; dem_aiainterpolresponse_sites,wavel,nout,logtempra,response_structure,logte_out, $
;            response_out,reserr[,reglog=reglog,dtemp=dtemp]
;
; INPUTS:
;   wavel = user's choice of required wavelength channels and their order e.g. [94,131,171,193,211,304,335]
;   nout = required number of output temperature bins e.g. 41
;   logtempra = two-element vector giving [min,max] of required logarithmic temperature range e.g. [5.7,7.0]
;   response_structure = the structure returned by aia_get_response function 
;
; OUTPUTS:
;   logte_out = the logarithmic temperature at each bin of the output response, vector of size Nout
;   response_out = the output temperature response curves, of size [Nout, Nchannels]
;   reserr = the estimated relative error of each temperature response curve, vector of size Nchannels
;   
;
; OPTIONAL KEYWORD
;   reglog = flag, if set (/reglog), sets a regular logarithmic temperature scale, else regular in
;             non-logarithmic space.
;   dtemp = returns the width of each temperature bin (in non-logrithmic scale, in Kelvin)
;
; PROCEDURE:
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
pro dem_aiainterpolresponse_sites,wavel,nout,logtempra,response_structure,logte_out, $
            response_out,reserr,reglog=reglog,dtemp=dtemp


;extract required temperature range and in wavelength/channel order set by user
nwl=n_elements(wavel)
indtemp=where(response_structure.logte ge logtempra[0] and response_structure.logte le logtempra[1],ntemp, $
                  comp=nindtemp,ncomp=ncomp)
ntin=n_elements(response_structure.logte)
dt0=deriv(dindgen(ntin),10^response_structure.logte)
logte_in=response_structure.logte[indtemp]
response_in=dblarr(ntemp,nwl)
reserr_in=dblarr(nwl)
wlall=long(strmid(response_structure.channels,1))
nwl=n_elements(wavel)
for iwl=0,nwl-1 do begin
  void=min(abs(wlall-wavel[iwl]),indwl)
  response_in[*,iwl]=response_structure.all[indtemp,indwl]
  if ncomp gt 0 then $
      reserr_in[iwl]=total(response_structure.all[nindtemp,indwl]*dt0[nindtemp])/ $
                       total(response_structure.all[indtemp,indwl]*dt0[indtemp])
endfor


wl0=[94,131,171,193,211,304,335]
ru0=[0.5,0.5,0.25,0.25,0.25,0.5,0.25]
resuncert=fltarr(nwl)
for i=0,nwl-1 do begin
  ind=where(wavel[i] eq wl0)
  resuncert[i]=ru0[ind]
endfor
reserr=sqrt(resuncert^2+reserr_in^2)

;and interpolate to required number of temperature bins
logte_out=keyword_set(reglog)? $
   interpol([min(logte_in),max(logte_in)],nout): $
    alog10(interpol([min(10^logte_in),max(10^logte_in)],nout))
response_out=dblarr(nout,nwl)
ntin=n_elements(logte_in)
dt=deriv(dindgen(ntin),10^logte_in)
ntout=n_elements(logte_out)
dt2=deriv(dindgen(ntout),10^logte_out)
for i=0,nwl-1 do $
  response_out[*,i]=interpol(response_in[*,i]/dt,logte_in,logte_out,/spline)*dt2

dtemp=dt2

end