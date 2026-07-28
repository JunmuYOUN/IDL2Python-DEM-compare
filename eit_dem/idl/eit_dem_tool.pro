;+
; NAME:
;       EIT_DEM_TOOL
;
;  A user guide is included in the Solar Software distribution.  For further information
;  go to your local directory $SSW/soho/eit/idl/response/dem_tool/user_guide.pdf 
;
;
; PURPOSE:
;       Create a Differential Emission Measure (DEM) Map, i.e. computed for
;       each pixel individually, based upon the 4 EIT channels. Optionally,
;       also returns a line/bandpass map or the irradiance for chosen line.
;       DEM maps can be from pre-computed daily map database
;
; CATEGORY:
;       Analysis
;
; CALLING SEQUENCE:
;       eit_dem_tool,files,dem_map,temp,[mk_line=mk_line],[line_map=line_map],$
;        [ion=ion], [wmin=wmin],[wmax=wmax],[irradiance=irradiance],$
;        [no_dem_map=no_dem_map], [date = date],[leak_284=leak_284],$
;        [coefs_only=coefs_only],[nrl=nrl]
;
; INPUTS:
;       files                    - Names of 4 EIT raw files in order 171,195,284,304
;                                    or a processed (1024,1024,4) array
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       mk_line               - set to produce output line map
;       ion                       - set to ion name for line maps
;       wmin                   - minimum wavelength to consider for line maps
;       wmax                  - maximum wavelength to consider for line maps
;       irradiance             - set if only want line irradiance
;       no_dem_map       - set if do not wish calculation of output DEM
;                                     map i.e. only wish eit_line_map output
;       date                     - read in pre-computed DEM for set date
;       leak_284             - set to use in eit_prep for the 284 A channel light leak
;       coefs_only          - set if return correction coefficients only
;       nrl	    	    	     - REQUIRED to set for FITS files from NRL archive (to run the
;                                    NRL version of eit_prep)
;
; OUTPUTS:
;       DEM_MAP        - differential emission measure map at
;                                     20 predefined temperatures
;       TEMP                 - temperatures (log)
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;       LINE_MAP        - output line map or irradiance
;
; COMMON BLOCKS: none.
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;   Modifies starter DEM to fit observed EIT data.
;   Starter DEM = CHIANTI V3 active_region for logT = 4.6 - 6.5 and
;     CHIANTI V5 quiet_sun for logT =  4.0 - 4.6
;
; SUBROUTINES:
;   E_INTERP, E_FIT2, EIT_KCORR : included in this file, automatically compiled
;   EIT_LINE_MAP : Separate procedure
;   MK_EIT_SPEC: Included in EIT_LINE_MAP
;
; MODIFICATION HISTORY:
;       Written by:  J. Newmark         October 2001
;       Modified:  J. Newmark           May 2006
;                        N. Rich                   October 2010  Add /NRL keyword, sample data
;
;-
;-----------------------------------------------------------
function e_interp, img1, img2, t1,t2,t3
; Simple linear interpolation
slope = (img2-img1)/(t2-t1)
return, slope*(t3-t1) + img1
end
;-----------------------------------------------------------
function e_fit2,x,y
; matix inversion to find array coefficents
; input either 2,3 or 5 x's independent of number of y's
; have y input as [nfits,nx]
; 2 elements solves for linear
; 3 elements solves for quadratic
; 5 elements solves for linear (2 terms) + quadratic (3 terms)
case n_elements(x) of
   2: mat = [[x(0),1],[x(1),1]]
   3: mat = [[x(0)^2,x(0),1],[x(1)^2,x(1),1],[x(2)^2,x(2),1]]
   5: mat = [[x(0),1,0,0,0],[x(1),1,0,0,0],[0,0,x(2)^2,x(2),1],$
            [0,0,x(3)^2,x(3),1],[0,0,x(4)^2,x(4),1]]
   else:
endcase
return, invert(mat,/double)##y
end
;-----------------------------------------------------------
pro eit_kcorr,eitimgs,dem,temp,x,dn_maps,iter=iter,kcorr171=kcorr171,$
      kcorr195=kcorr195,kcorr284=kcorr284,kcorr304=kcorr304

lowsig = [0.1,0.1,0.1,0.1]

if keyword_set(iter) then begin
    dn_pred_171 = (*dn_maps).dn171
    dn_pred_195 = (*dn_maps).dn195
    cool_284 = (*dn_maps).cool_284
    hot_284 = (*dn_maps).hot_284
    cool_304 = (*dn_maps).cool_304
    warm_304 = (*dn_maps).warm_304
    hot_304 = (*dn_maps).hot_304
    dn_pred_304 = cool_304 + warm_304 + hot_304
endif else begin
    dn_pred_171 = eit_line_map(171,dem,temp,x)
    dn_pred_195 = eit_line_map(195,dem,temp,x)
    cool_284 = eit_line_map(284,dem,temp,x,trange=[0,14])
    hot_284 = eit_line_map(284,dem,temp,x,trange=[14,19])
    cool_304 = eit_line_map(304,dem,temp,x,trange=[0,9])
    warm_304 = eit_line_map(304,dem,temp,x,trange=[9,15])
    hot_304 = eit_line_map(304,dem,temp,x,trange=[15,19])
    dn_pred_304 = cool_304 + warm_304 + hot_304
endelse

kcorr171 = alog10( ((*eitimgs).eit171>lowsig(0))/(dn_pred_171>lowsig(0)) > 1e-5)
kcorr195 = alog10( ((*eitimgs).eit195>lowsig(1))/(dn_pred_195>lowsig(1)) > 1e-5)
kcorr304 = alog10( ((*eitimgs).eit304>lowsig(3))/(dn_pred_304>lowsig(3)) > 1e-5)

k284_cool = e_interp(kcorr304,kcorr171,x(0),x(1),x(4))
cool_284 = (10^k284_cool) * cool_284
kcorr284 = alog10( (((*eitimgs).eit284>lowsig(2))- cool_284) / (hot_284 > 1e-5) > 1e-5)

k304_warm = e_interp(kcorr171,kcorr195,x(1),x(2),x(5))
warm_304 = (10^k304_warm) * warm_304
k304_hot = e_interp(kcorr195,kcorr284,x(2),x(3),x(6))
hot_304 = (10^k304_hot) * hot_304
hot_304 = hot_304 + warm_304
kcorr304 = alog10( (((*eitimgs).eit304>lowsig(3)) - hot_304) / cool_304 > 1e-5)

end
;-----------------------------------------------------------

pro eit_dem_tool,files,dem_map,temp,mk_line=mk_line,line_map=line_map,ion=ion,$
        wmin=wmin,wmax=wmax,irradiance=irradiance,no_dem_map=no_dem_map,$
        date=date,leak_284=leak_284,coefs_only=coefs_only,nrl=nrl

print,"_________________________________________________"
print, "A user guide is included in the Solar Soft distribution.  For further " $
+ "help go to your local directory" $
+ " $SSW/soho/eit/idl/response/dem_tool/user_guide.pdf"
print,"_________________________________________________"

; definitions
defsysv,'!xuvtop',exist=yes_chianti
if not yes_chianti then ssw_packages,/chianti
; EIT specific temperature response function peaks
x = [4.89,5.98,6.1,6.3,5.85,5.9,6.23]   ;LogT values
temp = indgen(20)*0.1+4.6
;starter DEM from CHIANTI v3 active_region.dem
dem = [23.3700, 23.0000,  22.6600,  22.3700,  22.1600,  22.0500, $
      21.9800,  21.9700,  21.9600,  21.9800,  22.0400,  22.1000, $
      22.1800,  22.2500,  22.3500,  22.4300,  22.4800,  22.4500, $
      22.3400 , 22.1500]

he_fac = 9.5    ; Assumes CHIANTI v3.0 Feldman abundance
;
if n_elements(date) eq 0 then begin
;process images if necessary
  if datatype(files) eq 'STR' then eit_prep,files,hdrs,imgs,leak_284=leak_284,nrl=nrl $
            else imgs = temporary(files)
  missing = fltarr(1024,1024)
  for i=0,3 do begin
    nm = where(imgs(*,*,i) le 0,mcnt)
    if mcnt gt 0 then missing(nm) = -99999
  endfor
  nm = where(missing lt 0,mcnt)
  eitimgs = ptr_new({eit171:imgs(*,*,0),eit195:imgs(*,*,1),eit284:imgs(*,*,2),$
                   eit304:imgs(*,*,3)})
  delvarx,imgs  ; save memory

;compute Correction Coefficient maps, one map at each Temp. peak- iterate once
  eit_kcorr,eitimgs,dem,temp,x,kcorr171=kcorr171,kcorr195=kcorr195,$
           kcorr284=kcorr284,kcorr304=kcorr304
  sz = size(kcorr171)
  sz = sz(1)*sz(2)
  y=fltarr(sz,5)
  y(*,0)=reform(kcorr304,sz)
  y(*,1)=reform(kcorr171,sz)
  y(*,2)=reform(kcorr171,sz)
  y(*,3)=reform(kcorr195,sz)
  y(*,4)=reform(kcorr284,sz)
  xt = [x(0:1),x(1),x(2:3)]
  coefs = e_fit2(xt,y)
  coefs = reform(coefs,sqrt(sz),sqrt(sz),5)

  dn_maps = ptr_new({dn171:fltarr(1024,1024),dn195:fltarr(1024,1024),$
              cool_284:fltarr(1024,1024),warm_304:fltarr(1024,1024),$
              cool_304:fltarr(1024,1024),hot_304:fltarr(1024,1024),$
              hot_284:fltarr(1024,1024)})
  (*dn_maps).dn171 = eit_line_map(171,dem,temp,x,coefs)
  (*dn_maps).dn195 = eit_line_map(195,dem,temp,x,coefs)
  (*dn_maps).cool_284 = eit_line_map(284,dem,temp,x,coefs,trange=[0,14])
  (*dn_maps).hot_284 = eit_line_map(284,dem,temp,x,coefs,trange=[14,19])
  (*dn_maps).cool_304 = eit_line_map(304,dem,temp,x,coefs,trange=[0,9])
  (*dn_maps).warm_304 = eit_line_map(304,dem,temp,x,coefs,trange=[9,15])
  (*dn_maps).hot_304 = eit_line_map(304,dem,temp,x,coefs,trange=[15,19])

  eit_kcorr,eitimgs,dem,temp,x,dn_maps,/iter,kcorr171=k171,kcorr195=k195,$
            kcorr284=k284,kcorr304=k304

  ptr_free,eitimgs  ; save memory
  ptr_free,dn_maps  ; save memory

  kcorr171 = temporary(kcorr171) + k171
  kcorr195 = temporary(kcorr195) + k195
  kcorr284 = temporary(kcorr284) + k284
  kcorr304 = temporary(kcorr304) + k304
  szi = size(kcorr171)
  sz = szi(1)*szi(2)
  y=fltarr(sz,5)
  y(*,0)=reform(kcorr304,sz)
  y(*,1)=reform(kcorr171,sz)
  y(*,2)=reform(kcorr171,sz)
  y(*,3)=reform(kcorr195,sz)
  y(*,4)=reform(kcorr284,sz)
  xt = [x(0:1),x(1),x(2:3)]
  coefs = e_fit2(xt,y)
  coefs = reform(coefs,sqrt(sz),sqrt(sz),5)
  if mcnt gt 0 then for i=0,4 do coefs(0,0,i) = coefs(*,*,i)+missing
endif else begin
; read in pre-computed coefficients
  coefs = fltarr(1024,1024,5,/nozero)
  ut = anytim2utc(date)
  dt = strtrim(ut.mjd,2)
  openr,lun,'eit_dem_'+dt+'.bin',/get_lun
  readu,lun,coefs
  free_lun,lun
  szi = size(coefs)
endelse

; compute DEM maps from coefficients - expansion to linear+quadratic fit
if not keyword_set(no_dem_map) then begin
  dem_m = fltarr(szi(1),szi(2),n_elements(temp))
  for i=0,n_elements(temp)-1 do begin
     this_t = temp(i)
     if this_t le x(1) then kcorr_t = coefs(*,*,0)*this_t + coefs(*,*,1) $
     else kcorr_t = coefs(*,*,2)*this_t*this_t+coefs(*,*,3)*this_t + $
                           coefs(*,*,4)
     if this_t eq 6.3 then k63 = kcorr_t
     if this_t eq 6.4 then kcorr_t = temporary(kcorr_t) > (k63 - 3) < (k63 + 0.6)
     if this_t eq 6.5 then kcorr_t = temporary(kcorr_t) > (k63 - 6) < (k63 + 1.)
     dem_m(*,*,i) = (dem(i) + (kcorr_t < 4)) >0
  endfor

; Add in cool continuum extension to DEM down to logT=4.0
; starter DEM = CHIANTI v5 quiet_sun.dem
; scale starter so matches at logT=4.6, CHIANTI shape not changed
  dem_scale = dem_m(*,*,0)/21.897
  dem = [26.391, 25.383, 24.469, 23.654, 22.946, 22.357]
  temp = indgen(26)*0.1+4.0 ;LogT values
  dem_map = fltarr(szi(1),szi(2),n_elements(temp),/nozero)
  dem_map[*,*,6:*] = temporary(dem_m)
  for i = 0,5 do dem_map[*,*,i] = dem(i)*dem_scale
endif
if keyword_set(coefs_only) then dem_map = temporary(coefs)

; Compute line map or irradiance
if n_elements(mk_line) ne 0 then begin
  line_map = eit_line_map(1,he_fac=he_fac,instr='dum',dem_map,temp,$
    ion=ion,wmin=wmin,wmax=wmax,irradiance=irradiance)
endif

end
