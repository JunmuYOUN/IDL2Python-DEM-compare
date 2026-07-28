.compile simple_reg_dem.pro

nx = 100 ; Size of test image (nx*ny MC random DEM inversions)
ny = 100
logtmin = 5.5 ; Min and Max (log10) temperature
logtmax = 7.5

; Get temperature response, take out 304 channel, and limit temperatures to logtmin...logtmax:
if(n_elements(tr_struct) eq 0) then tr_struct = aia_get_response(timedepend_date=timedepend_date,/temp,/dn)
ichan = where(tr_struct.channels ne "A304")
ilogt = where(tr_struct.logte ge logtmin and tr_struct.logte le logtmax)
nt0 = n_elements(ilogt)
nchan = n_elements(ichan)
channels = tr_struct.channels(ichan)
tresp = tr_struct.all(ilogt#(1+lonarr(nchan)),(1+lonarr(nt0))#ichan)
logt=tr_struct.logte[ilogt]
exptimes = 1.0+dblarr(nchan) ; 1 second exposures assumed

; Make a lognormal test input DEM:
dem_amp = 1.0e29 ; Amplitude in units cm^(-5) per unit log10(T)
dem_cen = 6.5 ; DEM center in log10(T)
dem_sig = 0.1 ; half-width (actually lognormal sigma)
dem_in = dem_amp*exp(-0.5*(logt-dem_cen)^2.0/dem_sig^2.0)

; Get DN for input DEM in each channel:
dn0 = dblarr(nchan)
for i=0,nchan-1 do dn0[i] = int_tabulated(logt,dem_in*tresp[*,i])

; Get errors for input DEM (stripping leading 'A' from channel with strmid):
errchans = strmid(channels,1)
errs0 = AIA_BP_ESTIMATE_ERROR(dn0, errchans)

; propogate dn0 onto array (assume gaussian errors for simplicity):
dn = dblarr(nx,ny,nchan)
errs = dblarr(nx,ny,nchan)
for i=0,nchan-1 do dn[*,*,i] = dn0[i]+errs0[i]*randomn(seed,nx,ny)
for i=0,nchan-1 do errs[*,*,i] = AIA_BP_ESTIMATE_ERROR(dn[*,*,i] > 0,replicate(errchans[i],nx,ny))

; Compute DEM. Replace the synthetic AIA dn with your own data to compute DEM for that data. Make sure to 
; update the exposure times, and set timedepend_date to the time of your observations:
dems_out = simple_reg_dem(dn, errs, exptimes, logt, tresp, chi2)

plot,logt,dem_in,xtitle='Log!D10!N(T) (Kelvin)', ytitle='DEM (cm!U-5!N per unit Log!D10!N(T)'
oplot,logt,dems_out[0,0,*],linestyle=1
ssw_legend,['Input test DEM','AIA DEM inversion'],linestyle=[0,1]