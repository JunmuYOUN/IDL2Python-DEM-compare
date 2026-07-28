; This script runs the firdem code on synthetic AIA data made from test DEMS composed
; of the sum of several random log-normal distributions, whose statistical properties are
; described by npeaks, emrange, sigrange, and trange, below. Comment out the 'delvar,seed'
; lines to obtain a different set of random peaks with each run of the code.
delvar,seed ; Ensures that we obtain the same set of random peaks each time...
.compile ../estimate_quantile.pro
.compile ../firdem.pro
.compile aia_firdem_wrapper.pro
.compile ../firdem_get_pixel.pro
.compile aia_firdem_test_4.pro
set_plot,'ps'
!p.font=0
device,file='plots/aia_dem_multimod_test_plots.eps',/encapsulated,/times
device,/inches,xsize=8.0,ysize=4.5
;set_plot,'x'
;ni=1000
;npeaks=1
;emrange=[5.0e31,5.0e31]
;sigrange=[4.0,4.0]
;trange=[6.5,6.5]
npeaks=5
emrange=[5.0e27,5.0e28]
sigrange=[0.05,0.15]
trange=[5.75,7.0]
;window,0
tempsin = [5e5,1e6,5e6,1e7]
npx=3
npy=3
!p.multi = [0,npx,npy]
for i=0,npx*npy-1 do aia_firdem_test_4, npeaks, emrange, sigrange, trange, temps_in=temps_in, $
		dem_in=dem_in, temps_out=temps_out, dem_out=dem_out, niter=ni, tr_struct=tr_struct, $
		its=its, chi2s=chi2s, seed=seed
device,/close
set_plot,'x'
delvar,seed

