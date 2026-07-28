; This script tests the fast DEM code on example quiet sun, active region, and flare DEMs
; provided with CHIANTI.
.compile ../estimate_quantile.pro
.compile ../firdem.pro
.compile aia_firdem_wrapper.pro
.compile ../firdem_get_pixel.pro
.compile aia_firdem_curve_test.pro

;set_plot,'x'
;device,set_character_size=[12,15]
set_plot,'ps'
device,file='plots/aia_dem_chianti_test_plots.eps',/encapsulated
device,/inches,xsize=9,ysize=3

!p.multi = [0,3,1]

filename='/ssw/packages/chianti/dbase/dem/quiet_sun.dem'
qsdem_struc = read_ascii(filename)
qstemps = qsdem_struc.field1(0,*)
validtemps = where(finite(qstemps) and qstemps ge 0)
qstemps = qstemps(validtemps)
qsdem = 10^qsdem_struc.field1(1,validtemps)*10^qstemps*alog(10)

aia_firdem_curve_test, qstemps, qsdem, qstemps_out, qsdem_out, niter=niter, tr_struct=aia_tr_struct,its=qsits,chi2s=qschi2s, data=data,seed=seed


filename='/ssw/packages/chianti/dbase/dem/active_region.dem'
ardem_struc = read_ascii(filename)
artemps = ardem_struc.field1(0,*)
validtemps = where(finite(artemps) and artemps ge 0)
artemps = artemps(validtemps)
ardem = 10^ardem_struc.field1(1,validtemps)*10^artemps*alog(10)

aia_firdem_curve_test, artemps, ardem, artemps_out, ardem_out, niter=niter, tr_struct=aia_tr_struct,its=arits,chi2s=archi2s, data=data,seed=seed


filename='/ssw/packages/chianti/dbase/dem/flare.dem'
fldem_struc = read_ascii(filename)
fltemps = fldem_struc.field1(0,*)
validtemps = where(finite(fltemps) and fltemps ge 0)
fltemps = fltemps(validtemps)
fldem = 10^fldem_struc.field1(1,validtemps)*10^fltemps*alog(10)

aia_firdem_curve_test, fltemps, fldem, fltemps_out, fldem_out, niter=niter, tr_struct=aia_tr_struct,its=flits,chi2s=flchi2s, data=data,seed=seed

device,/close
set_plot,'x'
