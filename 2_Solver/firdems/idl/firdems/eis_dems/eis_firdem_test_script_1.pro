; This code computes and plots the response of the firdem code to a log-normal DEM at
; several temperatures, using synthetic EIS data:
.compile ../estimate_quantile.pro
.compile ../firdem_get_pixel.pro
.compile ../firdem.pro
.compile eis_firdem_wrapper.pro
.compile eis_firdem_test_1.pro
set_plot,'ps'
!p.font=0
device,file="plots/eis_dem_iterative_test_plots.eps",/encapsulated,/times
device,/inches,xsize=10.0,ysize=6
!p.multi = [0,2,2]
;window,0
tempin=5e5
eis_firdem_test_1, tempin, temps_in, dem_in, temps_out, dem_out, niter=niter
tempin=10e5
eis_firdem_test_1, tempin, temps_in, dem_in, temps_out, dem_out, niter=niter
tempin=50e5
eis_firdem_test_1, tempin, temps_in, dem_in, temps_out, dem_out, niter=niter
tempin=100e5
eis_firdem_test_1, tempin, temps_in, dem_in, temps_out, dem_out, niter=niter
;device,/close
!p.multi=[0,0,0]
;window,1
;device,file="eis_response_basis_plot_1.ps",/encapsulated
;plot_basis,dem_out
device,/close
set_plot,'x'
