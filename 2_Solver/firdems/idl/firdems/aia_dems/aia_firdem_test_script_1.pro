; This code computes and plots the response the firdem code to a log-normal distribution at
; several temperatures, using AIA data.
.compile ../estimate_quantile.pro
.compile ../firdem.pro
.compile aia_firdem_wrapper.pro
.compile ../firdem_get_pixel.pro
.compile aia_firdem_test_1.pro
set_plot,'ps'
!p.font=0
device,file='plots/aia_dem_iterative_test_plots.eps',/encapsulated,/times
device,/inches,xsize=10.0,ysize=6
!p.multi = [0,2,2]

tempsin = [5e5,1e6,5e6,1e7]
tempin=5e5
aia_firdem_test_1, tempsin[0], temps_in, dem_in, temps_out, dem_out, niter=ni, tr_struct=tr_struct,its=its1,chi2s=chi2s1
tempin=10e5
aia_firdem_test_1, tempsin[1], temps_in, dem_in, temps_out, dem_out, niter=ni, tr_struct=tr_struct,its=its2,chi2s=chi2s2
tempin=50e5
aia_firdem_test_1, tempsin[2], temps_in, dem_in, temps_out, dem_out, niter=ni, tr_struct=tr_struct,its=its3,chi2s=chi2s3
tempin=100e5
aia_firdem_test_1, tempsin[3], temps_in, dem_in, temps_out, dem_out, niter=ni, tr_struct=tr_struct,its=its4,chi2s=chi2s4
device,/close
set_plot,'x'

;window,1
;device,file='plots/aia_dem_iterative_test_chi2_plots.eps',/encapsulated
!p.multi = [0,2,2]
hist = histogram(chi2s1,nbins=40,locations=locs,min=0,max=10)
plot,locs,hist,psym=10,title=string("Chi Squared histogram for T=",tempsin[0]),xtitle="Reduced Chi Squared",ytitle = "Counts"
hist = histogram(chi2s2,nbins=40,locations=locs,min=0,max=10)
plot,locs,hist,psym=10,title=string("Chi Squared histogram for T=",tempsin[1]),xtitle="Reduced Chi Squared",ytitle = "Counts"
hist = histogram(chi2s3,nbins=40,locations=locs,min=0,max=10)
plot,locs,hist,psym=10,title=string("Chi Squared histogram for T=",tempsin[2]),xtitle="Reduced Chi Squared",ytitle = "Counts"
hist = histogram(chi2s4,nbins=40,locations=locs,min=0,max=10)
plot,locs,hist,psym=10,title=string("Chi Squared histogram for T=",tempsin[3]),xtitle="Reduced Chi Squared",ytitle = "Counts"
;device,/close
set_plot,'x'
!p.multi=0
