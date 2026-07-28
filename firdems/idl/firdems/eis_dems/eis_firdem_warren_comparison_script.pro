; This script computes the DEM using EIS data from Warren et al 2011 (ApJ 734, 90) and plots
; it:
.compile ../estimate_quantile.pro
.compile ../firdem_get_pixel.pro
.compile ../firdem.pro
.compile eis_firdem_wrapper.pro
.compile eis_firdem_warren_comparison.pro

;set_plot,'x'
!p.multi=0
set_plot,'ps'
!p.font=0
device,file='plots/eis_dem_warren_comparison.eps',/encapsulated,/times
device,/inches,xsize=7.2,ysize=4
eis_firdem_warren_comparison,temp,dem
device,/close
set_plot,'x'

