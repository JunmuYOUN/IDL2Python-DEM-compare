; This is a simple script which plots true-color images from the DEM using 3 'virtual' channels.
; The channels are made by convolving the DEM with 3 log-normal temperature response functions
; whose centers and sigmas are determined by tr_centers and tr_sigs, below. This script is
; intended to be run after DEMs are computed from a set of images - by calling
; aia_firdem_fromfits_script.pro, for example.
.compile ../firdem_vchan_image.pro
.compile ../firdem_imagefile_plot.pro

logtmin=5.5
logtmax=7.5
nt=200
logt = logtmin+(logtmax-logtmin)*findgen(nt)/(nt-1.0)
tr_centers = [5.85,6.3,6.75]
tr_sigs = [0.24,0.08,0.24]
;tr_centers = [5.9,6.4,7.0]
;tr_sigs = [0.3,0.1,0.4]
tresps = dblarr(3,nt)

pi=2.0*acos(0)

for i=0,2 do tresps(i,*) = exp(-(tr_centers(i)-logt)^2/(2.0*tr_sigs(i)^2))/sqrt(2*pi);*tr_sigs(i)^2)

firdem_vchan_image, dem_struc, image, vtr_logt=logt, vtrs=tresps
firdem_imagefile_plot,image,"vchan_test.tiff",/sbtiff
