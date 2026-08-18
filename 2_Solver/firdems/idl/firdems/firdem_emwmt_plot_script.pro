; This is a simple script which plots true-color images from the DEM using the emission
; weighted median temperature (EMWMT) as hue and total EM as intensity. This script is
; intended to be run after DEMs are computed from a set of images - by calling
; aia_firdem_fromfits_script.pro, for example.
.compile ../dist_median.pro
.compile ../firdem_rgbmap_triangle.pro
.compile ../firdem_emwtemp_image_color.pro
.compile ../firdem_imagefile_plot.pro

logtmin = 5.85 ; Median temperatures below logtmin will show as 100% red
logtmax = 6.75 ; Median temperatures above logtmax will show as 100% blue

if(n_elements(emwmtplotname) eq 0) then emwmtplotname = "emwmt_test.tiff"

firdem_emwtemp_image_color, dem_struc, logtmin, logtmax, image
firdem_imagefile_plot,image,emwmtplotname,/sbtiff
