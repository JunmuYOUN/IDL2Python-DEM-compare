; This is a simple script which computes DEMs from a set of AIA image files.
; After it is finished (which may take a while depending on the resolution of the images),
; The DEMs are stored in memory in a structure called dem_struc (see readme and the comments
; for firdem.pro for a description of the structure's contents and how to work with it).
.compile ../estimate_quantile.pro
.compile ../firdem.pro
.compile aia_firdem_wrapper.pro

; This search string should be such that findfile finds the set of 6 AIA images to invert:
if(n_elements(search_string) eq 0) then search_string = "../../aia_test_images/aia*.fits" 
filenames=findfile(search_string)
indices_used = [0,1,2,3,4,5]
aia_prep,filenames,indices_used,indout,datout,/uncomp_delete
print,filenames(indices_used)
channels='A'+string(format='(%"%d")',indout.wavelnth)
; Computing DEMs for a full-resolution AIA image takes a while (~ 1 hour). For testing,
; we rebin the data to 1024x1024. Comment out this line for full resolution:
datout=rebin(datout,1024,1024,n_elements(channels))
aia_firdem_wrapper,channels,indout.exptime,datout,dem_struc
