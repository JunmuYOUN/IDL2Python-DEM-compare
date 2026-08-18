; This is a simple, direct script to compute the DEM with the
; included example files. There is also an example wrapper which
; provides most of the same functionality (sans loading files).
; Load example files. These are L1.5 prepped files.
files = ['AIA20180203_153211_0094.fits','AIA20180203_153206_0131.fits', $
		'AIA20180203_153209_0171.fits','AIA20180203_153204_0193.fits', $
		'AIA20180203_153209_0211.fits','AIA20180203_153212_0335.fits']
read_sdo,files[0],index,data
nx = n_elements(data[*,0]) ; Get dimension of images
ny = n_elements(data[0,*])
nchan = n_elements(files)
dn = dblarr(nx,ny,nchan)
exptimes = dblarr(nchan)
for i=0,n_elements(files)-1 do begin &$
	read_sdo,files[i],index,data &$
	dn[*,*,i] = data[0:nx-1,0:ny-1] &$
	exptimes[i] = index.exptime &$
	print,index.xcen,index.ycen &$
endfor

; Inputs -- dn: array of the set of 6 images, one from each of the 
; optically thin AIA channels, that you want to use to compute the DEMs. 
; Same order as in aia_get_response -- i.e., increasing wavelength
; (94, 131, 171, 193, 211, 335). Cube should be numerical
; array of dimensions nx*ny*nchan (nchan=6 in this case). Assign
; the exposure times to 'exptimes', in same order.

logtmin = 5.5 ; Temperature range
logtmax = 7.0
; Get AIA temp response
tr_struct = aia_get_response(/temp,/dn,timedepend_date=index.date,/evenorm)
; Omit 304 angstrom channel and limit temperature range
ichan = where(tr_struct.channels ne "A304") ; omit 304 angstrom channel
ilogt = where(tr_struct.logte ge logtmin and tr_struct.logte le logtmax)
nt0 = n_elements(ilogt)
channels = tr_struct.channels[ichan]
tresp = tr_struct.all(ilogt#(1+lonarr(nchan)),(1+lonarr(nt0))#ichan)
logt=tr_struct.logte[ilogt]

; Compute errors using aia_bp_estimate_error
errs = dblarr(nx,ny,nchan)
for i=0,nchan-1 do errs[*,*,i] = AIA_BP_ESTIMATE_ERROR(dn[*,*,i] > 0,replicate(strmid(channels[i],1),nx,ny))

; Compute DEMs
dems = simple_reg_dem(dn, errs, exptimes, logt, tresp, chi2)
