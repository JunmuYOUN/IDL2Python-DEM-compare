;+
;
; NAME:
; dem_openaiafiles_sites
;
; PURPOSE:
; Reads in AIA files and returns datacube of dimensions [nx,ny,nchannels]. This routine is useful
; if user wants to combine multiple exposures, for example combining 3 consecutive observations 
; across 6 channels. Also returns estimate of error based on Poisson and other measurement errors.
;
; CALLING SEQUENCE:
; dem_openaiafiles_sites,files,out,hdrmain,hdrsave,noise_estimate,npixuser=npixuser, $
;                   clip=clip,synoptic=synoptic,maxht=maxht
;
; INPUTS:
;   files = user can provide variable files as:
;     (i) string array e.g. 6 filenames in array corresponding to 6 channels taken at same observation date
;     (ii) Pointer array e.g. 6 pointers corresponding to 6 channels, each pointer points to a string or string array
;     containing one or more filenames. If a pointer contains multiple filenames for a channel,
;     this procedure will combine these observations by taking the mean.
;     In either case, output will be a datacube of [nx,ny,nchannels]
;
; OUTPUTS:
;   out = output datacube in units DN/s of size [nx,ny,nchannels]
;   hdrmain = pointer array of size [nchannels]. Each pointer contains that channel's header for each individual file
;   hdrsave = a convenient master header for the datacube. The observation time is set as the mean observation time
;             over all files contained in datacube.
;   noise_estimate = estimate of noise (in DN/s)
;
; OPTIONAL KEYWORDS
;   npixuser = required size of full image (prior to any region clipping)
;               E.g. if npixuser set to 1024 then full-res 4096x4096 images will be
;               first rebinned to 1024x1024. Any region clipping (keyword clip) is applied after.
;               Important that npixuser is integer multiple or factor of AIA image e.g. 2048, 1024, 512
;               If npixuser is set to 1024, and no clipping is set, then output datacube is [1024,1024,nchannels]
;   clip = four-element vector that defines required output region [xleft,ybottom,xright,ytop] in pixels. 
;           The output datacube and noise cube are clipped to this region.
;   synoptic = if set, then the files are AIA synoptic data
;   maxht = ;sets the field of view's maximum heliocentric distance in solar radii. Pixels above this
;             distance set to NAN. Default 1.15Rs.
;   clean = replace pixels with zero values with local median. Effective, but slow
;
; OPTIONAL INPUTS
;   none
;   
; OPTIONAL OUTPUTS
;   none
;
; PROCEDURE:
; Users aia_prep for full-res images, read_sdo for synoptic. Uses Boerner's routines for noise estimates
; (e.g. Boerner et al (2012, Sol Phys)). 
; 
; EXAMPLE:
  ; This example reads, for multiple channels, several files for each channel. So e.g. 6 files all taken
  ; consecutively over ~1 minute for 171, same for other channels. Can have different number of files for 
  ; different channels since we use pointers.
  ; 
  ; ;read in channel files from synoptic data directory
  ; directory='path_to_your_AIA_files'
  ; fileptr=ptrarr(nchannels)
  ; ;repeat following for all required channels
  ; f171=file_search(directory+'/AIA*_0171.fits')
  ; fileptr[0]=ptr_new(f171)
  ; 
  ; ;synoptic images are 1024x1024. Let's rebin to 512
  ; npixuser=512
  ; 
  ; ;define region near disk center
  ; clip=[200,200,300,350];note these are pixel coordinates in the rebinned 512x512 images
  ; 
  ; ;let's set maximum field of view height to 1.3Rs
  ; maxht=1.3
  ; 
  ; dem_openaiafiles_sites,files,out,hdrmain,hdrsave,noise_estimate,npixuser=npixuser, $
  ;                      clip=clip,/synoptic,maxht=maxht
;
; USE & PERMISSIONS
; If you use this code for DEM inversions, please cite https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract
; Any problems/queries, or suggestions for improvements, please email Huw Morgan, hmorgan@aber.ac.uk
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 07/2019 - Huw Morgan hmorgan@aber.ac.uk
;
;
;-

pro dem_openaiafiles_sites,files,out,hdrmain,hdrsave,noise_estimate,npixuser=npixuser, $
                clip=clip,synoptic=synoptic,maxht=maxht,clean=clean

;user can provide variable files as:
; (i) string array e.g. 6 filenames in array corresponding to 6 channels taken at same observation date
; (ii) Pointer array e.g. 6 pointers corresponding to 6 channels, each pointer points to a string or string array
;     containing one or more filenames. If a pointer contains multiple filenames for a channel,
;     this procedure will combine these observations by taking the mean.
; In either case, output will be a datacube of [nx,ny,nchannels]
fptr=size(files,/type) eq 10;is input files a string array or pointer array?

;nx and ny define required size of full image (prior to any region clipping) 
;E.g. if npixuser set to 1024 then full-res 4096x4096 images will be
;first rebinned to 1024x1024. Any region clipping (keyword clip) is applied after.
;Important that npixuser is integer multiple or factor of AIA image e.g. 2048, 1024, 512
nx=keyword_set(npixuser)?npixuser:1024
ny=keyword_set(npixuser)?npixuser:1024

;number of channels (Nwl=Nwavelength)
nwl=n_elements(files)

cdeltmain=0.60000002*4096/nx

;set the field of view's maximum heliocentric distance in solar radii. Pixels above this
;distance set to NAN
maxht=keyword_set(maxht)?maxht:1.15

hdrsave=-1

;read estimated errors
;DN per photon, table 2 of Boerner et al (2012, Sol Phys)
;gwl=[94,131,171,193,211,304,335]
;g_dn_phot=[2.128, 1.523, 1.168, 1.024, 0.946, 0.658, 0.596]
terr=AIA_BP_READ_ERROR_TABLE()

;set up arrays for reading in files
hdrmain=ptrarr(nwl)
tai=dblarr(nwl)
cdelt=dblarr(nwl)
crpix1=dblarr(nwl)
crpix2=dblarr(nwl)
for iwl=0,nwl-1 do begin
  
  filesnow=fptr?*files[iwl]:files[iwl]
  nfiles=n_elements(filesnow)
  
  ;if fatal problem with current file, return
  ;catch,error_status
  error_status=0
  if error_status ne 0 then begin
    print,'dem_openaiafiles_sites: error occurred ',error_status,filesnow
    print,!error_state.msg
    hdrsave=-1
    catch,/cancel
    return
  endif
  
  print,'Reading in files (patience!)'
  if keyword_set(synoptic) then $
  read_sdo,filesnow,hdr,im,/use_hdr_pnt,/UNCOMP_DELETE else $
  aia_prep,filesnow,indgen(n_elements(filesnow)),hdr,im;,/use_hdr_pnt          

  nimage=n_elements(hdr)
  exptime=hdr.exptime
  hdr=hdr[0]
  if hdr.img_type eq 'DARK' or hdr.exptime lt 0.1 then begin
    print,'Dark image,skipping (open_files_temp_maps_aia)'
    hdrsave=-1
    return
  endif

  if hdr.naxis1 ne 4096 then begin;problem with synoptic header? This not satisfactory but works for me
    hdr.cdelt1=0.60000002*4096/hdr.naxis1
    hdr.cdelt2=0.60000002*4096/hdr.naxis2
  endif

  im=float(im)
  indsat=where(im ge 1.4e4,cntsat)
  if cntsat gt 0 then im[indsat]=!values.f_nan


  if keyword_set(clean) then begin
    ;replace negative intensities with local medians
    ;increase local median width incrementally if negatives still present
    for i=0,nimage-1 do begin
      indneg=where(im[*,*,i] le 0,cntneg)
      if cntneg eq 0 then continue
      imnow=im[*,*,i]
      md=median(imnow,3)
      imnow[indneg]=md[indneg]
      indneg=where(imnow le 0,cntneg)
      im[*,*,i]=imnow
      if cntneg eq 0 then continue
      md=median(imnow,5)
      imnow[indneg]=md[indneg]
      indneg=where(imnow le 0,cntneg)
      im[*,*,i]=imnow
      if cntneg eq 0 then continue
      md=median(imnow,7)
      imnow[indneg]=md[indneg]
      im[*,*,i]=imnow
    endfor
  endif
  
  
  print,'Estimating count noise'
  ;estimate of poisson and pixel-to-pixel noise
  noisenow=fltarr(hdr.naxis1,hdr.naxis2,nimage)
  sumfactor=round(hdr.cdelt1/0.6)*round(hdr.cdelt2/0.6)
  indterr=where(terr.wavelnth eq hdr[0].wavelnth)
  for i=0,nimage-1 do begin
    dnperpht=terr[indterr].dnperpht
    ;despite chunk above to deal with negative intensities, I use
    ;the absolute value here just to avoid NANs. Not a perfect solution, but convenient
    ns=sqrt(abs(im[*,*,i]/dnperpht))*dnperpht
    noisenow[*,*,i]=ns
  endfor
  
  
  if hdr.naxis1 ne nx or hdr.naxis2 ne ny then begin
    print,'Rebinning images to required pixel size'
    im=rebin(im,nx,ny,nimage)
    
    ;and rebinning the noise estimate. This just rebins the values
    ;and the reduction in noise from rebinning is calculated later using
    ;the sumfactor variable
    noisenow=rebin(noisenow,nx,ny,nimage)
 
    sumfactor=sumfactor*(hdr.naxis1/float(nx))*(hdr.naxis2/float(ny))
    
    ;adjust header values to account for rebinning
    shrink=float(nx)/hdr.naxis1
    hdr.naxis1=nx
    hdr.naxis2=ny
    hdr.cdelt1=hdr.cdelt1/shrink
    hdr.cdelt2=hdr.cdelt2/shrink
    hdr.crpix1=hdr.crpix1*shrink
    hdr.crpix2=hdr.crpix2*shrink
  
  endif

  print,'Calculating image geometry'
  wcs=fitshead2wcs(hdr[0])
  c=wcs_get_coord(wcs)
  ht=sqrt(total(c^2,1))
  rsun=(pb0r(hdr[0].date_obs,/arcsec,/earth))[2]
  ht=ht/rsun
  ;...and set all pixels above maxht to NAN
  indnan=where(ht gt maxht,cntnan)
  one2n,indnan,ht,ixnan,iynan
  for i=0,nimage-1 do begin
    im[ixnan,iynan,lonarr(cntnan)+i]=!values.f_nan
    noisenow[ixnan,iynan,lonarr(cntnan)+i]=!values.f_nan
  endfor
  
  if n_elements(clip) eq 4 then begin
    print,'Clipping image (if set)'
    ;clip image, noise and image coordinate array to user-defined region [xleft,ybottom,xright,ytop] in pixels
    im=im[clip[0]:clip[2],clip[1]:clip[3],*]
    noisenow=noisenow[clip[0]:clip[2],clip[1]:clip[3],*]
    if iwl eq 0 then c=c[*,clip[0]:clip[2],clip[1]:clip[3]]
  endif
  
  print,'Estimating uncertainties'
  sumfactor=sumfactor*nimage
  noisenow=noisenow/sqrt(sumfactor)
  
  ; DN uncertainty due to dark subtraction
  darknoise = 0.18
  ; DN uncertainty due to read noise
  readnoise = 1.15/sqrt(sumfactor)
  ; DN uncertainty due to quantization
  quantnoise = 0.288819/sqrt(sumfactor)
  ; DN uncertainty due to onboard compression
  compressratio = terr[indterr].compress
  compressnoise = (noisenow / compressratio) > 0.288819
  lowcounts = WHERE(im lt 25, numlow) ; Linear portion of lookup table
  if numlow gt 0 then compressnoise[lowcounts] = 0.
  compressnoise = compressnoise/sqrt(sumfactor)

  noisenow = SQRT(noisenow^2. + darknoise^2. + readnoise^2. + quantnoise^2. + compressnoise^2.)

  print,'Normalizing by exposure time'
  ;normalize by exposure time
  for i=0,nimage-1 do begin
    im[*,*,i]=im[*,*,i]/exptime[i]
    noisenow[*,*,i]=noisenow[*,*,i]/exptime[i]
  endfor
  
  
  ;combine images (if more than one)
  if nimage gt 1 then begin
    print,'Combining ',nimage,' images'
    masknan=total(long(~finite(im)),3) gt 0
    im=mean(im,dim=3,/nan)
    noisenow=mean(noisenow,dim=3,/nan)
    ind=where(masknan,cnt)
    if cnt gt 0 then begin
      im[ind]=!values.f_nan
      noisenow[ind]=!values.f_nan
    endif
  endif
 
  hdr.exptime=1
  
  ;create main output datacube at first file
  if n_elements(out) eq 0 then begin
      sz=size(im)
      nxout=sz[1] & nyout=sz[2]
      out=fltarr(nxout,nyout,nwl)
      noise_estimate=fltarr(nxout,nyout,nwl)
  endif
  
  out[*,*,iwl]=im
  noise_estimate[*,*,iwl]=noisenow
  hdrmain[iwl]=ptr_new(hdr)
  
  tai[iwl]=anytim2tai(hdr.date_obs)
  cdelt[iwl]=cdeltmain
  
  if n_elements(clip) eq 4 then begin
    crpix1[iwl]=hdr.crpix1-clip[0]
    crpix2[iwl]=hdr.crpix2-clip[1]
  endif

endfor
  
print,'Adjusting header'
sz=size(im)
hdrsave=hdr
hdrsave.naxis1=sz[1]
hdrsave.naxis2=sz[2]
hdrsave.wavelnth=0
hdrsave.wave_str='tempmap'
hdrsave.date_obs=anytim2cal(mean(tai),form=11)   
hdrsave.cdelt1=median(cdelt)
hdrsave.cdelt2=median(cdelt)
hdrsave.crpix1=median(crpix1)
hdrsave.crpix2=median(crpix2)
hdrsave.r_sun=hdrsave.r_sun*double(nx)/hdr.naxis1

end
