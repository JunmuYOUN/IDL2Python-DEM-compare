pro dem_example_sites

synoptic=1b;work with synoptic AIA images
startdate='2011/01/01 00:00'
enddate='2011/01/01 00:30'
npixin=512;synoptic images are 1024x1024, we will rebin to 512x512
clip=[220,200,400,400];Define small region near image center, note in rebinned pixel values
wl=[171,193,211,131,94,335];define required channels

nwl=n_elements(wl)

;define top path level of where you keep AIA synoptic files.
;below this level, directory structure should be yyyy/mm/dd/www where www gives the channel wavelength e.g. 171 or 304 
topdir=getenv('AIA_SYNOPTIC')
dirs=topdir+'/'+anytim2cal(startdate,/da,form=11)

files=ptrarr(nwl);set up pointer array to store filenames
for iwl=0,nwl-1 do begin
  wlstr=strcompress(string(wl[iwl]),/remove_all)
  if strlen(wlstr) lt 3 then wlstr='0'+wlstr;so '94' goes to '094'
  f=file_search(dirs+'/'+wlstr+'/AIA*.fits')
  ;you can write in chunk here to limit files to certain time range (e.g. using startdate and enddate above)
  ;for convenience, I just use first 3 files for each channel.
  f=f[0:2]
  print,f
  print,''
  files[iwl]=ptr_new(f)
endfor

;open the datacube and associated errors
dem_openaiafiles_sites,files,datacube,hdrmain,header,errors,npixuser=npixin, $
                clip=clip,synoptic=synoptic

;now set up values for DEM inversion
ntemp=43
tempra=[4.e5,1.e7]
sz=size(datacube)
nx=sz[1] & ny=sz[2] & nwl=sz[3]
date=header.date_obs

;read in response file
s=aia_get_response(/temp,/dn,/evenorm,timedepend_date=date,/chiantifix,phot=phot)

;and set up arrays ready for inversion, based on response structure
dem_aiainterpolresponse_sites,wl,ntemp,alog10(tempra),s,logte,all,allerr,dtemp=dtemp,/reglog


convergence=0.04;4% tolerance, see Sol Phys paper

print,'Using SITES'
time0=systime(1)

dem=dblarr(nx,ny,ntemp);output DEM array
demerr=dblarr(nx,ny,ntemp);output DEM error array
obsmod=fltarr(nx,ny,nwl);output model measurements so you can check residuals

;loop through all pixels
for ix=0l,nx-1 do begin
  if ix mod ceil(nx/10.) eq 0 then print,long(ix*100/float(nx)),'%'
  
  for iy=0l,ny-1 do begin
    dem[ix,iy,*]=dem_sites(reform(datacube[ix,iy,*]),reform(errors[ix,iy,*]),all,allerr,dtemp, $
                    conv=convergence,demerr=demerrnow,obsmod=obsmodnow,ker=ker,res2=res2,totres2=totres2)
    demerr[ix,iy,*]=demerrnow
    obsmod[ix,iy,*]=obsmodnow
  endfor

endfor
tproc=systime(1)-time0
print,'Time taken ',tproc
print,'Processing rate ',nx*ny/tproc,' pix per s'

stop
;Add your own bits here to save DEM output!

end

