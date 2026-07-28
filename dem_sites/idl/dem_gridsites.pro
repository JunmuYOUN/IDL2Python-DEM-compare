;+
;
; NAME:
; dem_gridsites
;
; PURPOSE:
; Given (i) a set of multi-pixel input measurements, (ii) associated measurement uncertainties,
; (iii) temperature response curves, (iv) response relative uncertainties and
; (v) the temperature bin widths, this code implements the Grid-SITES wrapper for inversion described in
; https://ui.adsabs.harvard.edu/abs/2019SoPh..294..136P/abstract and returns the estimated 
; differential emission measure (DEM) for each pixel. This code uses the SITES method for inversion, described in 
; https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract.  

;
; CALLING SEQUENCE:
; dem=dem_gridsites(obs_in,err_in,response,response_err,delta_temp,wavel, $
;               [convergence=convergence, obsmodmain=obsmodmain,demerrmain=demerrmain, $
;               goodoffit_main=goodoffit_main,grid=grid,nprocess=nprocess])
;
; INPUTS:
;   obs_in = multi-pixel, multi-channel measurement. E.g. for a thousand-pixel, 6-channel AIA/SDO measurement, 
;             this would be a [1000,6] array in DN/s units. The order of the channels is not important, 
;             as long as the same ordering is used for all relevant inputs (i.e. response functions and errors).
;   err_in = the errors or uncertainties associated with the measurements,
;               so array same size (and in same channel order) as obs_in.
;               These are absolute errors, not relative. Same units as obs_in (e.g. DN/s)
;   response = the temperature response function for each channel, e.g.  as provided by the AIA/SDO Solarsoft routines.
;               E.g. for the 1000 pixel, 6-channel AIA/SDO example, this would be a [N-temperature,N-channel] array e.g. [43,6]
;   response_err = the RELATIVE uncertainties of the temperature response functions, one uncertainty for each channel.
;               For the 6-channel AIA example, response_err will be a 6-element vector. Fractional
;               uncertainties (not percentage!).
;               I use, for channels [94,131,171,193,211,304,335], uncertainties of [0.5,0.5,0.25,0.25,0.25,0.5,0.25]
;   delta_temp = the width of each temperature bin, in Kelvin. From trial and error, this code works best on
;               fixed logarithmic temperature bin width, so that delta_temp in Kelvin increases with temperature
;               (bins are wider with increasing temperature). This is a vector with temperature
;               bins corresponding to the input response functions, e.g [43 bins]
;   wavel = the wavelengths of each channel e.g. [94,131,171,193,211,335]
;
; OUTPUTS:
;   This function returns the DEM at each pixel and temperature bin
;
; OPTIONAL KEYWORD
;   convergence = sets the convergence threshold for the core SITES inversion, default 1.e-2 (1%)
;
; OPTIONAL INPUTS
;   grid = user can provide a structure containing the Grid-SITES grid and DEM results from a previous call
;           to Grid-SITES in order to speed up inversions considerably. See paper, particularly section 6 on
;           application to time series. If this structure is provided, Grid-SITES will only process measurement
;           bins that have not been previously recorded in the grid structure. The grid structure will be 
;           updated with the newly-processed measurements, enabling improved efficiency with each subsequent call. 
;           
;
; OPTIONAL OUTPUTS
;   obsmodmain = returns the model measurements, so user can compare the input and model measurements
;   demerrmain = gives the (absolute) errors for the output DEM
;   goodoffit_main = the goodness of fit of measurement to model at each pixel (averaged over all channels)
;   grid = outputs a structure containing the Grid-SITES grid and DEM results. Subsequent calls to Grid-SITES
;           can provide this grid to speed up inversions considerably. See paper, particularly section 6 on
;           application to time series.  
;
; PROCEDURE:
; See https://ui.adsabs.harvard.edu/abs/2019SoPh..294..136P/abstract
; 
; EXAMPLE:
; For time series of Nt images at Nwl channels, each of size 1024x1024 pixels here is an pseudocode/code loop. Assume 
; you have temperature response curves and associated errors. Number of temperature bins is Ntemp
; 
; wavel=[94,131,171,193,211,304,335]
; response=response curves as extracted from e.g. aia_get_response.pro
; response_errors = [0.5,0.5,0.25,0.25,0.25,0.5,0.25]
; delta_temp=width of each temperature bin
; for i=0,Nt-1 do begin
; 
;   current_image=read in current file and errors using read_sdo. Image is [nx,ny,nwl]
;   current_errors = appropriate user-written function for calculating errors. Error is [nx,ny,nwl]
;   
;   ;CALL GRID-SITES
;   current_dem=dem_gridsites(current_image,current_errors,response,response_errors,delta_temp,wavel,grid=grid)
;   ;AT FIRST ITERATION, GRID IS RETURNED AS OUTPUT. SUBSEQUENT CALLS TO DEM_GRIDSITES USES
;   ;GRID TO IMPROVE EFFICIENCY.
;   
;   ;SAVE CURRENT DEM
;   save,current_den, filename='appropriate path and filename.dem'
; 
; endfor
;
; USE & PERMISSIONS
; If you use this code for DEM inversions, please cite both SITES and Grid-SITES papers:
;   https://ui.adsabs.harvard.edu/abs/2019SoPh..294..136P/abstract
;   https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract
; Any problems/queries, or suggestions for improvements, please email Huw Morgan, hmorgan@aber.ac.uk
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC Studentship to Aberystwyth University (Pickering)
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 07/2019 - James Pickering & Huw Morgan hmorgan@aber.ac.uk
;
;
;-


function dem_gridsites,obs_in,err_in,response,response_err,logt,delta_temp,wavel, $
                          convergence=convergence,demerrmain=demerrmain,obsmodmain=obsmodmain, $
                          goodoffit_main=goodoffit_main,grid=grid,nprocess=nprocess

;query size of arrays
sz=size(obs_in)
np=sz[1]
nwl=sz[2]
sz=size(response)
nt=sz[1]

;work in log base 10 of intensities
m=alog10(obs_in)

if n_elements(grid) eq 0 then begin
  print,'Calculating GRID-SITES parameters'
  ;find minimum,maximum ranges in each channel
  ;also find median value in each channel (not log). This is used to set number of bins for each channel
  mn=dblarr(nwl)
  mx=dblarr(nwl)
  md=dblarr(nwl)
  for iwl=0,nwl-1 do begin
    mn[iwl]=robust_min(m[*,iwl],0.02,max=mxnow);min(m[*,iwl],max=mxnow);
    mx[iwl]=mxnow
    md[iwl]=median(m[*,iwl])
  endfor
  
 
  ;here we calculate the number of bins for each channel, based on  
  ;interpolation of median intensities. So highest median intensity channel
  ;will have mxbin bins, and lowest median intensity channel gets mnbin bins.
  mxbin=40;40
  mnbin=30;30
  nbin=round(mnbin+(md-min(md))*(mxbin-mnbin)/(max(md)-min(md)),/l64)
  
  dbin=(mx-mn)/nbin;width of bins for each channel

  
endif else begin
  
  print,'Reading GRID-SITES parameters from existing grid'
  nbin=grid.nbin
  mn=grid.mn
  mx=grid.mx
  mxbin=max(nbin)
  mnbin=min(nbin)
  dbin=(mx-mn)/nbin;width of bins for each channel
  
endelse


print,'Calculating pixel grid indices'
;calculate indices, as if creating histogram
;imagine a nwl-dimensional hypercube or histogram
;certain combinations of intensities will populate different bins in this histogram
;I avoid creating histogram (since inefficient), just need to record which histogram bins are populated
;following chunk achieves this very quickly, without strain on memory
;also calculate mean measurement intensity and error for each bin (and standard deviations)
index=ulon64arr(np)
mmn=fltarr(mxbin,nwl);mean measurement
emn=fltarr(mxbin,nwl);mean error
mst=fltarr(mxbin,nwl);std dev measurement
est=fltarr(mxbin,nwl);std dev error
indok=bytarr(np)+1b;assume all pixels within grid
for iwl=nwl-1,0,-1 do begin
  ind=floor((m[*,iwl]-mn[iwl])*nbin[iwl]/(mx[iwl]-mn[iwl]),/l64);maps input intensities into grid bins for this channel
  
  ;pixels which don't fit within bin range (due to using robust min and max in calculating ranges for each channel)
  ;if any channel has values outside the bin range, that bin is discarded using indok (eq 0)
  indok=indok and (ind ge 0 and ind le (nbin[iwl]-1))
  
  for ibin=0,nbin[iwl]-1 do begin
    ;which data points contribute to each bin?
    ind2=where(ind eq ibin,cnt2)
    if cnt2 eq 0 then continue
    emn[ibin,iwl]=mean(err_in[ind2,iwl])
    mmn[ibin,iwl]=mean(obs_in[ind2,iwl])
    est[ibin,iwl]=stddev(err_in[ind2,iwl])
    mst[ibin,iwl]=stddev(obs_in[ind2,iwl])
  endfor
  index=((index*nbin[iwl])+ind)*indok;main index. Note product with indok to prevent runaway indices for pixels outside grid
endfor

;identify which pixels are within grid, and small number (hopefully, in principle) outside 
indin=where(indok,cntin)
indout=where(~indok,nout)
print,nout,' pixels outside of grid (will be processed individiually)' 

;now find uniq combinations of channel intensities included in grid
if n_elements(grid) gt 0 then begin
  
  indexuniq=(index[indin])[uniq(index[indin],sort(index[indin]))]
  cnt=n_elements(indexuniq)
  
  isnew=-1
  cntnew=0
  for i=0,cnt-1 do begin
    ind2=where(grid.index eq indexuniq[i],cnt2)
    if cnt2 gt 0 then continue
    isnew=[isnew,indexuniq[i]]
    cntnew++
  endfor
  if cntnew gt 0 then isnew=isnew[1:*]
  indexuniq=isnew
  
  ;hgrid=histogram(grid.index,min=0,max=nindex-1,/l64)
  ;isnew=~(hgrid gt 0) and (h gt 0)
  ;indexuniq=xh[where(isnew,cntnew)]
  
endif else begin
  
  ;unique grid bin indices
  indexuniq=(index[indin])[uniq(index[indin],sort(index[indin]))]
  cntnew=n_elements(indexuniq)

endelse

if cntnew gt 0 then begin
  
  ndem=n_elements(indexuniq)
  
  ;and convert index into individual bin indices for each channel
  one2n,indexuniq,nbin,i0,i1,i2,i3,i4,i5,i6,/dim
  
  ;now need to calculate DEMs for each populated bin (elements of index)
  
  griddem=dblarr(ndem,nt)
  griddemerr=dblarr(ndem,nt)
  gridobsmod=dblarr(ndem,nwl)
  gridgof=dblarr(ndem)
  gridcntdata=lonarr(ndem)
  ;loop through each populated grid bin, creating DEMs and errors etc for each
  print,'Calculating ',ndem,' DEMs for grid points'
  
  for idem=0,ndem-1 do begin
  
      if idem mod 2000 eq 0 then print,idem,' out of ',ndem-1 
  
      ;for this populated bin, extract bin intensity values for each channel (and error etc)
      ;indgrid=[i0[idem],i1[idem],i2[idem],i3[idem],i4[idem],i5[idem],i6[idem]];make index
      
      ;following chunk a bit clumsy, but accounts for different number of channels
      indgrid=i0[idem]
      for iwl=1,nwl-1 do begin
        com='indgrid=[indgrid,i'+strcompress(string(iwl),/remove_all)+'[idem]]'
        void=execute(com)
      endfor
      mnow=mmn[indgrid,indgen(nwl)];channel intensities
      snow=mst[indgrid,indgen(nwl)];standard deviation of channel intensities
      enow=emn[indgrid,indgen(nwl)];channel errors
  
      ;calculate DEM for this bin
      dem=dem_sites(mnow,enow,response,response_err,delta_temp, $
                      convergence=convergence,demerr=demerr,obsmod=obsmod, $
                      ker=ker,res2=res2,totres2=totres2,gres=gres)
      
      ;the standard deviation of intensities in this bin adds additional error to the DEM
      ;do this calculation here
      gerr=snow/mnow
      demgerr=sqrt(total(rebin(reform(gerr^2,1,nwl),nt,nwl)*gres,2))   
      demerr=sqrt(demgerr^2+demerr^2)*dem
      
      ;store result for grid structure                
      griddem[idem,*]=dem
      griddemerr[idem,*]=demerr
      gridobsmod[idem,*]=obsmod
      gridgof[idem]=sqrt(total(((mnow-obsmod)^2)/enow^2))/nwl
  
  endfor
endif

if n_elements(grid) eq 0 then begin
  ;create structure to record DEM grid. 
  grid={nbin:nbin,logt:logt,wl:wavel,mn:mn,mx:mx, $
    index:indexuniq, $
    dem:griddem, $
    demerr:griddemerr, $
    obsmod:gridobsmod, $
    gof:gridgof, $
    cntdata:lonarr(ndem)};gof is goodness of fit
endif else begin
  if cntnew gt 0 then begin
    ;update existing grid with new gridded values
    grid=replace_structure_field(grid,'index',[grid.index,indexuniq])
    grid=replace_structure_field(grid,'dem',[grid.dem,griddem])
    grid=replace_structure_field(grid,'demerr',[grid.demerr,griddemerr])
    grid=replace_structure_field(grid,'obsmod',[grid.obsmod,gridobsmod])
    grid=replace_structure_field(grid,'gof',[grid.gof,gridgof])
    grid=replace_structure_field(grid,'cntdata',[grid.cntdata,lonarr(ndem)])
  endif
endelse

;map grid-calculated results to each pixel of input data
print,'Mapping GRID-SITES DEMs into output result'
st=systime(1)
demind=-1l
dataind=-1l
for i=0,n_elements(grid.index)-1 do begin
  ind=where(index eq grid.index[i],cnt)
  if cnt eq 0 then continue
  dataind=[dataind,ind]
  demind=[demind,replicate(i,cnt)]
endfor

dataind=dataind[1:*]
demind=demind[1:*]
demmain=dblarr(np,nt)
demerrmain=dblarr(np,nt)
obsmodmain=dblarr(np,nwl)
goodoffit_main=dblarr(np)
demmain[dataind,*]=grid.dem[demind,*]
demerrmain[dataind,*]=grid.demerr[demind,*]
obsmodmain[dataind,*]=grid.obsmod[demind,*]
goodoffit_main[dataind,*]=sqrt(total(((obs_in[dataind,*]-grid.obsmod[demind,*])^2)/(err_in[dataind,*]^2),2))/nwl
print,'Mapping time = ',systime(1)-st

;then loop through pixels which don't fit in grid
if nout gt 0 then print,'Calculting DEM for ',nout,' non-gridded pixels'
for idem=0,nout-1 do begin

    if idem mod 2000 eq 0 then print,idem,' out of ',nout-1 
    inddata=indout[idem]
    mnow=obs_in[inddata,*]
    enow=err_in[inddata,*]
    dem=dem_sites(mnow,enow,response,response_err,delta_temp, $
                    convergence=convergence,demerr=demerr,obsmod=obsmod, $
                    ker=ker,res2=res2,totres2=totres2)
    demmain[inddata,*]=dem
    demerrmain[inddata,*]=demerr
    obsmodmain[inddata,*]=obsmod
    goodoffit_main[inddata,*]=sqrt(total(((obs_in[inddata,*]-obsmod)^2)/(err_in[inddata,*]^2),2))/nwl
    
endfor


nprocess=nout+cntnew

return,demmain

end