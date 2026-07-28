pro mcmc_plot,logt,simdem,demerr,simprb,schme,slect = slect,col_rng=col_rng,$
ps_fil=ps_fil,col_tabl=col_tabl, col_shft = col_shft, clev = clev, _extra = e 

;+ 
;procedure   mcmc_plot
;
;       This procedure plots Monte Carlo Markov-Chain DEM
;       reconstruction results. One sigma uncertanties are 
;       shown as a grey scale (or otherwise specified) with a gradient
;       determined by how the N simulated DEMs (which  ordain the one sigma
;       confidence intervals) are distributed about the best-fit DEM. 
;       Alternatively, the gradient may be determined by the p(D|M) of 
;       each DEM.
; 
;       NOTE: ALL REQUIRED INPUTS ARE STANDARD MCMC_DEM() OUTPUTS
;
;parameters
;       logt    [INPUT; required] array which determines the mid-bin values 
;               of the temperature grid on which DEM is defined
;       simdem  [INPUT; required] two dimensional array [NT, NSIM+1]
;               simulated arrays, the last of which is the best fit. 
;       demerr  [INPUT; required] confidence bounds on MAP estimates of DEM
;		* DEMERR(*,0) are lower bounds, DEMERR(*,1) are upper bounds
;       simprb  [INPUT; required] the p(D|M) of each DEM
;            
;       schme    [INPUT] the scheme with which plot colors are determined and
;                        the plotting method employed. SCHME can be:
;                         
;                   'DIFFERENCE'  - Colors are determined by
;                                   difference of each DEM simulation
;                                   bin value from the best fit DEM
;                                   values. Colors are polyfilled
;                                   between simulated values within
;                                   each  logT bin. 
;                   'INDEX'       - Colors are determined by the order in which
;                                   they appear (starting from the best
;                                   DEM outwards). Colors are polyfilled
;                                   between simulated values within each
;                                   logT bin. 
;                   'SPLINE INDEX'- Colors are determined by the order in which
;                                   they appear (starting from the best
;                                   DEM outwards). Colors are polyfilled
;                                   between simulated values within each
;                                   logT bin. Spline interpolated
;                                   curves are plotted in lieu of
;                                   histogram style bins
;                   'PROB'        - Colors are determined by the p(D|M)
;                                   of each DEM. Values are oploted
;                                   with thick bars between simulated
;                                   values within each  logT bin.
;                   'SPLINE'      - Simulated DEMS are spline
;                                   interpolated and colors determined
;                                   by the p(D|M) associated with
;                                   each. Each spline is oploted onto
;                                   a best-fit dem plot. 
;                   'NICE'        - For each temperature bin, a median simulation
;                                   value is calculated. The median values will be
;                                   plotted (dotted line) together with the DEM
;                                   simulation (solid line) that best fits them.
;                                   Uncertainty will be displayed by plotting the MIN and
;                                   MAX simulation values for each
;                                   bin. (dashed lines) Use the SLECT and CLEV
;                                   keywords to limit which simulations to include.
;                   'NICE DIFFERENCE' - A combination of DIFFERENCE AND NICE
;                                   Colors are determined by
;                                   difference of each DEM simulation
;                                   bin value from the best fit DEM
;                                   simulation to the mean value as
;                                   described in the 'NICE'
;                                   definition.  Colors are polyfilled
;                                   between simulated values within
;                                   each  logT bin. 
;                   'ERRBAR'      - Just plot error bars using DEMERR
;         
;keywords
;       slect    [INPUT] choose which DEMS to include default=[2]
;            
;                   1:ALL     -Include all simulations
;                   2:CHI^2   -Limit simulations to best 50% p(D|M) DEMS
;                              (default = 50, use CLEV keyword to toggle)
;                   3:Errors  -Limit simulations to values in each bin
;                              that lie within the specified confidence bounds (demerr)
;       clev    [INPUT] To be used in conjunction with slect=2,
;                        determines what percentage by which to limit
;                        sims in p(D|M) 
;       col_rng  [INPUT] index range in the color table to cover. 
;                        gradient will taper
;                        off from 0 (darkest) to 
;                        col_rng (brightest). (maximum is 255)
;       col_tabl [INPUT] color table on which to base the gradient       
;       col_shft [INPUT] number of color indices by which to shift
;                        brighter(or darker if black background), col_rng will be adjusted accordingly
;       ps_fil   [I/O]   name of ps file to output. If not specified 
;                        plot will be sent to current device
;
;restrictions
;       The SPLINE and SPLINE INDEX plotting schemes may not be used with slect = 3
;subroutines 
;       mid2bound
;       peasecolr
;       stample
;history
;    
;    (LL 6/03) 
;    ADDED 'Nice' and 'Nice Difference' option for schme parameter (LL 2/04)
;    ADDED call to STAMPLE (LL 2/04)
;    ADDED call to PEASECOLR to restore settings after loadct (LL 2/04)
;    ADDED recognition of weather plot background is black or white
;          and and adjustment to col_shft accordingly (LL 2/04) 
;    FIXED SPLINE INDEX schme so that it actually polyfills and doesn't just
;          oplot i.e. no more 'speckled' look. Also added a NOERASE
;          plot statement to prevent polyfills from hiding tickmarks (LL 2/04) 
;    FIXED color schemes so that the background color recognition
;          scheme would work with both psuedo color and true color (LL 2/04)
;    REMOVED PEASECOLR because IDL dynamically updates the x windows
;          display whith the current color table (LL 2/04) 
;    BUGFIX variables CUPBND and CLOBND should be of size nT not nkpt (LL 3/04)
;          if this affected you the routine would have crashed. 
;-
ok='ok' & np = n_params() & sze= size(simdem) & nT = n_elements(logt)
nERR = sze(demerr)

if keyword_set(col_tabl) then col_tabl = col_tabl else col_tabl = 0
loadct, col_tabl
if keyword_set(ps_fil) then begin
my_device='x'          ;assuming prior output is 'X'
set_plot, 'PS' 
device, filename =ps_fil, /color 
endif 
if keyword_set(col_rng)  then col_rng  = col_rng<255L  else col_rng  = 255L
if keyword_set(col_shft) then col_shft =fix(col_shft) else col_shft = 0
;if !p.background eq 0 then col_shft = -col_shft+2*255L ;reverse colors for black background
; this line doesn't work for psuedo_color only for true_color :( 
if !p.background eq 0 then alt_cols=reverse(indgen(256L)) 
if keyword_set(schme) then schme = schme else schme = 'PROB'
if keyword_set(slect) then slect = slect else slect = 2


if np lt 4 then ok = 'insufficient parameters' else $ 
 if (schme eq 'SPLINE') and (slect eq 3) then ok = 'SPLINE may not be used with slect = 3' else $
  if (schme eq 'SPLINE INDEX') and (slect eq 3) then ok = 'SPLINE INDEX may not be used with slect = 3' else $
   if (slect gt 3) or (slect le 0) then ok = 'slect keyword selection not valid' else $
    if (col_tabl lt 0 ) then ok =  'col_tabl keyword selection not valid' 
     
if ok ne 'ok' then begin 
   print, 'Usage:  mcmc_plot,logt,simdem,demerr,simprb,schme=schme,slect = slect,$'
   print, 'col_rng=col_rng,ps_fil=ps_fil,col_tabl=col_tabl, col_shft = col_shft'
   if np ne 0 then message, ok, /info 
   return
endif 


TBB = mid2bound(logt)  ;LOGT grid bin boundaries

 nsim = n_elements(simdem(0,*))-1
 demerrl = demerr[*,0]
 demerru = demerr[*,1]

;initialize container arrays for SPLINE INDEX 
if schme eq 'SPLINE INDEX' then begin 
 splndxu=[simdem-simdem]
 splndxl=[simdem-simdem]
endif

;initialize container arrays for NICE 
if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then begin 
 meanofeachbin = fltarr(nT) 
 maxofeachbin  = fltarr(nT) 
 minofeachbin  = fltarr(nT) 
endif 


;if not opting for SPLINE interpolation or 'ERRBAR' then  do this
if (schme ne 'SPLINE') and (schme ne 'ERRBAR') then begin 

 ;initial DEM plot to lay basis for polyfill bars
 ; _extra keywords are caught here
 if schme ne 'NICE' and schme ne 'NICE DIFFERENCE' then begin
     plot,logt,simdem(*,nsim),psym=10, /ylog,/ynoz,xtitle='logT'$
       ,ytitle = 'DEM', title = 'Best DEM + error bands  ', thick = 5,$
       yrange = [0.8*min(simdem), 1.2*max(simdem)],$
       xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1, /nodata,_extra = e
 endif

;loop through each bin in logT space 
for dd = 0, nT-1 do begin 
 ;initialize some variables for polyfill/bar loop according to selection mechanism
     if slect eq 1 then begin 
         upbnd     = max(simdem(dd,*))
         lobnd     = min(simdem(dd,*)) 
         sortbyem  = sort(simdem(dd,*))
         sortdem   = simdem(dd,sortbyem)
         sortprb   = simprb(sortbyem)
     endif
     if slect eq 3 then begin 
         upbnd     = demerru(dd)
         lobnd     = demerrl(dd)
         sortbyem  = sort(simdem(dd,*))
         sortdem   = simdem(dd,sortbyem)
         sortprb   = simprb(sortbyem)
     endif
     if slect eq 2 then begin 
         if keyword_set(clev) then clev = clev else clev = 50
         cutnumber = fix(((clev/100d)<1.0)*nsim)
         cutarray  = simprb(sort(simprb))
         cutvalu   = cutarray(cutnumber)         
         sortndx   = where(simprb le cutvalu)
         tmpo      = simdem(dd,sortndx)      
         tmpoprb   = simprb(sortndx)
         sortbyem  = sort(tmpo) 
         sortdem   = tmpo(sortbyem)
         sortprb   = tmpoprb(sortbyem)
         upbnd     = max(tmpo(sortbyem))   
         lobnd     = min(tmpo(sortbyem)) 
     endif
 ;store initialized arrays   
     ;on first iteration initialize container arrays 
     if dd eq 0 then begin 
         nkpt = n_elements(sortprb) & cupbnd= dblarr(nT) & clobnd= dblarr(nT) 
         csortdem = dblarr(nT,nkpt) & csortprb=csortdem 
     endif 
     ;populate the container arrays 
     cupbnd(dd) = upbnd & clobnd(dd) = lobnd 
     csortdem(dd,*) = sortdem & csortprb(dd,*) = sortprb 
     ;if NICE or NICE DIFFERENCE set 
     if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then begin 
         meanofeachbin(dd)=median(sortdem)
         maxofeachbin(dd)=upbnd 
         minofeachbin(dd)=lobnd  
     endif ;the NICE/NICE DIFFERENCE setm if 
endfor ; the loop through each bin in logT space loop

;do the nice and nice difference stuff difference stuff
if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then begin     
  if slect eq 2 then ll=sortndx else ll=findgen(nsim) & cutsamp = simdem(*,ll) 
  ;find best fit of SIMDEM to MEANOFEACHBIN.
  ;use a statistic that takes the squared residuals weighted by (UPBOUND-LOBND)
    mystat = fltarr(n_elements(ll))
  for j = 0, n_elements(ll)-1 do mystat(j)= $
          total(((meanofeachbin-simdem[*,ll[j]])^2)/(maxofeachbin-minofeachbin))
    jnk = min(mystat,kk) 
  if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then begin 
       plot, logt, simdem[*,kk],linestyle=0, /ylog,/ynoz,xtitle='logT',$
       ytitle = 'DEM', title = 'Best DEM + error bands  ', thick = 5,$
       yrange = [0.8*min(simdem), 1.2*max(simdem)],$
       xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1 , _extra = e;, psym = 10
    if schme eq 'NICE' then begin 
        oplot, logt, meanofeachbin, linestyle = 1, thick = 3 ;, psym = 10
        oplot, logt, maxofeachbin, linestyle = 2, thick = 3 ;, psym = 10
        oplot, logt, minofeachbin, linestyle = 2, thick = 3 ; , psym = 10
    endif
  endif
endif ; the NICE/NICE DIFFERENCE if

;loop through each bin in logT space (again) and do the pollyfill et.al)
for dd = 0, nT-1 do begin 
 upbnd = cupbnd(dd) & lobnd = clobnd(dd) 
 sortdem = transpose(csortdem(dd,*)) & sortprb = transpose(csortprb(dd,*)) 
 ;prepare the simulated DEM arrays to be passed through the polyfill mechanism
 ;if schme 1 or to and then loop through sim values to polyfill
  if (schme eq 'DIFFERENCE') or (schme eq 'INDEX') or (schme eq 'NICE DIFFERENCE') then begin
         bestndx = where(sortdem eq simdem[dd, nsim]) ;  
         if schme eq 'NICE DIFFERENCE' then bestndx = where(sortdem eq simdem[dd,kk])
         jnk1=min(abs(1-upbnd/sortdem(*)), undx)      ;
         jnk2=min(abs(1-lobnd/sortdem(*)), lndx)      ; 
         erorintu =sortdem(undx)-sortdem(bestndx(0))  ;
         erorintl =sortdem(bestndx(0))-sortdem(lndx)  ;
 ;first loop for sim values below best value
     for qq = lndx, bestndx(0)-1 do begin 

     ;determine the color for each sim according to the schme selected
         if schme eq 'DIFFERENCE' or schme eq 'NICE DIFFERENCE' then clr = fix( (col_rng-col_shft)*$
                                             ((sortdem(bestndx(0))-sortdem(qq)))/erorintl ) + col_shft
         if schme eq 'INDEX'      then clr = fix( (col_rng-col_shft)*$
                                             ((bestndx(0) - qq))/(bestndx(0)-lndx) ) + col_shft
         if !p.background eq 0 then clr = alt_cols(clr) 
         polyfill, [TBB(dd),TBB(dd+1),TBB(dd+1),TBB(dd)] $
           , [sortdem(qq),sortdem(qq), sortdem(qq+1), sortdem(qq+1)] $
           , color = clr(0),noclip = 0
     endfor; end loop through sim values -below

   ;then loop for sim values above best value
     for qq = bestndx(0), undx-1 do begin 
         if schme eq 'DIFFERENCE' or schme eq 'NICE DIFFERENCE' then clr = fix( (col_rng-col_shft)*$
                                             ((sortdem(qq)-sortdem(bestndx)))/erorintu ) + col_shft 
         if schme eq 'INDEX'      then clr = fix( (col_rng-col_shft)*$
                                             ((qq-bestndx(0)))/(undx-bestndx(0)) ) + col_shft 
         if !p.background eq 0 then clr = alt_cols(clr) 
         TBB = mid2bound(logt)  ; LOGT grid bin boundaries
         polyfill, [TBB(dd),TBB(dd+1),TBB(dd+1),TBB(dd)] $
           , [sortdem(qq),sortdem(qq), sortdem(qq+1), sortdem(qq+1)] $
           , color = clr(0), noclip = 0
     endfor ; end loop though sim values -above 

  endif; the schme eq 'DIFFERENCE' or 'INDEX' if 

  if schme eq 'PROB' then begin
      TBB= mid2bound(logt) 
      ;must plot lower confidence (high prb number) values 
      ;first so we sort again in prb space
      jnk1=min(abs(1-upbnd/sortdem(*)), undx)      
      jnk2=min(abs(1-lobnd/sortdem(*)), lndx)       
      sortdem = sortdem[lndx: undx]
      sortprb = sortprb[lndx: undx]
      sortndx = sort(sortprb)
      sortprb = sortprb(sortndx)
      sortdem = sortdem(sortndx)
      icols=fix((col_rng-col_shft)*(sortprb-min(sortprb))/(max(sortprb)-min(sortprb)))+col_shft
      if !p.background eq 0 then icols = icols+1
      os=reverse(sort(icols))
      for i=0L,n_elements(sortndx)-1L do oplot,[TBB(dd),TBB(dd+1)]$
      ,[replicate(sortdem(os(i)),2)],col=icols(os(i)) , _extra=e
  endif; the schme 'PROB' if 

  if schme eq 'SPLINE INDEX' then begin 
      bestndx = where(sortdem eq simdem[dd, nsim]) 
      jnk1=min(abs(1-upbnd/sortdem(*)), undx)      
      jnk2=min(abs(1-lobnd/sortdem(*)), lndx)       
      sortdem = sortdem[lndx: undx]
      splndxu[dd,findgen(undx-bestndx(0)+1)] = sortdem[bestndx(0):undx]
      splndxl[dd,findgen(bestndx(0)-lndx+1)] = sortdem[lndx:bestndx(0)]
  endif ;the SPLINY INDEX if 

endfor ;endfor loop through LOGT BINS

;oplot the 'best' DEM again over the painted stuff
if schme eq 'NICE' or schme eq 'NICE DIFFERENCE' then $
    oplot, logt, simdem[*,kk],linestyle=0, thick = 5 else begin 
    if schme ne 'SPLINE INDEX' then begin 
        oplot,logt,simdem(*,nsim),psym=10,thick=5
    endif
endelse

endif ;endif the not spliny if

if schme eq 'SPLINE INDEX' then begin 

    TBB= mid2bound(logt) & Tn = [logt, TBB] & sT = sort(Tn) & Tn=Tn(sT)

    ;get the  array w/ largest number of non zero elements (in DEM space)
    tst_gt0u=fltarr(nt)
    tst_gt0l=fltarr(nt)
    for qq=0, nt-1 do begin ;loop through temperature bins
        tst_gt0u(qq) = n_elements(where(splndxu[qq,*] ne 0))
        tst_gt0l(qq) = n_elements(where(splndxl[qq,*] ne 0))
    endfor
    lrgstu=max(tst_gt0u) & lrgstl=max(tst_gt0l)

    ;create array of size [nt,nlargest] and interpolate 
    ;the smaller arrays to an nlargest component array   
    nusplndxu=fltarr(nt,lrgstu) & nusplndxl=fltarr(nt,lrgstl)
    for ww=0, nt-1 do begin 
        nusplndxu[ww,*]=interpol(splndxu[ww,indgen(tst_gt0u(ww))],lrgstu)
        nusplndxl[ww,*]=interpol(splndxl[ww,indgen(tst_gt0l(ww))],lrgstl)
    endfor

    ;SPLINE interpolate to grid NT and polyfill plot
    splndu = 10^( spline(logt, alog10(nusplndxu[*,0]),Tn) )
    for ee=1,lrgstu-1 do begin 
        colru = fix( (col_rng-col_shft)*(ee/lrgstu) ) + col_shft
        if !p.background eq 0 then colru = alt_cols(colru)
        splndudn1=splndu
        splndu=10^( spline(logt, alog10(nusplndxu[*,ee]),Tn) )
          ;need this because Spline interpolation fuzzy:
          tmp = [[splndudn1],[splndu]]
          for j = 0,n_elements(Tn)-1 do splndu(j)=max(tmp(j,*)) 
          for j = 0,n_elements(Tn)-1 do splndudn1(j)=min(tmp(j,*))            
          polyfill,[Tn, reverse(Tn)],[splndudn1,reverse(splndu)],  color = colru,noclip=0
    ; oplot, Tn, splndu, color = colru
    endfor
    splndl = 10^( spline(logt, alog10(nusplndxl[*,0]),Tn) )
    for rr=1,lrgstl-1 do begin
        colrl = fix( (col_rng-col_shft)*((lrgstl-rr)/lrgstl) ) + col_shft
        if !p.background eq 0 then colrl =  alt_cols(colrl)
        splndldn1=splndl
        splndl=10^( spline(logt, alog10(nusplndxl[*,rr]),Tn) )
          ;need this because Spline interpolation fuzzy:
          tmp = [[splndldn1],[splndl]]
          for j = 0,n_elements(Tn)-1 do splndl(j)=max(tmp(j,*)) 
          for j = 0,n_elements(Tn)-1 do splndldn1(j)=min(tmp(j,*))      
        polyfill,[Tn, reverse(Tn)],[splndldn1,reverse(splndl)], color = colrl,noclip=0
;    oplot, Tn, splndl, color = colrl
    endfor

endif; the SPLNIY INDEX if

if schme eq 'SPLINE' then begin 

    if slect eq 1 then sortndx = findgen(nsim) 
    if slect eq 2 then begin 
        if keyword_set(clev) then clev = clev else clev = 50
        cutnumber= fix(((clev/100d)<1.0)*nsim)-1
        cutarray = simprb(sort(simprb))
        cutvalu = cutarray(cutnumber)         
        sortndx = where(simprb le cutvalu)
    endif 

    TBB= mid2bound(logt) & Tn = [logt, TBB] & sT = sort(Tn) & Tn=Tn(sT)
    sortprb=simprb(sortndx)
    sortdem=simdem[*,sortndx]
    icols=fix((col_rng-col_shft)*(sortprb-min(sortprb))/(max(sortprb)-min(sortprb)))+col_shft
    if !p.background eq 0 then icols = alt_cols(icols)
    os=reverse(sort(icols))

    ;initial DEM plot to lay basis for bars
    ; _extra keywords are caught here
    plot,Tn,10^( SPLINE(logt,alog10(simdem(*,nsim)),Tn )), /ylog, /ynoz,xtitle='logT'$
      ,ytitle='DEM', title = 'Best DEM + error bands  ', thick = 5,$
      yrange = [0.8*min(simdem), 1.2*max(simdem(*))],$
      xrange = [min(logt), max(logt)], xstyle = 1, _extra = e

    for i = 0L, n_elements(sortndx)-1L do $
      oplot, Tn,10^( SPLINE(logt,alog10(sortdem(*,os(i))),Tn) ),col=icols(os(i)), thick = 2

endif ;the SPLINY if 

if schme eq 'ERRBAR' then begin 

    plot,logt,simdem(*,nsim),psym=10, /ylog,/ynoz,xtitle='logT'$
      ,ytitle='DEM', title = 'Best DEM + error bands  ', thick = 5,$
      yrange = [0.8*min(simdem), 1.2*max(simdem)],$
      xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1 , _extra = e
    for j = 0,n_elements(logt)-1 do oplot, [logt(j),logt(j)], [demerrl(j), demerru(j)] 

endif ; the ERRBAR if

;do some /noerase /nodata plots to ensure tick marks arn't painted over
     plot,logt,simdem(*,nsim),/yl,/ynoz,$
        yrange = [0.8*min(simdem), 1.2*max(simdem)],$
        xrange = [min(TBB), max(TBB)], xstyle = 1, ystyle = 1, /nodata,/noerase,_extra=e

stample,_extra=e;put PoA stamp on corner of plot
;peasecolr,_extra=e;restore default PoA colors.. no don't do
;this.. because the x window plot is dynamically updated to these
;colors and mahem ensues
if keyword_set(ps_fil) then begin 
    device, /close_file
    set_plot, my_device
endif 

end 
