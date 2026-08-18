;***************************************************************************
;***************************************************************************
;
; README for <xrt_dem_iterative2.pro>
; 2008-Oct-05
;
; This is a complicated program with many subroutines that are mostly 
; contained within this file. The end-user should proceed directly from 
; this README section to the final program in this file, which is the one 
; named <xrt_dem_iterative.pro>. That program has an extensive header to
; explain its usage. It is not necessary for end-users to understand the
; the details of the various subroutines in order to run 
; <xrt_dem_iterative.pro>. (To find the program quickly, search for the 
; following word: wrapperstartshere .)
;
; DEM estimation is a tricky business. This solver represents one
; algorithmic approach, and the default XRT temperature responses must
; necessarily make some assumptions (e.g., about the atomic physics, etc.).
; This file cannot possibly be a primer on DEM estimation, and the end-user
; is expected to know something about the topic. However, an effort has
; been made to keep the interface as straightforward and useful as 
; possible. Furthermore, the temperature response codes 
; (<calc_xrt_temp_resp.pro> and associated routines) are designed to be
; modular, allowing end-users to "roll their own".
;
;
; This file contains the following routines:
;
; XRT_ITER_DEMSTAT
; MPDEMFUNCT
; XRT_DEM_ITER_SOLVER   --> MPFIT (external routine),
;                           MPDEMFUNCT,
;                           XRT_ITER_DEMSTAT
; MP_PREP
; COMPACT
; XRT_DEM_ITER_ESTIM    --> COMPACT
; XRT_DEM_ITER_NOWIDGET --> XRT_DEM_ITER_ESTIM, 
;                           MP_PREP,
;                           XRT_DEM_ITER_SOLVER
; XRT_DEM_ITERATIVE2    --> XRT_DEM_ITER_NOWIDGET
;
; The parent routine is XRT_DEM_ITERATIVE2. See that routine's
; header for a description of how to use this software.
;
;***************************************************************************
;***************************************************************************

pro xrt_iter_demstat,line_data,dem,t,weights=weights,abunds=abunds, $
  i_mod=i_mod, di=di,text=text,chisq=chisq

;==========================================================================
;+
; NAME: xrt_iter_demstat
;
; PURPOSE:  Return statistics given a set of contribution functions
; and a DEM curve.  Calculates modeled intensity as
;     Integral(EMIS(T)*DEM(T)*DT)
;
; CATEGORY:  DEM
;
; INPUTS:  line_data - array of structures with minimum fields of:
;                            emis   emissivity for the line
;                            t   temperature points for each bin
;                            nt     number of temperature points for each line
;                            i_obs   observed intensity of line
;                            i_err   error in observed intensity
;                            elem   element name for line
;                            ion    which ion
;                            wave   wavelength of line peak
;           dem - the modeled dem curve
;           t - log t for each dem bin
;
; OPTIONAL INPUTS:
;           weights - vector, 1 for line being used, 0 to ignore
;           abunds - abundance factor, defaults to 1.0
;
; OPTIONAL OUTPUTS:
;           i_mod - modeled intensity
;           di - difference between modelled and observed intensity
;           chisq - total chi-squared
;
; EXAMPLE:  IDL> help,line_data
;           LINE_DATA           STRUCT    = -> <Anonymous> Array[9]
;           IDL> help,line_data,/st
;           ** Structure <5c5e68>, 10 tags, length=1688, refs=1:
;               T               FLOAT     Array[200]
;               EMIS            FLOAT     Array[200]
;               NT              INT             26
;               I_OBS           FLOAT       1.90185e-07
;               I_ERR           FLOAT       4.39380e-08
;               I_UNITS         STRING    'ergs/sec/pix'
;               ELEM            STRING    'Thin Al Mesh'
;               ION             STRING    ''
;               WAVE            FLOAT           0.00000
;               IPF             STRING '/data/sdata2/XRT/spectra/apec_09_c.genx'
;           IDL> xrt_iter_demstat,line_data,dem,t,text=text
;           IDL> print,text
;
;
; MODIFICATION HISTORY:
;          2002  - created by PSH
progver = 'v2007-May-19' ;--- (MW) Adapted <demstat.pro> to
;                             <xrt_iter_demstat.pro>.
;
;-
;==========================================================================


   if(NOT KEYWORD_SET(weights)) then weights = $
     replicate(1.0,n_elements(line_data))
   if(NOT KEYWORD_SET(abunds)) then abunds = $
     replicate(1.0,n_elements(line_data))
   dt = t[1]-t[0]
   n_line = N_ELEMENTS(line_data)
   nt = N_ELEMENTS(t)

   ;  p array will contain t*emis*dlnt (when multiplied by dem and summed
   ;  it gives the intensity)

   p = FLTARR(n_line,nt)
   emis = FLTARR(n_line,nt)


   for i=0,n_line-1 do begin
     line_nt = line_data[i].nt
     emis[i,*] = interpol(line_data[i].emis[0:line_nt-1], $
                          line_data[i].t[0:line_nt-1],t    ) > 0.0
     p[i,*] = emis[i,*]*10^t[*]*alog(10^dt)
   endfor

   ; calculate the intensities
   ldem = 10^dem
   i_mod = (ldem##p)*abunds
   di = (i_mod-line_data.i_obs)
   chisq = di^2 * weights / (line_data.i_err)^2
   idx = where(weights eq 0.0,count)
   if(count gt 0) then chisq[idx] = 0.0

   names = line_data.elem+' '+line_data.ion
   text = ' Elem  | Wave | Wgt | Abd | Mod I  | Obs I  |  Err I |   ' $
          + 'di   | chisq '
   fmmt='(a8,T10,f5.1,T16,f5.2,T22,f5.2,T28,G8.3,T37,G8.3,T46,G8.3,' $
        + 'T55,G8.2,T64,G7.2)' 
   for i=0,n_line-1 DO BEGIN
     s = string(names[i],line_data[i].wave,weights[i],abunds[i], $
                i_mod[i],line_data[i].i_obs,line_data[i].i_err,di[i],$
                chisq[i],format=fmmt)
     text = [text,s]
   ENDFOR
   text = [text,'******** Total Chi-squared = '+trim(total(chisq))+ $
           '**********']

   chisq = total(chisq)
   RETURN


END ;======================================================================


;***************************************************************************
;***************************************************************************


FUNCTION mpdemfunct,p,pm=pm,n_line=n_line,i_obs=i_obs,i_err=i_err, $
        i_mod=i_mod,spl_t=spl_t,n_spl=n_spl,t=t,weights=weights, $
        abunds = abunds

;==========================================================================
;+
;
; arguments
; p             dem values at spline knots
; pm[i,*]       emis*T*d(ln t) array for the ith line
; i_obs         observed intensities (scaled)
; i_err         errors in observed intensity
; i_mod         modeled intensities (scaled)
; spl_t         log T values of spline knots
; n_spl         number of spline knots
; t             points to evaluate DEM on
; weights       relative weighting of each line, default = 1.0
; abunds        relative abundance of each line with respect to
;               original value
;
; HISTORY:
progver = 'v2007-May-19' ;--- (M.Weber) Just tweaked a bit to fit into the
;                             organization.
;
;-
;==========================================================================

   dem = 10^(spline(spl_t[0:n_spl-1],p,t))

   i_mod = (dem##pm)*abunds
   chisq = float((i_mod-i_obs)*weights/i_err)
   idx = where(weights eq 0.0,count)
   if(count gt 0) then chisq[idx] = 0.0
   return,chisq


END ;======================================================================


;***************************************************************************
;***************************************************************************


PRO xrt_dem_iter_solver, fit_info, restart=restart, maxiter=maxiter

;==========================================================================
;+
;
;; --- <xrt_dem_iter_solver> calculates the new DEM curve.
progver = 'v2007-May-19' ;--- (M.Weber) Adapted from <dem_display.pro>.
;
;-
;==========================================================================

;=== Preparations ========================================

  ;; Set weights in the function arguments for mpfit  
  fit_info.fa.weights = fit_info.weights * fit_info.use_line
  fit_info.fa.abunds  = fit_info.abundances
  x = fit_info.t
  y = fit_info.dem
  n_spl = fit_info.fa.n_spl

  ;; If /restart set, use existing spline - fit_info must be defined
  if keyword_set(restart) then begin
    ;; Use the current curve as the inital guess for mpfit
    spl_dem = fit_info.spl_dem[0:n_spl-1]
  endif else begin
    ;; Use the estimated (aka flat) DEM as an inital guess for mpfit
    spl_dem = interpol(fit_info.est_dem, x, fit_info.fa.spl_t[0:n_spl-1]) 
  endelse


;=== Do the fit ==========================================

  ;; -- Limit the spline knots to -20 log DEM.

  parinfo = replicate({limited:[1,0], limits:[-20.0,0.D]}, n_spl)


  spl_dem = mpfit('mpdemfunct', spl_dem, functargs=fit_info.fa,    $
                  status=status, errmsg=errmsg, bestnorm=bestnorm, $
                  /quiet, parinfo=parinfo, MAXITER=maxiter)

  ;; Store DEM values for spline knots.
  fit_info.spl_dem[0:n_spl-1] = spl_dem

  ;; Calculate new curve and store.
  fit_info.dem = spline(fit_info.fa.spl_t[0:n_spl-1], spl_dem, x)

  ;; Show results.
  xrt_iter_demstat, fit_info.input_data, fit_info.dem, x, text=text, chisq=chisq, $
          weights=fit_info.weights*fit_info.use_line, $
          abunds=fit_info.abundances
  fit_info.chisq = chisq


;=== Finish ==============================================

  return


END ;======================================================================


;***************************************************************************
;***************************************************************************

pro mp_prep,line_data,weights,abunds,n_spl,spl_dem,fa,$
            est_t=est_t,est_dem=est_dem,min_t=min_t,$
            max_t=max_t,dt=dt

; =========================================================================
;+
;
; Project     : XRT
;
; Name        : mp_prep
;
; Purpose     : To prepare the function arguments needed for the
;               mp_fit call in dem_widget
;
; Explanation : Sets up the initial dem, temperature arrays.  Computes
;               an array of values that when multiplied by the DEM
;               gives intensity
;
; Use         : IDL> mp_prep,line_data,est_t=est_t,$
;                            est_dem=est_dem,dt=dt,weights,abunds,$
;                            n_spl,spl_dem,fa
;
; Inputs      : line_data - an array of structures containing
;                           information about observed lines
;               est_t - log T values for an estimated DEM
;               est_dem - estimated log DEM
;
; Opt. Inputs : weights - initial weights, must have same number of
;                         elements as line_data
;               abunds - intial abundances, must have same number of
;                        elements as line_data
;               n_spl - initial number of spline knots
;               min_t - minimum t value for dem curve
;               max_t - maximum t value for dem curve
;               dt - spacing of t values for dem curve
; Outputs     : weights - the initial weights for each line
;               abunds - initial abundances for each line
;               n_spl - initial number of spline knots
;               spl_dem - initial dem values of spline knots
;               fa - function arguments for use with mp_fit
;
; Opt. Outputs: None
;
; Keywords    : None
;
; Calls       : None
;
; Restrictions: None
;
; Side effects: None
;
; Category    : Util, DEM
;
; Prev. Hist. : None
;
; Written     : P S Hamilton, 30-Oct-02
;
; Modified    : None
;
; Version     : Version 1, 30-Oct-02
;
progver = 'v2007-May-19' ;--- (M.Weber) Not changed much-- same name.
;                             Just a few tweaks to make it fit.
;
;-
; =========================================================================

; Set defaults
  n_line = N_ELEMENTS(line_data)
  if(N_ELEMENTS(n_spl) eq 0) then n_spl = n_line-1 < 7

  if(N_ELEMENTS(weights) eq 0) then weights = replicate(1.0,n_line)
  if(N_ELEMENTS(abunds) eq 0) then abunds = replicate(1.0,n_line)

  if(NOT KEYWORD_SET(dt)) then $
    dt = 0.1

  if(NOT KEYWORD_SET(max_t)) then begin
    max_t = 0.0
    for i=0,n_line-1 do $
      max_t = MAX(line_data[i].t[0:line_data[i].nt-1]) > max_t
  endif
  max_t = round_off(max_t,dt)

  if(NOT KEYWORD_SET(min_t)) then begin
    min_t = 10000.0
    for i=0,n_line-1 do $
      min_t = MIN(line_data[i].t[0:line_data[i].nt-1]) < min_t
  endif
  min_t = round_off(min_t,dt)

  nt = round((max_t - min_t)/dt + 1)
  new_t = min_t + findgen(nt)*dt
  if(N_ELEMENTS(est_dem) eq 0) then $
    new_dem = replicate(22.0,nt) $
  else $
    new_dem = spline(est_t,est_dem,new_t)
  est_dem = new_dem
  est_t = new_t

; p array will contain t*emis*dlnt (when multiplied by dem and summed
; it gives the intensity)
  p = FLTARR(n_line,nt)

  emis = fltarr(n_line,nt)

  for i=0,n_line-1 do begin
    line_nt = line_data[i].nt
    emis[i,*] = interpol(line_data[i].emis[0:line_nt-1], $
                         line_data[i].t[0:line_nt-1],est_t) > 0.0
    p[i,*] = emis[i,*]*10^est_t[*]*alog(10^dt)
  endfor

; Scale the dem to 0-100
;  demscale = MAX(dem) - 2

;  i_obs = line_data.iobs
  i_mod = FLTARR(n_line)

  ; scale everything
;  pscale = max(p) / 10.
;  p = p / pscale
;  i_obs = i_obs / pscale / 10^demscale

; intially set splines at equally spaced intervals
  spl_t = FLTARR(100)
  spl_dem = FLTARR(100)
  spl_t[0:n_spl-1] = min(est_t) + findgen(n_spl) * $
       (max(est_t)-min(est_t)) / (n_spl-1)
  spl_dem[0:n_spl-1] = spline(est_t,est_dem,spl_t[0:n_spl-1])

; set up the function arguments for mpdemfunct
  fa = {pm:p, n_line:n_line, i_obs:line_data.i_obs, i_err:line_data.i_err,$
        i_mod:i_mod, t:est_t, spl_t:spl_t, n_spl:n_spl, weights:weights, $
        abunds:abunds}
return


END ;======================================================================


;***************************************************************************
;***************************************************************************

PRO COMPACT,x,y,newx,newy

; =========================================================================
;+
;
; Project     : XRT
;
; Name        : compact()
;
; Purpose     : Given arrays of x and y values, average y values
;               over any repeated x values
;
;
; Explanation : None
;
; Use         : IDL> compact,x,y,newx,newy
;
; Inputs      : x - array of x values
;               y - array of y values
;
;
; Opt. Inputs : None
;
; Outputs     : newx - compacted x values
;               newy - averaged y values
;
; Opt. Outputs: None
;
; Keywords    : None
;
; Calls       : None
;
; Restrictions: X and Y should have the same number of elements
;
; Side effects: None
;
; Category    : Util, Numerical
;
; Prev. Hist. : None
;
; Written     : P S Hamilton, 28-Oct-02
;
; Modified    : None
;
; Version     : Version 1, 28-Oct-02
;
progver = 'v2007-May-19' ;--- (M.Weber) Added into 
;                             <xrt_dem_iter_nowidget.pro> .
;
;-
; =========================================================================


  done = x & done[*] = 0
  tempx = x & tempy = y
  num = 0
  nx = N_ELEMENTS(x)
  FOR i=0,nx-1 DO BEGIN
    IF(done[i] EQ 0) THEN BEGIN
      idx = WHERE(x EQ x[i])
      tempx[num] = x[i]
      tempy[num] = TOTAL(y[idx])/N_ELEMENTS(idx)
      num = num + 1
      done[idx] = 1
    ENDIF
  ENDFOR
  newx = tempx[0:num-1] & newy = tempy[0:num-1]
  idx = sort(newx)
  newx = newx(idx) & newy = newy(idx)

RETURN


END ;======================================================================


;***************************************************************************
;***************************************************************************


PRO xrt_dem_iter_estim,line_data,t,dem,cutoff = cutoff,verbose=verbose,$
                  nospline=nospline,min_t=min_t,max_t=max_t,dt=dt

; =========================================================================
;+
;
; Project     : XRT
;
; Name        : xrt_dem_iter_estim
;
; Purpose     : To calculate a rough estimate of the DEM assuming
;               given a set of observed line intensities.
;
; Explanation : Calculates a DEM at the peak of emissivity for each
;               line.  Assumes that the DEM is constant over the width
;               of the line and there is no blending.
;
; Use         : IDL> help,line_data,/st
;                 ** Structure MS_089761536002, 10 tags, length=1660:
;                    T               FLOAT     Array[200]
;                    EMIS            FLOAT     Array[200]
;                    NT              INT             11
;                    IOBS            FLOAT           0.00000
;               IDL> xrt_dem_iter_estim,line_data,t,dem
;
; Inputs      : line_data - an array of structures of the form seen above
;
;
; Opt. Inputs : min_t - minimum temperature value of range for DEM to
;                       be spliced onto
;               max_t - maximum temperature value "
;               dt - spacing of temperature values "
; Outputs     : t - An array of the temperature values for these points
;               dem - An array of estimated log DEM values
;
; Opt. Outputs: None
;
; Keywords    : cutoff - the factor at which contributions from the
;                        emissivity should be cutoff, default is 1/e
;               verbose - plot the results
;               nospline - if set do not spline the results onto a
;                          new temperature range, return points only
;
; Calls       : COMPACT, ROUND_OFF
;
; Restrictions: The specified t values in line_data are assumed to
;               be evenly spaced
;
; Side effects: None
;
; Category    : Util, DEM
;
; Prev. Hist. : None
;
; Written     : P S Hamilton, 28-Oct-02
;
; Modified    : None
;
; Version     : Version 1, 28-Oct-02
;
progver = 'v2007-May-19' ;--- (M.Weber) Adapted <xrt_dem_iter_estim.pro>
;                             from <dem_estimator.pro>.
;
;-
; =========================================================================


  IF(NOT KEYWORD_SET(cutoff)) THEN cutoff = 1 / exp(1.0)
  idx = where(line_data.i_obs ne 0.0)
  n_line = N_ELEMENTS(idx)
  t = FLTARR(n_line)
  dem = FLTARR(n_line)

  ;; Assuming that the bulk of the contribution comes from the peak of
  ;; the emissivity, estimate a DEM value for each line

  FOR j=0,n_line-1 DO BEGIN

     l = line_data[idx[j]]
     dlnt = ALOG(10^(l.t[1]-l.t[0]))

     ;; Find the peak temperature of the emissivity
     max_e = MAX(l.emis,max_idx)
     t[j] = round_off(l.t[max_idx],0.1)

     ;; Use values down to where the emissivity = max_e * cutoff
     good = WHERE(l.emis GT max_e * cutoff)

     ;; Integrate over the peak
     eintgrl = TOTAL(10^(l.t[good]) * l.emis[good] * dlnt)

     ;; Determine an average DEM at the peak temperature
     dem[j] = ALOG10(l.i_obs / eintgrl)

  ENDFOR

  ;; remove any duplicate temp points by averaging
  compact,t,dem,t,dem

  if(NOT KEYWORD_SET(nospline)) then BEGIN

    ;; spline onto desired temp values
    if(NOT KEYWORD_SET(dt)) then $
      dt = 0.1

    if(NOT KEYWORD_SET(max_t)) then BEGIN

      ;; Find the max temp of all the lines
      max_t = 0.0
      for i=0,n_line-1 do $
        max_t = MAX(line_data[i].t[0:line_data[i].nt-1]) > max_t

    endif
    max_t = round_off(max_t,dt)

    if(NOT KEYWORD_SET(min_t)) then BEGIN

      ;; Find the min temp for all lines
      min_t = 10000.0
      for i=0,n_line-1 do $
        min_t = MIN(line_data[i].t[0:line_data[i].nt-1]) < min_t
    endif
    min_t = round_off(min_t,dt)

    nt = round((max_t - min_t)/dt + 1)
    new_t = min_t + findgen(nt)*dt
;    new_dem = spline(t,dem,new_t)
;    dem = new_dem
    dem = 0.0*findgen(nt) + 1.0 ; Use flat dem for initial guess
    t = new_t
  endif

  IF(KEYWORD_SET(verbose)) THEN $
    plot,t,dem,xtit='log Temperature', $
       ytit='log DEM',tit='Estimator DEM'
RETURN


END ;======================================================================

;***************************************************************************
;***************************************************************************


PRO xrt_dem_iter_nowidget, input_data, fit_info, max_t=max_t, $
            min_t=min_t, dt=dt, out_t=out_t, out_dem=out_dem, $
            restart=restart, maxiter=maxiter, nosolve=nosolve

; =========================================================================
;+
; NAME: XRT_DEM_ITER_NOWIDGET
;
; PURPOSE:  Calculates DEM given observed intensities, errors, and
; contribution functions.  Was originally written assuming that
; intensity had units of ergs/cm^2/s/ster, DEM had units of 1/cm^5/K,
; and the contribution functions had units of ergs*cm^3/s/ster.
; But, it can be used to solve for DEM in any equation of the form:
;    I = Integral(G(T)*DEM(T)*dT)
;
; Based on dem_widget this is designed to run without user input.
;
; CATEGORY:  Analysis
;
; CALLING SEQUENCE:  
;
;       XRT_DEM_ITER_NOWIDGET, input_data, fit_info, max_t=max_t, $
;                              min_t=min_t, dt=dt, out_t=out_t,   $
;                              out_dem=out_dem, restart=restart
;
; INPUTS:   
;    input_data   - array of structures with the following
;                   fields:
;                     t -   array of temperatures in log T
;                     emis -  emissivity (contribution fxn) at these
;                             temperatures
;                     nt - number of valid t points 
;                     i_obs - observed intensity
;                     elem - element name (for display purposes)
;                     ion - ion Roman numeral (for display purposes)
;                     wave - wavelength of the line (for display only)
;                     ierr - error in i_obs
;
; OPTIONAL INPUTS:
;    max_t  -  maximum temperature (log T) for DEM calculation
;    min_t  -  minimum temperature (log T) for DEM calculation
;    dt     -  temperature step (log T) for DEM calculation
;    restart - Don't use flat initial estimate, start from existing spline
;    nosolve - Don't solve DEM, just fill other values.
;
; OPTIONAL OUTPUTS:
;    out_t  - temperature array (log T) used in DEM calculation
;    out_dem - calculated DEM (log DEM)
;
; MODIFICATION HISTORY:
;    2001   - created by PSH, first major IDL program so I apologize
;             for the disorganization
;    2003-Aug-01 - EED increased maxiter to 400 from default of 200
progver = 'v2007-May-19' ;--- (MW) Adapted <xrt_dem_iter_nowidget.pro>
;                             from <aia_dem_nowidget.pro>.
;
;-
; =========================================================================


;=== Preparations ========================================

  n_ch = n_elements(input_data)

  lines = strarr(n_ch)
  for ii = 0, (n_ch-1) do lines[ii] = string(format='(A," ",A,"-",I0.0)', $
    input_data[ii].elem, input_data[ii].ion, round(input_data[ii].wave))

  uvalue  = strarr(n_ch) + 'LINES'
  weights = fltarr(n_ch) + 1.0
  abundances = fltarr(n_ch) + 1.0

  ;; Only use channels which have a non-zero signal.
  idx = where(input_data.i_obs ne 0.0, cnt_idx)
  use_line = lonarr(n_ch)
  if (cnt_idx ge 1) then use_line[idx] = 1


;=== Fitting info ========================================

  if (not keyword_set(restart)) then begin 

    ;; Get a rough estimate of the DEM.
    if (cnt_idx ge 1) then begin
      xrt_dem_iter_estim, input_data[idx], est_t, est_dem, min_t=min_t, $
                          max_t=max_t, dt=dt
    endif else begin
      est_t = input_data.t
      est_dem = est_t * 0
    endelse

    ;; Set up the function arguments for <mpfit>.
    mp_prep, input_data, weights, abundances, n_spl, spl_dem, fa, $
             est_t=est_t, est_dem=est_dem, min_t=min_t, max_t=max_t, dt=dt

    ;; --- structure to hold program info 
    fit_info = { input_data: input_data,  $ ; the line observations
                 lines:      lines,       $ ; names of the lines
                 t:          est_t,       $ ; temperature array
                 est_dem:    est_dem,     $ ; estimated DEM
                 dem:        est_dem,     $ ; actual DEM
                 spl_dem:    spl_dem,     $ ; DEM at spline knots
                 weights:    weights,     $ ;
                 abundances: abundances,  $ ;
                 use_line:   use_line,    $ ; current state of check 
                                            ;   boxes for the lines
                 fa:         fa,          $ ; function arguments for mpfit
                 chisq:      0.0          } ; Return the chisq of fit
  endif


;=== Do fitting ==========================================

  if (not keyword_set(nosolve))                                        $
    then xrt_dem_iter_solver, fit_info, restart=restart, maxiter=maxiter


;=== Finish ==============================================

  out_t = fit_info.t

  if (keyword_set(nosolve))         $
    then out_dem = fit_info.t * 0   $
    else out_dem = fit_info.dem

  return


END ;======================================================================


;***************************************************************************
;***************************************************************************

; wrapperstartshere

PRO xrt_dem_iterative2, obs_index, obs_val, tresp, logT_out, dem_out,  $
                        obs_err=obs_err, base_obs=base_obs,            $
                        mod_obs=mod_obs, chisq=chisq,                  $
                        min_T=min_T, max_T=max_T, dT=dT,               $
                        MC_iter=MC_iter, solv_factor=solv_factor,      $
                        maxiter=maxiter, quiet=quiet, verbose=verbose, $
                        qabort=qabort, qstop=qstop

; =========================================================================
;+
; PROJECT:
;       Solar-B / XRT
;
; NAME:
;
;       XRT_DEM_ITERATIVE2
;
; CATEGORY:
;
;       DEM (Differential Emission Measures)
;
; PURPOSE:
;
;       WARNING: This routine replaces <xrt_dem_iterative.pro>. 
;       See Note #6.
;
;       Estimate a DEM(T) curve, given some observations B_i in 
;       channels "i", and given the temperature response functions
;       in every channel R_i(T). These functions satisfy the equation:
;            B_i = integral{ DEM(T) * R_i(T) * dT }
;
;       The inversion is ill-posed and technically fraught with perils.
;       This routine employs a forward-fitting approach: A DEM is guessed
;       and folded through the R_i(T) to generate "model" observations.
;       This process is iterated to reduce the chi-square between the
;       actual and model observations. The DEM function is interpreted
;       from some spline points, which are directly manipulated by the 
;       chi-square fitting routine (MPFIT.pro). There are N_i - 1 
;       splines, representing the degrees of freedom for N_i observations.
;       (Note that the number of temperature bins requested for the DEM
;       solution are usually greater than N_i.)
;
;       To estimate errors on the DEM solution, this routine provides for
;       Monte-Carlo iteration. On each iteration, the observations are
;       varied normally by their sigma error, and then solved for a DEM.
;       According to Monte Carlo theory, the distribution of DEM solutions
;       is a measure of the error in DEM(T). 
;
;       This routine only works on one pixel x N channels, at a time.
;
; CALLING SEQUENCE:
;
;       XRT_DEM_ITERATIVE2, obs_index, obs_val, tresp, logT_out, dem_out
;               [,obs_err=obs_err] [,base_obs=base_obs] [,mod_obs=mod_obs]
;               [,chisq=chisq]     [,min_T=min_T]
;               [,max_T=max_T]     [,dT=dT]             [,maxiter=maxiter]
;               [,MC_iter=MC_iter] [,solv_factor=solv_factor]  [,/verbose]
;               [,/quiet]          [,/qstop]              [,qabort=qabort]
;
; INPUTS:
;
;       OBS_INDEX   - [Mandatory] (structure or string array, [Nobs])
;                     This parameter identifies the XRT x-ray channel
;                     for the corresponding element of the OBS_VAL array.
;                     OBS_INDEX may take one of two forms:
;                     a) The "index" or "catalog" structure for the
;                        image from which this pixel OBS_VALUE came from.
;                     b) A string with the name of the XRT x-ray channel
;                        or channel temperature response that corresponds
;                        to this OBS_VALUE. Here are the possible names
;                        for the default XRT channel names:
;                        {Al-mesh,  Al-poly, C-poly, Ti-poly, Be-thin, 
;                         Be-med, Al-med, Al-thick, Be-thick,
;                         Al-poly/Al-mesh, Al-poly/Ti-poly, 
;                         Al-poly/Al-thick, Al-poly/Be-thick          }
;       OBS_VAL     - [Mandatory] (float array, [Nobs])
;                     This is the set of XRT data values seen in the 
;                     x-ray channels identified by OBS_INDEX, respectively.
;                     Units = 'DN s^-1 pix^-1', where "pix" means a
;                     one-arcsecond,]            full-resolution XRT pixel. 
;       TRESP       - [Mandatory] ("temp_resp" structure array, [Nchn])
;                     This structure contains data and information
;                     about the temperature responses for a set of XRT
;                     x-ray channels. A "channel" is primarily specified
;                     by a particular setting of the x-ray analysis
;                     filters and a CCD contamination thickness. This
;                     structure is obtained using <make_xrt_wave_resp.pro>
;                     and then <make_xrt_temp_resp.pro>. More information
;                     can be found in the headers of those programs
;                     and in the XRT Analysis Guide (located in the
;                     Solarsoft directory $SSW_XRT/docs/XAG.pdf and at
;           http://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf ).
;
; KEYWORDS:
;
;       OBS_ERR     - [Optional] (float array, [Nobs]) 
;                     This input may provide the Gaussian one-sigma errors
;                     for the values in OBS_VAL. The solver requires errors
;                     in order to calculate the least-squares. If the
;                     user does not provide these, then the default will
;                     be used. Default = (OBS_VALUE * 0.03). (Min = 2 DN.)
;       MIN_T       - [Optional] (float scalar)
;                     Input the low end of the DEM temp. range.
;                     Units = 'log K'. Default = 5.5. See Note #2.
;       MAX_T       - [Optional] (float scalar)
;                     Input the high end of the DEM temp. range.
;                     Units = 'log K'. Default = 8.0. See Note #2.
;       DT          - [Optional] (float scalar)
;                     Input the bin-width of the DEM temp. range.
;                     Units = 'log K'. Default = 0.1. See Note #2.
;       MAXITER     - [Optional] (float scalar)
;                     This program works by iterating a least-squares
;                     search. This keyword may be used to specify the 
;                     maximum number of iterations for each DEM solution.
;                     Default = 2000.
;       MC_ITER     - [Optional] (float scalar)
;                     Use this keyword to cause the program to perform 
;                     Monte Carlo runs. The value of MC_ITER indicates
;                     how many runs. For each MC run, each of the OBS_VAL 
;                     values are modified from their original value by a 
;                     random normal amount using the respective OBS_ERR as 
;                     a Gaussian sigma. Then, the DEM is solved again using
;                     the new set of OBS_VAL. See the descriptions of the
;                     output variables for how the results are reported.
;                     No MC runs are performed if this keyword is not used.
;       SOLV_FACTOR - [Optional] (float scalar)
;                     The least-squares solver is not completely 
;                     insensitive to the order of magnitude of the numbers
;                     it is manipulating. SOLV_FACTOR is used to 
;                     normalize the inputs to move the solver into a 
;                     "sweet spot". The default choice is arbitrary but 
;                     seems to work well (default = 1e21). Of course,
;                     the outputs are un-normalized at the end. It is
;                     recommended that users use the default.
;       /VERBOSE    - [Optional] (Boolean) If set, print out extra
;                     information. Overrides "/quiet" (see Note #1).
;       /QUIET      - [Optional] (Boolean) If set, suppress messages
;                     (see Note #1).
;       /QSTOP      - [Optional] (Boolean) For debugging.
;
; OUTPUTS:
;
;       LOGT_OUT    - [Mandatory] (float array, [Ntemp])
;                     The temperatures corresponding to the DEM solution.
;                     Units = 'log K'. Runs from MIN_T to MAX_T with 
;                     bin-width = DT.
;       DEM_OUT     - [Mandatory] (float array, [Ntemp, 1+MC_ITER])
;                     A DEM solution is a 1D array of length = Ntemp.
;                     If Monte Carlo runs were performed, then the
;                     the 2nd dimension spans these runs.
;                     DEM_OUT[*,0] is the solution for the original
;                     OBS_VAL. DEM_OUT[*,1:MC_ITER] correspond to the
;                     solutions for the MC runs. Units = 'cm^-5 K^-1'.
;                     See Note #3 for more about units.
;       BASE_OBS    - [Optional] (float array, [Nobs, 1+MC_ITER])
;                     These are the observations that a DEM solution
;                     corresponds to. BASE_OBS[*,0] = OBS_VAL.
;                     BASE_OBS[*,1:MC_ITER] correspond to the observations
;                     produced by adding random OBS_ERR noise to OBS_VAL
;                     for each MC run. You may think of these as the 
;                     "actual" set of observations that each MC run was
;                     trying to solve the DEM for.
;       MOD_OBS     - [Optional] (float array, [Nobs, 1+MC_ITER])
;                     These are the model observations which are produced
;                     by the corresponding DEM solution. BASE_OBS is the
;                     actual set, and MOD_OBS is how close the DEM can get
;                     to matching it.
;       CHISQ       - [Optional] (float array, [1+MC_ITER])
;                     These are the chi-square values for the DEM 
;                     solutions. 
;               CHISQ[i] = total (BASE_OBS[*,i]-MOD_OBS[*,i]/OBS_ERR[*])^2
;       QABORT      - [Optional] (Boolean) Indicates that the program
;                     exited gracefully without completing. (Might be
;                     useful for calling programs.)
;                     0: Program ran to completion.
;                     1: Program aborted before completion.
;
; EXAMPLES:
;
;       Simplest possible usage:
;       IDL> help, obs_index, obs_val
;       OBS_INDEX       STRING    = Array[13]
;       OBS_VAL         DOUBLE    = Array[13]
;       IDL> xrt_dem_iterative, obs_index, obs_val, tresp, logT_out, dem_out     
;       IDL> help, logT_out, dem_out
;       LOGT_OUT        FLOAT     = Array[26]
;       DEM_OUT         DOUBLE    = Array[26]
;       IDL> plot, logT_out, dem_out, psym=10, xtit='log T [log K]', $
;                  ytit='DEM [cm^-5 K^-1]'
;
;       Provide errors (obs_err) and perform 100 Monte Carlo runs:
;       IDL> help, obs_index, obs_val, obs_err
;       OBS_INDEX       STRING    = Array[13]
;       OBS_VAL         DOUBLE    = Array[13]
;       OBS_ERR         DOUBLE    = Array[13]
;       IDL> xrt_dem_iterative, obs_index, obs_val, tresp, logT_out, dem_out, $
;       IDL> obs_err=obs_err, base_obs=base_obs, mod_obs=mod_obs, $
;       IDL> chisq=chisq, MC_iter=100
;       IDL> help, logT_out, dem_out
;       LOGT_OUT        FLOAT     = Array[26]
;       DEM_OUT         DOUBLE    = Array[26]
;       IDL> help, base_obs, mod_obs, chisq
;       BASE_OBS        DOUBLE    = Array[13, 101]
;       MOD_OBS         DOUBLE    = Array[13, 101]
;       CHISQ           FLOAT     = Array[101]
;
; COMMON BLOCKS:
;
;       none
;
; NOTES:
;
;       1) There are three levels of verbosity.
;          a) "verbose" = highest priority. All errors and messages are
;                         displayed. ("if q_vb")
;          b) "quiet" = lower priority. No errors or messages are
;                       displayed. ("if q_qt")
;          c) neither = lowest priority. All errors and some messages are
;                       displayed. ("if not q_qt")
;
;       2) The user may specify the temperature bins over which the
;          DEM will be solved. However, this solution temperature range 
;          must lie within all of the T-ranges of the temperature 
;          responses, to permit interpolations. (The default XRT responses
;          run logT = 5.5 to 8.0.) Also, the solver requires that the 
;          T-range be regular in 'log T' units. Therefore, the solution 
;          range is controlled with these three keywords: MIN_T, MAX_T, 
;          and DT. Although these temperatures are all presented in 
;          'log T', the underlying integrals are performed over 'T'.
;
;       3) A discussion of units.
;          It is easier to understand these units in the larger context.
;          For a spectral model S, such as is returned by
;          <get_xrt_spec_genx>:  S(wave, T) ~ [ph cm^3 s^-1 sr^-1 A^-1].
;          For a spectral response R, such as <calc_xrt_spec_resp.pro> 
;          returns:
;          R(wave) ~ [DN cm^2 sr ph^-1 pix^-1], which says something about
;          how many DataNumbers are generated in the camera for a photon of
;          a given wavelength. Then:
;          S(wave, T) * R(wave) ~ [DN A^-1 s^-1 pix^-1 (cm^-5)^-1]
;                               ~ [DN A^-1 s^-1 pix^-1 EM^-1],
;          where EM ~ [cm^-5] is the "line of sight" (or "column")
;          emission measure. For a channel's temperature response F{T):
;          F(T) = integrate{ S(wave, T) * R(wave) * d(wave)}
;               ~ [DN s^-1 pix^-1 EM^-1]
;          Note that one must assume a spectral model to calculate a
;          temperature response.
;
;          Going just a bit further, for an EM at T = T0, the signal G
;          observed is
;          G|T0 = F(T0) * EM|T0 ~ [DN s^-1 pix^-1].
;          For an emission measure continuously distributed over a range
;          of temperatures T, one uses the differential emission measure:
;          DEM(T) ~ [cm^-5 K^-1].
;          Now the net signal G is
;          G = integrate{F(T) * DEM(T) * dT} ~ [DN s^-1 pix^-1].
;
;          Note that this program <xrt_dem_iterative.pro> just solves
;          the equation:
;            B_i = integral{ DEM(T) * R_i(T) * dT }
;          So DEM(T) will satisfy the unit-equation:
;            U{DEM} = U{G} / U{F} / U{dT}.
;          As long as the units are all consistent, it is less important
;          what they are. So if the temperature response units were
;          'DN cm^5 s^-1 pix^-1', and the observations were in units of
;          'DN s^-1 pix^-1', then DEM will still have units of 'cm^-5 K^-1'.
;          The moral of the story is that a combination of temperature
;          responses and observations can be solved, given sufficient care
;          about units and cross-calibrations.
;
;       4) There is a program called <xrt_dem2obs.pro>, which may be used
;          to generate a set of XRT observations from a DEM(T) curve.
;          (See that program for an explanation of its usage.)
;
;          We can use these faux observations to illustrate the usage of 
;          the solver, since you can compare the "real" original DEM 
;          against the estimated DEM. Here is a simple tutorial using both 
;          programs. First, create a simple DEM curve.
;          IDL> demt = findgen(26)*0.1 + 5.5 ;; Temps are logarithms.
;          IDL> dem = demt * 0
;          IDL> dem[15] = 1.3e22
;          IDL> print, demt[15]
;                7.00000
;          IDL> wdef, 0, 1000
;          IDL> plot, demt, dem, psym=10, /ytype, yr=[1e20,1e23]
;
;          Now generate the set of observations XRT would see.
;          IDL> obs = xrt_dem2obs(demt, dem)
;          IDL> help, obs                      ;; 14 XRT channels
;          OBS             STRUCT    = -> <Anonymous> Array[14]
;          IDL> help, obs[0], /st
;          ** Structure <26648d4>, 3 tags, length=28, data length=28, refs=1:
;             CHANNEL_NAME    STRING    'Al-mesh'
;             OBS             FLOAT           8344.32
;             OBS_UNITS       STRING    'DN s^-1 pix^-1'
;
;          Just to be clear, prepare input variables for the solver.
;          IDL> obs_val = obs.obs
;          IDL> obs_index = obs.channel_name
;          IDL> help, obs_index, obs_val
;          OBS_INDEX       STRING    = Array[14]
;          OBS_VAL         DOUBLE    = Array[14]
;
;          Run the solver.
;          IDL> xrt_dem_iterative, obs_index, obs_val, logT_out, dem_out
;          IDL> help, logT_out, dem_out
;          LOGT_OUT        FLOAT     = Array[26]
;          DEM_OUT         DOUBLE    = Array[26]
;
;          Overplot the solution.
;          IDL> oplot, logT_out, dem_out, psym=10, linesty=2
;
;          As a follow-up exercise, re-run the solver for 10 Monte Carlo
;          iterations and overplot the results.
;
;       5) When the OBS_VAL values are smaller than the OBS_ERR 
;          uncertainties, there is a chance that the Monte Carlo
;          random variations to the data set will generate a set 
;          of all-zero observations. The program checks for this case,
;          and returns a DEM function equal to zero.
;
;       6) This routine deprecates and replaces the older routine
;          <xrt_dem_iterative.pro>. The older routine is not
;          compatible with the most recent response routines.
;
; CONTACT:
;
;       Comments, feedback, and bug reports regarding this routine may be
;       directed to this email address:
;                xrt_manager ~at~ head.cfa.harvard.edu
;
; MODIFICATION HISTORY:
;
progver = 'v2007-May-21' ;--- (MW) Written, (Based heavily on
;                             <aia_dem_err.pro> and associated routines.)
progver = 'v2007-Aug-09' ;--- (MW) Fixed bug with MAX_T.
progver = 'v2007-Aug-15' ;--- (MW) Fixed crash when Monte Carlo generates
;                             a set of all-zero observations. Also added
;                             a reform() on OBS_VAL to fix a bug.
progver = 'v2008-Oct-05' ;--- (MW) Developed <xrt_dem_iterative2> to 
;                             replace <xrt_dem_iterative>. See Note #6.
;
;-
; =========================================================================


; === Initial setup ======================================

  ;=== Initialize program constants.
  prognam = 'XRT_DEM_ITERATIVE2'
  prognul = '                  '

  ;=== Set Booleans which control print statements.
  ;=== Keyword "verbose" overrules "quiet".
  q_vb = keyword_set(verbose)
  q_qt = keyword_set(quiet) and (q_vb eq 0)

  ;=== Announce version of <read_xrt.pro>.
  if q_vb then box_message, prognam+': Running ' + progver + '.'
  if q_vb then print, prognam+': Performing initial setup...'

  ;=== Set some keyword Booleans and defaults.
  default, min_T,       5.5                       ;; log K
  default, max_T,       8.0                       ;; log K
  default, dT,          0.1                       ;; log K
  default, solv_factor, 1d21
  default, maxiter,     2000
  default, obs_err,     (obs_val*0.03) > 2d       ;; DN/s/pix
  q_MC     = keyword_set(MC_iter)
  qstop    = keyword_set(qstop)
  qabort   = 0B

  ;=== Other useful stuff.
  n_obs = n_elements(obs_val)

  if q_vb then print, prognul+'                          ...OK'


;=== Check inputs ========================================

  if q_vb then print, prognam+': Checking inputs...'

  ;=== Check for number of inputs.
  if (n_params() ne 5) then begin
    if (not q_qt) then box_message, prognam+': Incorrect number of ' $
                                    + 'parameters. Aborting.'
    qabort = 1B
    return
  endif

  ;=== Compare obs_index with obs_val.
  if (n_elements(obs_index) ne n_elements(obs_val)) then begin
    if (not q_qt) then box_message, [prognam+': OBS_INDEX does not ' $
          + 'have the same length as OBS_VAL.', prognul+'  Aborting.']
    qabort = 1B
    return
  endif 

  ;=== Check type of obs_index.
  ind_typ = datatype(obs_index)
  case 1 of
    (ind_typ eq 'STR') : begin
      obs_index1 = obs_index
    end

    (ind_typ eq 'STC') : begin
      obs_index1 = xrt_index2chname(obs_index, qabort=qabort)
      if qabort then begin
        if (not q_qt) then box_message, prognam+': Cannot resolve ' + $
                              'channel names from OBS_INDEX. Aborting.'
        qabort = 1B
        return
      endif
    end

    else               : begin
      if (not q_qt) then box_message,                                $
        [prognam+': OBS_INDEX must be a structure or string array.', $
         prognul+'  Please refer to the header description. Aborting.']
      qabort = 1B
      return
    end
  endcase

  ;=== Check type of obs_val.
  if (total(is_number(obs_val)) ne n_obs) then begin
    if (not q_qt) then box_message, $
      [prognam+': OBS_VAL must be a number array. Aborting']
    qabort = 1B
    return
  endif 
  obs_val = double(obs_val)

  ;=== Check T-response selection.
  trtyp = datatype(tresp)
  if (trtyp eq 'STC') then begin
    if not required_tags(tresp, 'TYPE,NAME,TEMP,TEMP_UNITS,' $ 
                        + 'TEMP_RESP,TEMP_RESP_UNITS,LENGTH') then begin
      if (not q_qt) then box_message, [prognam+': Do not recognize' $
          + ' the TRESP structure. Aborting. (See', $
          prognul+'  <make_xrt_temp_resp.pro> for more info ' $
          + 'about this structure.)']
      qabort = 1B
      return
    endif
  end
  if q_vb then begin
    print, prognam+': Here is the list of valid channel names from the'
    print, prognul+'  loaded temp. responses:'
    print, prognul+'  ', tresp.confg.name
  endif

  ;=== Sort tresp against obs_index1.
  ss_obs = 1L
  for ii = 0,(n_obs-1) do begin
    ss = where(tresp.confg.name eq obs_index1[ii], count)
    if (count eq 0) then begin
      if (not q_qt) then box_message, $
        [prognam+': OBS_INDEX[' + strcompress(ii,/rem) + '] = "' + $
                obs_index1[ii] + '" does not match any of the', $
         prognul+'  loaded temp. responses. Aborting.', $
         prognul+'  (Run again with /verbose to see list of valid names.)']
      qabort = 1B
      return
    endif
    ss_obs = [ss_obs, ss[0]]
  endfor
  ss_obs = ss_obs[1:*]
  tresp1 = tresp[ss_obs]

  ;=== Check temperature range.
  min_T = min(min_T)    ;; Intended solution range
  max_T = max(max_T)    ;; Intended solution range
  dT    = float(dT[0])  ;; Intended solution range
  ntemp = long(((max_T-min_T)/dT)+1)
  logT_out = findgen(ntemp) * dT + min_T

  ;; Find most constrictive limits in all temp_resps.
  tresp_minT = alog10(max(tresp.temp[0]))
  tresp_maxT = alog10(min(max(tresp.temp, dim=1)))
  tresp_range = '[' + strcompress(tresp_minT,/rem) + ', ' + $
                strcompress(tresp_maxT,/rem) + ']'

  ;; Can't interpolate beyond temp_resps.
  if ((min_T lt tresp_minT) or (max_T gt tresp_maxT)) then begin
    if (not q_qt) then box_message,                                       $
      [prognam+': MIN_T or MAX_T is outside the interpolation range for', $
       prognul+'  these temp. responses ' + tresp_range + '.',            $
       prognul+'  Please pick appropriate values. Aborting.']
    qabort = 1B
    return
  endif

  ;=== Make sure MC_iter is useful.
  if q_MC then MC_iter = long(MC_iter[0])

  if q_vb then print, prognul+'                 ...OK'


;=== Prepare input for xrt_dem_iter_nowidget =============

  if q_vb then print, prognam+': Preparing for the solver...'

  input_str = { t   : fltarr(200), emis   : dblarr(200), $
               nt   :   0,         i_obs  : 0d,          $
               i_err:  0d,         i_units: '',          $
               elem :  '',         ion    :  '',         $
               wave : 0.0,         ipf    : ''            }

  input_str.t[0:ntemp-1] = logT_out
  input_str.nt           = ntemp

  input1 = replicate(input_str, n_obs)

  for obscnt = 0,(n_obs-1) do begin
    tlen = tresp1[obscnt].length
    emis1 = interpol(tresp1[obscnt].temp_resp[0:tlen-1], $
                     alog10(tresp1[obscnt].temp[0:tlen-1]), logT_out)
    input1[obscnt].emis[0:ntemp-1] = emis1
    input1[obscnt].i_obs = obs_val[obscnt]
    input1[obscnt].i_err = obs_err[obscnt]
  endfor

  ;; The factor is necessary for the DEM to be of a good order
  ;; of magnitude for the least-squares solver.
  input1.i_obs = input1.i_obs / solv_factor
  input1.i_err = input1.i_err / solv_factor

  if q_vb then print, prognul+'                          ...OK'


;=== Run solver ==========================================

  if q_vb then print, prognam+': Run the solver for the base case...'

  base_obs = reform(obs_val)

  ;; If the number of observations "not equal to zero" is zero,
  ;; then all observations are equal to zero,
  ;; and so we just want to set "nosolve = 1" and dem_out = 0...
  ss_noz = where(input1.i_obs ne 0.0, cnt_noz)
  nosolve = (cnt_noz eq 0)

  xrt_dem_iter_nowidget, input1, fit_info, min_t=min_t, max_t=max_t, $
                         dT=dT, out_t=out_t, out_dem=out_dem,        $
                         maxiter=maxiter, nosolve=nosolve
  logT_out = out_t
  dem_out = 10.^out_dem * solv_factor

  chisq = mpdemfunct(fit_info.spl_dem, pm=fit_info.fa.pm,              $
                     spl_t=fit_info.fa.spl_t, n_spl=fit_info.fa.n_spl, $
                     t=fit_info.t, abunds=fit_info.abundances,         $
                     i_obs=fit_info.input_data.i_obs,                  $
                     weights=fit_info.weights,                         $
                     i_err=fit_info.input_data.i_err, i_mod=i_mod)
  chisq = total(chisq^2)
  mod_obs = i_mod * solv_factor

  if q_vb then print, prognul+'                                  ...OK'


;=== Run the Monte Carlo loops ===========================

  if q_MC then begin

    if q_vb then print, prognam+': Run the solver for M.C. iterations...'

    seed = systime(1)

    for ii = 1,MC_iter do begin
      if (q_qt eq 0) then print, prognam+': M.C. loop ' +              $
         strcompress(ii,/rem) + ' of ' + strcompress(MC_iter,/rem) + '.'
      input2 = input1
      input2.i_obs =                                           $
        (input1.i_obs + randomn(seed, n_obs)*input1.i_err) > 0.0
      base_obs = [[base_obs], [input2.i_obs * solv_factor]]

      fit_info = 0

      ;; If the number of observations "not equal to zero" is zero,
      ;; then all observations are equal to zero,
      ;; and so we just want to set "nosolve = 1" and dem_out = 0...
      ss_noz = where(input2.i_obs ne 0.0, cnt_noz)
      nosolve = (cnt_noz eq 0)

      xrt_dem_iter_nowidget, input2, fit_info, min_t=min_t, max_t=max_t, $
                             dT=dT, out_t=out_t, out_dem=out_dem,        $
                             maxiter=maxiter, nosolve=nosolve

      if nosolve then dem_out = [[dem_out], [out_dem]]                 $
                 else dem_out = [[dem_out], [10.^out_dem * solv_factor]]

      chisq2 = mpdemfunct(fit_info.spl_dem, pm=fit_info.fa.pm,            $
                        spl_t=fit_info.fa.spl_t, n_spl=fit_info.fa.n_spl, $
                        t=fit_info.t, abunds=fit_info.abundances,         $
                        i_obs=fit_info.input_data.i_obs,                  $
                        weights=fit_info.weights,                         $
                        i_err=fit_info.input_data.i_err, i_mod=i_mod)
      chisq   = [chisq, total(chisq2^2)]
      mod_obs = [[mod_obs], [i_mod * solv_factor]]

    endfor

    if q_vb then print, prognul+'                                    ...OK'

  endif


;=== Finish ==============================================

  if q_vb then print, prognam+': Finished.'

  if qstop then stop


END ;======================================================================
