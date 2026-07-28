FUNCTION VDEM,flux_parameters,ION_IN=ion_in,RETFLUX=retflux,CAXIX=CAXIX,$
	DVEL=dvel,WEIGHT=weight,XI=xi,VERBOSE=verbose,specs=specs,$
	width_sumer=width_sumer,ref_file=ref_file

; NAME:

;	vdem.pro

; PURPOSE:

;	This is the main routine in the package of routines associated with
;	calculation of Velocity Differential Emission Measure (VDEM)
;	[photons/s/(cm/s)] from observed spectral line.  See
;	Newton, Emslie, & Mariska,(1995, ApJ., 447:915) for derivation
;	of concept.  This routine can process both SUMER and BCS data.

; CATEGORY:

;	Data analysis

; CALLING SEQUENCE:

;	vdem_parameters=VDEM(flux_parameters[,ION_IN=ion_in,RETFLUX=retflux,$
;		DVEL=dvel,WEIGHT=weight,XI=xi,spec=specs,$
;		width_sumer=width_sumer,ref_file=ref_file,/VERBOSE,/CAXIX])

; INPUTS:

;	flux_parameters= structure returned from vdem processing program
;		(PREP_VDEM_SUMER or PREP_VDEM_BCS).  Includes flux,
;		error in flux, and wavelength array.

;	(The user will be prompted for the following if ION_IN keyword
;	is not set:)

;		mass= atomic weight from periodic table (in units of
;			proton mass)
;		rest wavelength= measured rest wavelength in lab frame
;			(Angstroms)
;		temperature= temperature where ion is formed (degrees
;			Kelvin)

; OPTIONAL INPUTS:

;	ion_in= array containing [mass (units of proton mass),
;		rest wavelength (Angstroms), temperature (K)]
;		For example, if deconvolving VDEM from an oxygen line with
;		rest wavelength of 1032.0 and formed at a temperature
;		of 3.0e5, the ion_in array would be:
;			ion_in=[16.0,1032.0,3.0e5]
;	dvel= velocity discretization, default is the observational limit:
;		dvel=c*(wavelength bin width/rest wavelength) (cm/s)
;	weight=smoothing parameter used to caculate VDEM.  Default uses
;		GCV technique to calculate optimal weighting.  Otherwise
;		suggested ranges are [0-1.e4]
;	xi=nonthermal velocity (cm/s). Default is 0.
;	specs= array containing [spec_first,spec_last] where spec_first
;		is the first spectrum the user wishes vdem to analyze
;		and spec_last is the last spectrum the user wishes vdem to
;		analyze.  If not set, VDEM will analyze all submitted spectra.
;	width_sumer= the width of instrumental broadening function for SUMER in Angstroms.
;		If specified, this width will be added to the kernel function
;		and deconvolved from the flux.  This input should not be
;		set for BCS data because the instrumental function is
;		non-gaussian.
;	ref_file=filename containing excitation coefficients.  Default
;		is $VDEM/excitcoef.dat.  This is only employed if
;		the CAXIX keyword is set.

; OUTPUTS:

;	vdem_parameters= structure containing returned vdem, error in
;		vdem, and velocity array.

; OPTIONAL OUTPUTS:

;	retflux= flux calculated from forward convolution of inverted vdem.

; KEYWORDS:

;	/CAXIX= Should be set when examining the Ca XIX line in Yohkoh data
;		(see Restrictions)
;	/verbose=user will be given options of viewing fluxes, vdems, and
;		returned fluxes.

; COMMON BLOCKS:

;	None.

; SIDE EFFECTS:

;	If /VERBOSE keyword is set, the program may set plot parameters
;	(such as !p.region) to default values.

;	If user specifies the optional input, weight=0, then routine
;	defaults to employing the GCV technique for determining the
;	optimal smoothing parameter.

; RESTRICTIONS:

; 	The input, flux_parameters, of this routine must be calculated
; 	using either PREP_VDEM_SUMER or PREP_VDEM_BCS.

;	When Ca XIX spectral data are used, routine reads an ascii
;	file, 'exitcoef.dat' or user-defined file, containing a
;	4 x 11 array ('coef') of the
;	excitation coefficient for the w line, d13 line, and
;	dielectric recombination contribution to the w line.  Data taken from
;	Bely-Dubau et al. 1982.  This file can be edited if user
;	wishes to use different values.

;	Users are cautioned that the deconvolved VDEM's behavior at its
;	endpoints is poorly constrained by the data, and hence its endpoint
;	values will exhibit artifacts which typically arise in inverse
;	problems.

;	VDEM routines have only been validated for SUMER data and the
;	BCS CaXIX and SXV data.  Inversions of BCS Fe lines should be
;	approached with caution because of the presence of large
;	satellite lines and low flux levels.

; ERRORS:

; 	There are several error messages that may occur during the
;	running of vdem.pro.  Most fundamentally, the data analysis
;	window should only contain points with substantial counts.
;	The flux and errors in flux must not be negative (or zero
;	in the case of the errors).

;	If the error states: 'inverse of
; 	the matrix is suspect', the inversion error is significant
;	and the returned vdem should not be trusted.
;	If the error states: 'the returned vdem is negative',
;	user should attempt to correct
;	situation by subtracting less background in the prep_vdem routines
;	 or narrowing the
; 	data analysis window (pixels or wavelength range).  Inversion
;	results in a  slightly oscillatory solution, so if the
;	analysis window is too large, the solution  may drift below 0  where
;	the vdem should instead be 0.  If user has specified a smoothing
;	parameter which results in a negative VDEM, the user may correct
; 	the problem by using a smaller smoothing parameter or using the
;	default values calculated from the GCV method.

; PROCEDURE:

;	Before running VDEM, user must add the vdem directory to the IDL path
;	in order to utilize the accompanying data file.  For example, the user
;	will add the following line to idl_startup.pro to designate the location
;	of the VDEM software installation:
;			setenv,'VDEM=/where_vdem_routines_are/'.

; 	First run either PREP_VDEM_SUMER or PREP_VDEM_BCS to arrive
;	at flux_parameters.  (For more information on these routines,
;	see documentation.)  VDEM may then be calculated.

; EXAMPLE: (after preparing flux_parameters for a data set containing
;	51 flux bins and 11 spectra.)

; IDL> help,flux_parameters,/str
;** Structure <426e90>, 3 tags, length=4692, refs=2:
;    FLUX            FLOAT     Array(51, 11)
;    EFLUX           FLOAT     Array(51, 11)
;    WAVE            FLOAT     Array(51)

; IDL> vdem_parameters=vdem(flux_parameters)

; Input atomic weight from periodic table (units of proton mass): 16.

; Input laboratory rest wavelength of line (Angstroms):1031.9261

; Input average temperature where line is formed (K):2.88e5

; Mass:      16.0000
; Rest Wavelength:      1031.93
; Temperature:      288000.
; Is this correct? (y/n): y

; IDL> help,vdem_parameters,/str
;** Structure <956558>, 3 tags, length=9200, refs=1:
;    VDEM            DOUBLE    Array(50, 11)
;    EVDEM           DOUBLE    Array(50, 11)
;    VELOCITY        DOUBLE    Array(50)

; To plot the returned VDEM of the 4th spectrum vs velocity with errors:

; IDL> plot_err,vdem_parameters.velocity,vdem_parameters.vdem(*,4),$
; IDL> vdem_parameters.evdem(*,4)

; MODIFICATION HISTORY:

;	E. Newton 1994:	Original version
;	A. Winebarger 8/97: Changed to accept SUMER data
;	Newton & Winebarger 6/98: Final preparation for SolarSoft software tree


;........................Get input data..............................


IF KEYWORD_SET(ION_IN) THEN BEGIN
	mass=ion_in(0)
	L0=ion_in(1)
	tbar=ion_in(2)
ENDIF ELSE BEGIN

START:
mass=0.
PRINT,''
READ,mass,$
    PROMPT="Input atomic weight from periodic table (units of proton mass):"
PRINT,''
L0=0.
PRINT,''
READ,L0,PROMPT="Input laboratory rest wavelength of line (Angstroms):"
PRINT,''
tbar=0.
PRINT,''
READ,tbar,PROMPT="Input average temperature where line is formed (K):"
PRINT,''

PRINT,'Mass:',mass
PRINT,'Rest Wavelength:',L0
PRINT,'Temperature:',tbar

answer=''
PRINT,''
READ,answer,PROMPT="Is this correct? (y/n): "
PRINT,''
IF (answer NE 'y') THEN GOTO,START

ENDELSE
;.......................Set Constants.................................

k=1.38054e-16   	;Boltzmann's constant(ergs/K)
m=mass*1.6713e-24      	;weight of specific atom(g)
c=2.997925e10	   	;speed of light (cm/sec)
au=1.495979e13 		;astronomical unit (cm)
sun_sat_dist=0.99*au	;Distance between sun and satellite (cm)

;......Rename arrays and Multiply flux by 4pi(sun-sat distance)^2.........

IF keyword_set(specs) THEN BEGIN
	spec_first=specs(0)
	spec_last=specs(1)
ENDIF ELSE BEGIN
	spec_first=0
	spec_last=n_elements(flux_parameters.flux(0,*))-1
ENDELSE

flux=flux_parameters.flux(*,spec_first:spec_last)*4.*!pi*sun_sat_dist^2
eflux=flux_parameters.eflux(*,spec_first:spec_last)

tempnum=n_elements(flux)
for j=0,tempnum-1 DO BEGIN
	IF (flux(j) LT 0.) THEN BEGIN
PRINT,'Flux must always be positive.  Suggest subtracting less background'
PRINT,'or narrowing the data analysis window.'
	GOTO,DONE
	ENDIF
	IF (eflux(j) EQ 0.) THEN BEGIN
PRINT,'Error in flux must always be non-zero.  Suggest narrowing analysis'
PRINT,'window.'
	GOTO,DONE
	ENDIF
ENDFOR
wave=flux_parameters.wave
nflux=N_ELEMENTS(flux(*,0))
specnum=n_elements(flux(0,*))

mag_flux=float(floor(alog10(max(flux/(4.*!pi*sun_sat_dist^2)))))
							;find magnitude
norm_flux=10^(mag_flux)					;for plotting

;......................Set plot parameters............................
IF KEYWORD_SET(VERBOSE) THEN BEGIN
!p.thick=1.8
!p.charsize=1.3
!p.charthick=1.8
!p.region=0
!p.multi=0
!x.range=0
xlabel_flux='Wavelength ('+string(197b)+')'
ylabel_flux='10!u'+strcompress(string(fix(mag_flux)))+ $
	'!n photons s!u-1!n cm!u-2!n ' +string(197b) +'!u-1!n'
xlabel_vdem='Velocity (km s!u-1!n)'
ENDIF

;......................Look at flux?................

IF KEYWORD_SET(VERBOSE) THEN BEGIN

answer=''
READ,answer,PROMPT="Would you like to view the flux? (y/n): "

IF answer EQ 'y' THEN BEGIN
PRINT,''
PRINT,'Reformatting, please wait...'
PRINT,''
	window,/free,title='FLUX',xsize=450,ysize=450,xpos=0,ypos=400
	data=bytarr(450,450,specnum)
	!p.region=[.05,.05,.95,.95]
	max_flux=max(flux/(4.*!pi*sun_sat_dist^2))+max(eflux)

	FOR j=0,specnum-1 DO BEGIN

	tempflux=flux(*,j)/(norm_flux*4*!pi*sun_sat_dist^2)
	tempeflux=eflux(*,j)/norm_flux

	plot_err,wave,tempflux,yerr=tempeflux,psym=10,xtitle=xlabel_flux,$
		ytitle=ylabel_flux,yrange=[0,max_flux/norm_flux],$
		title='Spectrum:'+strcompress(string(fix(j+spec_first)))
	data(*,*,j)=tvrd()

	ENDFOR

stepper,data
data=0
wdelete
!p.region=0

answer=''
READ,answer,PROMPT="Would you like to continue? (y/n): "
IF answer EQ 'n' THEN GOTO,DONE

ENDIF

ENDIF


;........................Define velocity bins.............................


diff=wave(1)-wave(0)

IF KEYWORD_SET(DVEL) THEN dvel=dvel ELSE dvel=c*diff/L0

max= -1.*c*(min(wave)-L0)/L0  	;get max and min possible velocities
min= -1.*c*(max(wave) -L0)/L0
nvel=fix((max-min)/dvel)		;get # of needed velocity bins


velocity=dblarr(nvel)
velocity(0)=min+dvel/2.			;set velocities
for i=1,nvel-1 do velocity(i)=velocity(i-1) + dvel

;..........Create base kernel matrix for use in inversion................

basekernel=dblarr(nvel,nflux)
center=dblarr(nvel)

IF keyword_set(width_sumer) THEN width_sumer=width_sumer ELSE width_sumer=0.
v_instr=width_sumer*c/L0

center=L0*(1.-velocity/c)      ;wavelengths corresponding to velocity bin values

IF KEYWORD_SET(XI) THEN xi=xi ELSE xi=0.	;Set nonthermal velocity

sig=L0*sqrt(k*tbar/m + xi^2+v_instr^2)/c

IF KEYWORD_SET(CAXIX) THEN BEGIN	;Compute basekernel for CAXIX data,
					;including d13 satellite features
d13L0=3.182				;d13 feature's rest wavelength
d13center=dblarr(nvel)
d13center=d13L0*(1.-velocity/c)
d13sig=d13L0*sqrt(k*tbar/m + xi^2)/c
IF keyword_set(ref_file) then ref_file=ref_file ELSE BEGIN
	dir=getenv("VDEM")
	ref_file=dir+'excitcoef.dat'
ENDELSE
coef=rd_tfile(ref_file,/nocomment,/auto,/convert)
n19=0.87				;Ca XIX abundance
n20=0.05				;Ca XX abundance
wcoef=interpol(coef(1,*),coef(0,*),tbar)	;w line coefficient
d13coef=interpol(coef(2,*),coef(0,*),tbar)	;d13 coefficient
direcoef=interpol(coef(3,*),coef(0,*),tbar)	;dielectric recombination coef.
d13scale=(n19*d13coef)/(n19*wcoef+ n20*direcoef);scale factor for d13 in kernel

FOR j=0,nflux-1 DO BEGIN
FOR i=0,nvel-1 DO BEGIN

basekernel(i,j)=dvel*1./(sqrt(2.*!pi)*sig)* $
        exp(-(wave(j)-center(i))^2/(2.*sig^2)) + $
        dvel*d13scale/(sqrt(2.*!pi)*d13sig)* $
        exp(-(wave(j)-d13center(i))^2/(2.*d13sig^2))

	IF basekernel(i,j) LT 1.e-15 THEN  basekernel(i,j)=0.

ENDFOR
ENDFOR

ENDIF ELSE BEGIN

FOR j=0,nflux-1 DO BEGIN
FOR i=0,nvel-1 DO BEGIN

	basekernel(i,j)=dvel*1./(sqrt(2.*!pi)*sig)* $
		exp(-(wave(j)-center(i))^2/(2.*sig^2))

	IF basekernel(i,j) LT 1.e-15 THEN $
		basekernel(i,j)=0.

ENDFOR
ENDFOR

ENDELSE

;.....Compute the third difference matrix, Ke3, and the regularization.......
;........................matrix, baseH=Ke3^T Ke3...........................

ke3=dblarr(nvel,nvel-3)
baseh=dblarr(nvel,nvel)

FOR n=0,nvel-4 DO BEGIN

	p=n
	ke3(p,n)=-1.
	ke3(p+1,n)=3.
	ke3(p+2,n)=-3.
	ke3(p+3,n)=1.

ENDFOR

baseh=ke3 # transpose(ke3)
h=baseh*1.e16*(1./6.)/1.e4 	;Incorporate constants of smoothing parameter
				;into H and scale to make magnitude of
				;smoothing parameter more convenient for user
				;input
TrH=0.
TrH=h(0,0) & for i=0,nvel-2 do TrH=TrH+h(i+1,i+1)


;*********** Invert spectrum to Compute VDEM *****************************

kernel=dblarr(nvel,nflux)
invertvdem=dblarr(nvel,specnum)
evdem=dblarr(nvel,specnum)
tempretflux=dblarr(nflux,specnum)

for time=0,specnum-1 do begin

;........Normalize data and kernel by data's errors and Rename arrays........
        for k=0,nflux-1 do $
		kernel(*,k)=basekernel(*,k)/eflux(k,time)
        A=kernel # transpose(kernel)
        inten=flux(*,time)/eflux(*,time)

;........ Compute initial guess for smoothing parameter....................
        TrA=A(0,0) & for i=0,nvel-2 do TrA=TrA+A(i+1,i+1)
        initial=TrA/TrH         ;initial guess for smoothing parameter

        data=inten # transpose(kernel)          ;data=K^T g

	IF KEYWORD_SET(WEIGHT) THEN q=weight*initial $
	ELSE BEGIN
					;call gcvcompute sub-routine to obtain
                                        ;optimal q (optq)
                tempflux=flux(*,time)
                q1=initial
                q2=initial+0.05

               optq=gcvcompute(q1,q2,tempflux,basekernel,A,H,data,kernel,$
			num=nvel)
                if (optq lt 0.) then optq=abs(optq)

                if (optq lt initial) then q=initial else q=optq

        ENDELSE

;........Compute VDEM using LU decomposition routines to solve.............

                matrix0=A + q*H                 ;matrix=K^T K + qH

if (fix(strmid(!version.release,0,1)) lt 5) then begin
        lu_matrix=transpose(matrix0)          ;nr_ludcmp requires column format
        lu_data=transpose(data)              ;nr_lubksb requires column format
        nr_ludcmp,lu_matrix,permut,/double    ;matrix is replaced by its
					      ;LU decomposition
        invertvdem(*,time)=nr_lubksb(lu_matrix,permut,lu_data)

endif else begin

	lu_matrix0=matrix0
	lu_data=transpose(data)
	invertvdem(*,time)=v5_invert(lu_matrix0,lu_data)
endelse

PRINT,'VDEM calculated for spectrum:'+strcompress(string(fix(time+spec_first)))

;...............Compute errors in VDEM.................................

        inverse=invert(matrix0,status)
        if (status ne 0) then print,'Inverse of matrix0 is suspect'
        covariance=inverse#A#inverse
        for i=0,n_elements(covariance(*,0))-1 do $
                evdem(i,time)=4.*!pi*sun_sat_dist^2*sqrt(covariance(i,i))

;...................Check for negative values........................

negative=''
FOR j=0,nvel-1 DO BEGIN
IF invertvdem(j,time) LT 0. THEN negative='y'
ENDFOR

IF negative EQ 'y' THEN BEGIN
PRINT,''
PRINT,$
'The inverted vdem has values which are negative for spectrum number'$
	+strcompress(string(fix(time+spec_first)))
PRINT,'Its VDEM should not be trusted.'
PRINT,'***If negative values are at the endpoints, suggest subtracting '
PRINT,'   less background or narrowing data analysis window.'
PRINT,'***If negative values are within VDEM and user has supplied'
PRINT,'   a smoothing parameter, suggest using a smaller smoothing '
PRINT,'   parameter value.'
PRINT,''
ENDIF

ENDFOR

;........................View VDEM?.................................

IF KEYWORD_SET(VERBOSE) THEN BEGIN

answer=''
PRINT,'Would you like to view the observed flux with the VDEM that was'
PRINT,'deconvolved from it?'
READ,answer, PROMPT="(y/n):"

IF answer EQ 'y' THEN BEGIN

window,/free,title='FLUX AND DECONVOLVED VDEM',xsize=450,ysize=450,$
		ypos=400,xpos=0
data=bytarr(450,450,specnum)
!p.multi=[0,0,2,0,0]
PRINT,''
PRINT,'Reformatting, please wait...'
PRINT,''

max_flux=max(flux/(4.*!pi*sun_sat_dist^2))+max(eflux)

diff_high=max(wave)-L0
diff_low=L0-min(wave)
diff_max=max([diff_high,diff_low])

min_vdem=min(invertvdem)
max_vdem=max(invertvdem)+max(evdem)
mag_vdem=float(floor(alog10(max_vdem))-1.)
norm_vdem=10^(mag_vdem)
ylabel_vdem='10!u'+strcompress(string(fix(mag_vdem+5)))+$
	'!n photons s!u-1!n (km s!u-1!n)!u-1!n'

FOR time=0,specnum-1 DO BEGIN

tempflux=flux(*,time)/(4.*!pi*sun_sat_dist^2*norm_flux)
tempeflux=eflux(*,time)/norm_flux


!p.region=[0.05,0.525,0.95,0.95]

	plot_err,wave,tempflux,yerr=tempeflux,psym=10,$
		xtitle=xlabel_flux, ytitle=ylabel_flux,$
		title='Flux, Spectrum:'+$
		strcompress(string(fix(time+spec_first))),$
		yrange=[0,max_flux/norm_flux],$
		xrange=[(L0-diff_max),(L0+diff_max)]

tempvel=velocity/1.e5
max_vel=max(abs(tempvel))
tempvdem=invertvdem(*,time)/norm_vdem
tempevdem=evdem(*,time)/norm_vdem

!p.region=[0.05,0.05,0.95,0.475]

	plot_err,tempvel,tempvdem,yerr=tempevdem,psym=10,$
		xtitle=xlabel_vdem,ytitle=ylabel_vdem,$
		title='Deconvolved VDEM',$
		yrange=[min_vdem/norm_vdem,max_vdem/norm_vdem],$
		xrange=[max_vel,-max_vel]
data(*,*,time)=tvrd()

ENDFOR
stepper,data
data=0
wdelete
!p.region=0
!p.multi=0
ENDIF
ENDIF

;...........................Calculate convolved flux..................


FOR time=0,specnum-1 DO $
	tempretflux(*,time)=(invertvdem(*,time) # basekernel)/(4.*!pi*$
	sun_sat_dist^2)
retflux=tempretflux

;.........................View convolved flux?........................

IF KEYWORD_SET(VERBOSE) THEN BEGIN

answer=''
PRINT,'Would you like to compare the observed flux with the flux'
READ,answer, PROMPT="computed from VDEM? (y/n):"

IF answer EQ 'y' THEN BEGIN

window,/free,title='FLUX AND VDEM-CONVOLVED FLUX',xpos=0,ypos=400,$
	xsize=450,ysize=450
data=bytarr(450,450,specnum)
PRINT,''
PRINT,'Reformatting, please wait...'
PRINT,''

!p.region=[.05,.05,.95,.95]
max_flux=max(flux/(4.*!pi*sun_sat_dist^2))

FOR time=0,specnum-1 DO BEGIN

tempflux=flux(*,time)/(4.*!pi*sun_sat_dist^2*norm_flux)

temp_retflux=tempretflux/norm_flux


	plot,wave,tempflux,psym=10,$
		xtitle=xlabel_flux, ytitle=ylabel_flux,$
		yrange=[0,max_flux/norm_flux],$
		title='Spectrum:'+strcompress(string(fix(time+spec_first)))
	oplot,wave,temp_retflux(*,time),linestyle=2
data(*,*,time)=tvrd()
ENDFOR
stepper,data
data=0
wdelete
!p.region=0
ENDIF
ENDIF


vdem_parameters={vdem:invertvdem,evdem:evdem,velocity:velocity}

RETURN, vdem_parameters

DONE:

;...................Return plot parameters to normal..................

IF KEYWORD_SET(VERBOSE) THEN BEGIN
!p.thick=1
!p.charsize=1
!p.charthick=1
ENDIF

END


