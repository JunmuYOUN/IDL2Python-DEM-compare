FUNCTION PREP_VDEM_BCS,index,data,CHAN=chan,WBINS=wbins,$
	OFFSET=offset, BCKGRD=bckgrd

; NAME:
;	prep_vdem_bcs.pro

; PURPOSE:
;	To extract spectral data from Yohkoh BCS spectrally calibrated
;	data files, subtract background, and define appropriate wavelength
;	range. Part of the package of routines associated with computation
;	of Velocity Differential Emission Measure (VDEM) from spectral lines.

; CATEGORY:
;	Data analysis

; CALLING SEQUENCE:
;	flux_parameters=PREP_VDEM_BCS(index,data,CHAN=chan,WBINS=wbins,$
;			OFFSET=offset,BCKGRD=bckgrd)

; INPUTS:
;	index=index array associated with YOHKOH BCS data set (from mk_bsc.pro)
;	data=data array associated with YOHKOH BCS data set (from mk_bsc.pro)
;	chan=Yohkoh BCS channel number: (1=Fe XXVI, 2=Fe XXV, 3=Ca XIX, 4=S XV)
;	wbins= an array containing [wfirst,wlast], the wavelength bin numbers
;		which define the wavelength interval over which spectral line
;		will be examined and inverted. (Suggest that user decide
;		by inspecting spectra with plot_bsc.)
;	offset=the amount by which the flare spectra's 'rest wavelength'
;		differs from the ion's laboratory rest wavelength
;		(Angstroms).  This is necessary because BCS does not
;		have absolute wavelength calibration, and offset will
;		vary with location of flare on the Sun. (Suggest that
;		user inspect, with plot_bsc, last spectra observed
;		in flare to determine ion's "rest wavelength"
;		for given flare location.)
;	bckgrd=vector or constant value employed for background subtraction
;		(units of flux)

; OPTIONAL INPUTS:
;	None

; OUTPUTS:
;	flux_parameters=structure containing the flux, errors in flux, and
;		wavelength array which will be input to VDEM.pro.

; OPTIONAL OUTPUTS:
;	None

; COMMON BLOCKS:
;	None

; SIDE EFFECTS:
;	None

; RESTRICTIONS:
;	Yohkoh BCS index and data arrays, generated from yodat.pro and
;	mk_bsc.pro, must be the inputs.

; PROCEDURE:
;	(1) Use yodat.pro to extract index & data arrays for selected
;	flare, channel, integration interval, etc.  (2) Then call mk_bsc
;	routine to do instrument corrections and calibrations
;	(mk_bsc,index,data,index,data,dp_sync=dp_sync).
;	(3) Examine spectra using plot_bsc to determine values for:
;	wfirst,wlast,offset,& bckgrd.  Use these values as inputs to
;	this routine.

;	After PREP_VDEM_BCS.PRO has been run, user can run VDEM.PRO to invert
;	spectral lines.

; MODIFICATIONS:
;	E. Newton 1994: Original Code
;	E. Newton 1998: Modification for SolarSoft submission with VDEM package

;=====================================================================

;....Extract flux, wavelength and error arrays from BCS spectrally
;....calibrated data

sel_bsc,index,data,newindex,newdata,wave,flux,eflux,chan=chan


;.......Adjust wavelength array to bin center and account for BCS pointing
;.......offset

wfirst=wbins(0)
wlast=wbins(1)
width=fltarr(1,(wlast-wfirst)+1)      ;shift wavelength values to center of bin
for p=1,(wlast-wfirst)+1 do width(p-1)=wave(wfirst+p,0)-wave(wfirst-1+p,0)
lambda=wave(wfirst:wlast,0)
midbin=(lambda+width/2.) + offset   ;correct for BCS wavelength pointing offset


;...............Background subtraction........................................

flux2=(flux(wfirst:wlast,*) - bckgrd)    ;photons cm-2 s-1 A-1
eflux2=eflux(wfirst:wlast,*)		 ;error in photons cm-2 s-1 A-1


;.............................................................................


flux_parameters={flux:flux2, eflux:eflux2, wave:midbin}

return, flux_parameters

end

