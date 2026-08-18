pro aia_teem_table,wave_,tsig_,te_range,teem_table
;+
; Project     : AIA/SDO
;
; Name        : AIA_TEEM_TABLE 
;
; Category    : AIA Temperature and Emission measure analysis   
;
; Explanation : calculates AIA fluxes for 6 wavelengths (WAVE_)
;		based on single-Gaussian DEM distributions  
;
; Syntax      : IDL>aia_dem_table,wave_,tsig_,te_range,teem_table
;
; Inputs      : wave_(newave)	AIA wavelength set [Angstroem]
;		tsig_(nsig)	Gaussian half width [0.01,...,1.00]
;		te_range(2)	temperature range [K]: [Tmin,Tmax,Tminpeak,Tmaxpeak]
;		teem_table	name of IDL savefile for output parameters:
;
; Outputs     : teem_table	creates IDL savefile containing:
; 			wave_(newave)	AIA wavelength set [Angstroem]
;			tsig_(nsig)	Gaussian half width [0.01,...,1.00]
;			resp(nte,nwave)	Response function
;			telog(nte)	logarithmic temperature range
;			dte(nte)	logarithmic temperature bins
;			flux(nte,nsig,nwave) AIA fluxes for unity emission measure
;
; History     : 23-Oct-2014, Version 1 written by Markus J. Aschwanden
;
; Contact     : aschwanden@lmsal.com
;-

;_____________________AIA RESPONSE FUNCTION________________________
tresp	=aia_get_response(/temp,/full,/dn)
telog_  =tresp.logte
telog1	=alog10(te_range(0))
telog2	=alog10(te_range(1))
ind_te	=where((telog_ ge telog1) and (telog_ le telog2),nte)
telog	=telog_(ind_te)
nwave	=n_elements(wave_)
ichan_  =fltarr(nwave)
resp    =fltarr(nte,nwave)
for iw=0,nwave-1 do begin
 filter ='A'+wave_(iw)
 if (wave_(iw) eq '094') or (wave_(iw) eq '94') then filter='A94'
 ichan  =where(tresp.channels eq filter)
 resp_  =tresp.tresp(*,ichan)
 resp(*,iw)=resp_(ind_te)
endfor

;_____________________CALCULATES LOOPUP TABLES_____________________
dte1	=10.^telog(1:nte-1)-10.^telog(0:nte-2)
dte	=[dte1,dte1[nte-2]]
em1	=1.			;unity emission measure
nsig	=n_elements(tsig_)
flux	=fltarr(nte,nsig,nwave)
for i=0,nte-1 do begin
 for j=0,nsig-1 do begin
  em_kelvin=em1*exp(-(telog-telog(i))^2/(2.*tsig_(j)^2))
  for iw=0,nwave-1 do flux(i,j,iw)=total(resp(*,iw)*em_kelvin*dte)
 endfor
endfor
print,'Lookup table created : ',teem_table
save,filename=teem_table,wave_,tsig_,resp,telog,dte,te_range,flux,dte
help,teem_table,wave_,tsig_,resp,telog,dte,te_range,flux,dte
end
