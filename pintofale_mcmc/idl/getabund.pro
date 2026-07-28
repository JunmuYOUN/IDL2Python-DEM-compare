function getabund,hint,elem=elem,source=source,norm=norm,metals=metals,$
	fipbias=fipbias,_extra=e
;+
;function	getabund
;		returns elemental abundances in an array
;
;syntax
;	abund=getabund(hint,elem=elem,source=source,norm=norm,$
;	metals=metals,fipbias=fipbias)
;
;parameters
;	hint	[INPUT] some hint as to what output is required.  examples:
;		'help' -- prints out usage
;		'anders', 'anders & grevesse', etc. -- Anders & Grevesse (1989)
;		'meyer', 'coronal', etc. -- Meyer (1985)
;		'allen', 'cosmic' -- Allen (1973)
;		'ross', 'aller', 'ross & aller', etc. -- Ross & Aller (1976)
;		'grevesse', 'grevesse et al.', etc. -- Grevesse et al (1992)
;			modifications to Anders & Grevesse (1989)
;		'feldman', etc. -- Feldman et al (1992)
;		'waljeski', etc. -- Waljeski et al. (1994)
;		'murphy', 'chromospheric', etc. -- Murphy (1985) chromospheric
;		'grevesse & anders' -- Grevesse & Anders (1991)
;		'grevesse & sauval' -- Grevesse & Sauval (1998) "standard
;			abundances"
;		'young', 'photo' -- photospheric abundances from Young et
;			al. (1997)
;		'fludra', 'schmelz' -- Fludra & Schmelz (1999) "hybrid coronal"
;		'path/file', etc. -- if string contains a "/" in it, then
;			read from disk file which contains one entry per
;			line, with atomic symbol first, and abundance next
;
;keywords
;	elem	[INPUT] element -- if given, return abundances for only
;		the specified element
;		* if string, converts to atomic number using SYMB2ZION
;		* if byte, integer, or float, assumed to be atomic number
;		* may be an array
;	source	[INPUT] another way to say "hint"
;		1: Anders & Grevesse (1989)
;		2: Meyer (1985)
;		3: Allen (1973)
;		4: Ross & Aller (1976)
;		5: Grevesse et al. (1992)
;		6: Feldman et al. (1992)
;		7: Waljeski et al. (1994)
;		8: Murphy (1985)
;		9: Grevesse & Anders (1991)
;		10: Grevesse & Sauval (1998)
;		11: Young et al. (1997)
;		12: Fludra & Schmelz (1999)
;		* overrides HINT iff HINT is not understood
;	norm	[INPUT] normalization
;		* use this to set output format
;		  NORM=1: output will be n(Z)/n(H) (this is the default)
;		  NORM>1: output will be NORM(0)+log10(n(Z)/n(H))
;		  	(e.g., NORM=12 is the standard log form)
;		  0<NORM<1: output will be NORM(0)*n(Z)/n(H)
;		  	(this can be useful when ELEM is set)
;		  -1<=NORM<0: output will be n(Z)/DEFAULT_n(Z)
;		  	(if you want to compare different abundance tables)
;		  NORM<-1: output will be log10(n(Z)/DEFAULT_n(Z))
;		  	(I have no idea why anyone would use this)
;	metals	[INPUT] metallicity.  if set, modifies the abundances of
;		all elements past He by the appropriate amount
;		* if not set, makes no changes
;		* assumed to be in log10 form if
;		  -- abs(NORM) > 1, or
;		  -- < 0
;		* assumed to be multiplicative factor otherwise
;		* NOTE: for example, if norm=12 then metals=0 means
;		  Solar abundances, while if norm=1 then metals=0
;		  means complete metal depletion!
;	fipbias	[INPUT] if set to a non-zero scalar, multiplies the
;		abundances of low-FIP elements by this factor.  be
;		careful not to use it for abundance lists that already
;		include the FIP effect.
;		* if not set, or set to 0, does nothing.
;	_extra	[INPUT] junk -- here only to prevent program from crashing
;
;subroutines
;	RDABUND
;	SYZE
;	SYMB2ZION
;
;history
;	vinay kashyap (Dec96)
;	avoid analyzing HINT if SOURCE is set (VK; Jan97)
;	added 'angr' and tightened 'Anders&Grevesse' v/s 'Grevesse' (VK; Jul97)
;	added option of reading from disk file (VK; Jun98)
;	added Murphy (1985) data, and allowed ELEM to be array (VK; Jul98)
;	corrected bug associated with non-array ELEMs (VK; Sep98)
;	added keyword METALS, deleted keywords ALTER, FORM, removed call to
;	  SETABUND, added call to RDABUND (VK; Dec98)
;	added Grevesse & Anders (1991), Grevesse & Sauval (1998), Young et
;	  al (1998) (VK; Apr99)
;	added keyword FIPBIAS (VK; Jul01)
;	added Fludra & Schmelz (1999); now HINT overrides SOURCE and program
;	  warns if there is a conflict (VK; Jun02)
;-

;	check input
np=n_elements(hint)
if np eq 0 or np gt 1 then todo='help' else todo=strlowcase(string(hint))
if not keyword_set(norm) then norm=1
;
if n_elements(metals) gt 0 then begin
  met=metals(0)
  if abs(norm) gt 1 or met lt 0 then lmet=1 else lmet=0
endif

codemax=13				;so far this many have been coded.
					;CODEMAX would be file input

code_s=0 & code_h=0

if keyword_set(source) then begin	;(SOURCE has been set
  szs=size(source) & nszs=n_elements(szs)
  code_s=0
  if szs(nszs-2) lt 7 then code_s=fix(source(0))
  if szs(nszs-2) eq 7 then begin
    code_s=codemax	;must be a filename
    hint=strtrim(source(0),2)	;the name of the file
  endif
endif					;SOURCE)

;	now analyze the HINT
  code_h=0					;default
  			;ANDERS & GREVESSE?
  i0=strpos(todo,'ander') & i1=strpos(todo,'gre')
  if i0 ge 0 or i1 ge 0 then code_h=1
  i0=strpos(todo,'angr') & if i0 ge 0 then code_h=1
  			;MEYER?
  i0=strpos(todo,'mey') & i1=strpos(todo,'mayer')
  i2=strpos(todo,'coro')		;coronal?
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=2
  			;ALLEN?
  i0=strpos(todo,'allen') & i1=strpos(todo,'allan')
  i2=strpos(todo,'cos')			;cosmic?
  if i0 ge 0 or i1 ge 0 or i2 ge 0 then code_h=3
  			;ROSS & ALLER?
  i0=strpos(todo,'ros') & i1=strpos(todo,'aller')
  if i0 ge 0 or i1 ge 0 then code_h=4
  			;GREVESSE, NOELS, SAUVAL?
  i0=strpos(todo,'gre') & i1=strpos(todo,'et. al') & i2=strpos(todo,'et al')
  i3=strpos(todo,'noe') & i4=strpos(todo,'sau')
  if i0 ge 0 and (i1 ge 0 or i2 ge 0) then code_h=5
  i5=strpos(todo,'ander')
  if i0 ge 0 and i5 ge 0 then code_h=1		;<-- Anders & Grevesse
  if i0 ge 0 and i5 lt 0 then code_h=5		;<-- no Anders w. Grevesse!
  if i5 gt i0 and i0 ge 0 then code_h=9		;<-- Grevesse & Anders
  if i4 ge 0 then begin
    if i0 ge 0 and i3 lt 0 then code_h=10	;<-- G&S, no N
    if i0 ge 0 and i3 ge 0 then code_h=5	;<-- G, N *and* S!
  endif
  			;FELDMAN, MANDELBAUM, SEELY, DOSCHEK, GURSKY?
  i0=strpos(todo,'fel') & i3=strpos(todo,'mand')
  i4=strpos(todo,'see') & i5=strpos(todo,'dos') & i6=strpos(todo,'gur')
  if i0 ge 0 and (i1 ge 0 or i2 ge 0) then code_h=6
  if i0 ge 0 or i3 ge 0 or i4 ge 0 or i5 ge 0 or i6 ge 0 then code_h=6
  			;WALJESKI, MOSES, DERE, SABA, STRONG, WEBB, ZARRO?
  i0=strpos(todo,'wal') & i3=strpos(todo,'mos') & i4=strpos(todo,'dere')
  i5=strpos(todo,'sab') & i6=strpos(todo,'stro') & i7=strpos(todo,'we')
  i8=strpos(todo,'zar')
  if i0 ge 0 and (i1 ge 0 or i2 ge 0) then code_h=7
  if i0 ge 0 or i3 ge 0 or i4 ge 0 or i5 ge 0 or i6 ge 0 $
	or i7 ge 0 or i8 ge 0 then code_h=7
			;MURPHY (1985)?
  i0=strpos(todo,'mur') & i1=strpos(todo,'chromo')
  if i0 ge 0 or i0 ge 0 then code_h=8
			;GREVESSE & ANDERS (1991)?
  i0=strpos(todo,'gre') & i1=strpos(todo,'ander')
  if i0 ge 0 and i1 ge 0 and i1 gt i0 then code_h=9
			;GREVESSE & SAUVAL (1998)?
  i0=strpos(todo,'gre') & i1=strpos(todo,'sau')
  if code_h ne 5 and i0 ge 0 and i1 ge 0 then code_h=10
			;YOUNG ET AL (1997)/PHOTOSPHERIC?
  i0=strpos(todo,'you') & i1=strpos(todo,'phot')
  if i0 ge 0 or i1 ge 0 then code_h=11
			;FLUDRA & SCHMELZ (1999)
  i0=strpos(todo,'flu') & i1=strpos(todo,'sch')
  if i0 ge 0 or i1 ge 0 then code_h=12
			;disk file?
  i0=strpos(todo,'/')
  if i0 ge 0 then code_h=codemax
  ;
  			;AIUTO!
  i0=strpos(todo,'help') & if i0 ge 0 then code_h=0

;	the above block used to be inside the first if block below,
;	until the possibility of a SOURCE/HINT conflict was addressed
;if not keyword_set(source) then begin			;{analyze the hint
;endif else begin				;}{hint, shint..
;  szs=size(source) & nszs=n_elements(szs)
;  code=0
;  if szs(nszs-2) lt 7 then code=fix(source(0))
;  if szs(nszs-2) eq 7 then begin
;    code=codemax	;must be a filename
;    hint=strtrim(source(0),2)	;the name of the file
;  endif
;endelse						;loud and clear}

;	catch some trivial errors

;	what if there is a conflict between SOURCE and HINT?
if not keyword_set(source) then code_s = code_h else $
 if np ne 1 then code_h = code_s
if code_s ne code_h then begin
  message,'keyword value incompatible with parameter value',/info
  if code_h ne 0 then begin
    message,'choosing to go with HINT="'+todo+'"',/info
    code=code_h
  endif else begin
    message,'choosing to go with the SOURCE='+strtrim(code_s,2),/info
    code=code_s
  endelse
endif else code=code_h
;
if code lt 0 then code=0 & if code gt codemax then code=0

;	usage
if code eq 0 then begin
  print,'Usage: abund=getabund(hint,elem=elem,source=source,norm=norm,$'
  print,'       metals=metals,fipbias=fipbias)'
  print,'  returns elemental abundances'
  print,''
  print,'  SOURCE=1: Anders & Grevesse (1989) <<<DEFAULT>>>'
  print,'  SOURCE=2: Meyer (1985)'
  print,'  SOURCE=3: Allen (1973)'
  print,'  SOURCE=4: Ross & Aller (1976)'
  print,'  SOURCE=5: Grevesse et al. (1992)'
  print,'  SOURCE=6: Feldman et al. (1992)'
  print,'  SOURCE=7: Waljeski et al. (1994)'
  print,'  SOURCE=8: Murphy (1985)'
  print,'  SOURCE=9: Grevesse & Anders (1991)'
  print,'  SOURCE=10: Grevesse & Sauval (1998)'
  print,'  SOURCE=11: Young et al. (1997)'
  print,'  SOURCE=12: Fludra & Schmelz (1999)'
  print,''
  print,'  0<NORM<=1:	ABUND(Z)=NORM*n(Z)/n(H) {DEFAULT: NORM=1}'
  print,'  NORM>1:	ABUND(Z)=log10(n(Z)/n(H))+NORM'
  print,'  -1<=NORM<0:	ABUND(Z)=ABUND(Z)/DEFAULT_ABUND(Z)'
  print,'  NORM<-1:	ABUND(Z)=log10(ABUND(Z)/DEFAULT_ABUND(Z))'
  print,''
endif

;	default abundances
defabu=[1., 0.0977, 1.45e-11, 1.15e-11, 3.98e-8, 3.63e-4, 1.12e-4,$
	8.51e-4, 3.63e-8, 1.23e-4, 2.14e-6, 3.8e-5, 2.95e-6, 3.55e-5,$
	2.82e-7, 1.62e-5, 3.16e-7, 3.63e-6, 1.32e-7, 2.29e-6, 1.26e-9,$
	9.77e-8, 1e-8, 4.68e-7, 2.45e-7, 4.68e-5, 8.32e-8, 1.78e-6,$
	1.62e-8, 3.98e-8] & nelem=n_elements(defabu)
;Anders, E. \& Grevesse, N.\ 1989, {\it Geochimica et Cosmochimica Acta},
;53, 197

;	database of abundances
case code of

  2: begin	;Meyer, J.P.\ 1985, ApJS, 57, 173
		;NOTE: CHIANTI lists abund(7,10)=[8.39,6.44]
    abund=fltarr(nelem)
    abund(0)=12.00 & abund(1)=10.99 & abund(5)=8.37 & abund(6)=7.59
    abund(7)=8.33 & abund(9)=7.55 & abund(11)=7.57 & abund(12)=6.44
    abund(13)=7.59 & abund(15)=6.94 & abund(17)=6.33 & abund(19)=6.47
    abund(25)=7.59 & abund(27)=6.33
    oo=where(abund ne 0) & abund(oo)=10.^(abund(oo)-abund(0))
    oo=where(abund eq 0) & abund(oo)=defabu(oo)
  end

  3: begin	;Allen, C.W.\ 1973, {it Astrophysical Quantities}, 3$^{rd}$
		;Ed., Athlone Press, London
    frac=[1., 0.87, 0.35, 0.89, 2.51, 0.91, 0.81, 0.78, 1.10, 0.68, 0.83,$
	0.69, 0.83, 0.93, 1.17, 0.98, 1.26, 1.74, 0.68, 0.87, 1.32, 1.38,$
	2.51, 1.51, 1.02, 0.85, 1.51, 1.12, 1.95, 0.4]
    abund=frac*defabu
  end

  4: begin	;Ross, J.E., \& Aller, L.H.\ 1976, {\it Science}, 191, 1223
    frac=[1., 0.65, 0.69, 1., 0.32, 1.15, 0.78, 0.81, 1., 0.3, 0.89, 1.02,$
	  1.12, 1.26, 1.12, 0.98, 1., 0.28, 1.1, 0.98, 0.87, 1.15, 1.05,$
	  1.1, 1.07, 0.68, 0.95, 1.07, 0.71, 0.71]
    abund=frac*defabu
  end

  5: begin	;Grevesse, N., Noels, A., \& Sauval, A.J.\ 1992, in {\it
		;Coronal Streamers, Coronal Loops and Coronal and Solar
		;Wind Composition}, Proc.\ 1st SOHO Workshop, ESA SP-348, ESA
		;Publ.\ Div., ESTEC, Noordwijk, p.305
    abund=defabu
    abund(1)=0.095 & abund(5)=3.55e-4 & abund(6)=9.33e-5 & abund(7)=7.41e-4
    abund(8)=3.63e-8 & abund(9)=1.2e-4 & abund(16)=3.16e-7 & abund(17)=3.31e-6
    abund(20)=1.58e-9 & abund(21)=1.1e-7 & abund(25)=3.24e-5
  end

  6: begin	;Feldman, U., Mandelbaum, P., Seely, J.L., Doschek, G.A., \&
		;Gursky H.\ 1992, ApJS, 81, 387
    abund=[-3.41, -4.00, -3.11, -7.4, -3.92, -5.07, -3.85, -4.96,$
	-3.90, -6.48, -4.73, -6.4, -5.42, -7.05, -5.07, -8.78, -6.87, -7.6,$
	-6.15, -6.60, -3.90, -6.9, -5.16, -7.5, -7.8]
    abund=10.^(abund) & abund=[1.,10.^(-1.1),defabu(2:4),abund]
  end

  7: begin	;Waljeski, K., Moses, D., Dere, K.P., Saba, J.L., Strong, K.T.,
		;Webb, D.F., \& Zarro, D.M.\ 1994, ApJ, 429, 909
    abund=fltarr(nelem)
    abund(0)=12. & abund(1)=10.90 & abund(5)=9.28 & abund(6)=8.50
    abund(7)=9.30 & abund(9)=8.45 & abund(11)=8.48 & abund(12)=7.35
    abund(13)=8.50 & abund(15)=7.84 & abund(17)=7.24 & abund(19)=7.38
    abund(25)=8.50 & abund(27)=7.24
    oo=where(abund ne 0) & abund(oo)=10.^(abund(oo)-abund(0))
    oo=where(abund eq 0) & abund(oo)=defabu(oo)
  end

  8: begin	;Murphy, R.J.\ 1985, Ph.D.\ Thesis, Univ. of Maryland
    abund=fltarr(nelem)
    abund(1-1)=1. & abund(2-1)=7.57e-2 & abund(6-1)=1.14e-4
    abund(7-1)=4.82e-5 & abund(8-1)=1.72e-4 & abund(10-1)=7.84e-5
    abund(12-1)=2.63e-5 & abund(14-1)=3.92e-5 & abund(16-1)=1.88e-5
    abund(20-1)=6.27e-6 & abund(26-1)=2.47e-5
    oo=where(abund eq 0) & abund(oo)=defabu(oo)
  end

  9: begin	;Grevesse, N. \& Anders, E.\ 1991, in Solar Interior and
		;Atmosphere, ed. A.N.Cox, W.C.Livingston, & M.S.Matthews,
		;(Tucson: Univ. Arizona Press), 1227.
    abund=fltarr(nelem)
    abund(1-1)=12. & abund(2-1)=10.90
    abund(6-1:30-1)=[8.60,8.00,8.93,4.56,8.09,6.33,7.58,6.47,7.55,5.45,$
	7.21,5.5,6.56,5.12,6.36,3.10,4.99,4.00,5.67,5.39,7.67,4.92,$
	6.25,4.21,4.60]
    oo=where(abund eq 0,moo)
    abund=10.^(abund-12.) & if moo gt 0 then abund(oo)=defabu(oo)
  end

  10: begin	;Grevesse, N. \& Sauval, A.J.\ 1998, in Space Science Reviews,
		;85, 161-174
		;CHIANTI says they are same as Grevesse, Noels, & Sauval 1996.
    abund=[12.00,10.93,1.10,1.40,2.55,8.52,7.92,8.83,4.56,8.08,6.33,7.58,$
	6.47,7.55,5.45,7.33,5.50,6.40,5.12,6.36,3.17,5.02,4.00,5.67,5.39,$
	7.50,4.92,6.25,4.21,4.60]
    abund=10.^(abund-12.)
  end

  11: begin	;Young, P.R., Mason, H.E., Keenan, F.P., \& Widing K.G.
		;1997, A&A, 323, 243
		;also see CHIANTI's photospheric_may97.abund
    abund=fltarr(nelem)
    abund(1-1)=12. & abund(2-1)=10.99
    abund(6-1:30-1)=[8.55,7.97,8.87,4.56,8.08,6.33,7.58,6.47,7.55,5.45,$
	7.21,5.5,6.47,5.12,6.36,3.20,5.02,4.00,5.67,5.39,7.51,4.92,6.25,$
	4.21,4.60]
    oo=where(abund eq 0,moo)
    abund=10.^(abund-12.) & if moo gt 0 then abund(oo)=defabu(oo)
  end

  12: begin	;Fludra, A. & Schmelz, J. T., 1999, A&A, 348, 286
		;CHIANTI's fludra_schmelz_hybrid.abund
    abund=fltarr(nelem)
    iz=[1,    2,    6,   7,   8,   10,  11,  12,  13,  14,  15,  16,$
	17,  18,  19,  20,  22,  24,  26,  28,  30]
    xz=[12.00,10.80,8.41,7.81,8.74,7.95,6.63,7.90,6.80,7.87,5.44,7.32,$
	5.08,6.36,5.46,6.66,5.25,6.00,7.83,6.56,4.09]
    abund(iz)=xz
    oo=where(abund eq 0,moo)
    abund=10.^(abund-12.) & if moo gt 0 then abund(oo)=defabu(oo)
  end

  13: begin	;read from disk file
    abund=rdabund(hint,norm=1, _extra=e)
  end

  else: abund=defabu

endcase

;	size check
nab=n_elements(abund)
if nab lt nelem then abund=[abund,defabu(nab:*)]

;	include metallicity
if n_elements(metals) gt 0 then begin
  if keyword_set(lmet) then met=10.^(met)
  abund(2:*)=met*abund(2:*)
endif

;	include FIP bias
if keyword_set(fipbias) then begin
  fip=abs(float(fipbias[0]))
  ;	the following are elements with FIP < 10 eV
  zlofip=[3,4,5,11,12,13,14,19,20,21,22,23,23,25,26,27,28,29,30]-1L
  if fip gt 0 then abund[zlofip]=abund[zlofip]*fip else message,$
	'FIPBIAS of '+strtrim(fip,2)+' not understood; ignoring',/info
endif

;	extract element(s)
sze=size(elem) & nsze=n_elements(sze) & icode=sze(nsze-2) & nn=sze(nsze-1)
case icode of
  1: atno=((elem > 1) < nelem)
  2: atno=((elem > 1) < nelem)
  3: atno=((elem > 1) < nelem)
  4: atno=((elem > 1) < nelem)
  7: begin
     atno=intarr(nn)
     for i=0,nn-1 do begin
       symb2zion,elem(i),x,y & atno(i)=x
     endfor
  end
  else: atno=lindgen(nelem)+1
endcase
;
if icode gt 0 then atno=atno(uniq(atno,sort(atno)))
abund=[abund(atno-1)]

;	renorm
if norm(0) gt 0 and norm(0) le 1 then abund=norm(0)*abund
if norm(0) gt 1 then abund=alog10(abund)+norm(0)
if norm(0) lt 0 and norm(0) ge -1 then abund=abund/defabu
if norm(0) lt -1 then abund=alog10(abund/defabu)

return,abund
end
