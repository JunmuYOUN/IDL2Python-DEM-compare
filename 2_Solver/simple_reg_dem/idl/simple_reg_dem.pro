;+
; function simple_reg_dem, data, errors, exptimes, logt, tresps, chi2, $
;		kmax=kmax, kcon=kcon, steps=steps, drv_con=drv_con, chi2_th=chi2_th, tol=tol
; 
; Computes a DEM given a set of input data, errors, exposure times, temperatures, and
; temperature response functions. Inputs:
;	data: image data (dimensions n_x by n_y by n_channels)
;	errors: image of uncertainties in data (dimensions n_x by n_y by n_channels)
;	exptimes: exposure times for each image (dimensions n_channels)
;	logt: Array of temperatures (dimensions n_temperatures)
;	tresps: Arrays of temperature response functions (dimensions n_temperatures by n_channels)
;
; Returns:
;	array of DEMs (dimensions n_x by n_y by n_temperatures). Dimensions depend on the units
;	of logt (need not be log temperature, despite the name) and tresps. For instance,
;	if tresps has units cm^5 and the logt values are base 10 log of temperature in Kelvin,
;	the output DEM will be in units of cm^{-5} per unit Log_{10}(T).
;
; Outputs:
;	chi2: Array of reduced chi squared values (dimensions n_x by n_y)
;
; Optional Keywords:
;	kmax: The maximum number of iteration steps (default - 100).
;	kcon: Initial number of steps before terminating if chi^2 never improves (default - 5).
;	steps: Two element array containing large and small step sizes (default - 0.1 and 0.5).
;	drv_con: The size of the derivative constraint - threshold delta_0 limiting the change in
;			 log DEM per unit input temperature (default - 8; e.g., per unit log_{10}(T)).
;	chi2_th: Reduced chi^2 threshold for termination (default - 1).
;	tol: How close to chi2_th the reduced chi^2 needs to be (default - 0.1).
;
; Author: Joseph Plowman -- 11-07-2019
; See: https://arxiv.org/abs/2006.06828
;-
function simple_reg_dem, data, errors, exptimes, logt, tresps, chi2, $
		kmax=kmax, kcon=kcon, steps=steps, drv_con=drv_con, chi2_th=chi2_th, tol=tol

	if(n_elements(kmax) ne 1) then kmax = 100
	if(n_elements(kcon) ne 1) then kcon = 5
	if(n_elements(steps) lt 2) then steps = [0.1,0.5]
	if(n_elements(drv_con) ne 1) then drv_con = 8.0
	if(n_elements(chi2_th) ne 1) then chi2_th = 1.0
	if(n_elements(tol) ne 1) then tol = 0.1

	nt = n_elements(logt)
	nx = n_elements(data[*,0,0])
	ny = n_elements(data[0,*,0])
	dT = logt[1:nt-1]-logt[0:nt-2]
	Bij = (diag_matrix([0,dT]+[dT,0])*2.0 + shift(diag_matrix([0,dT]),-1) + $
		shift(diag_matrix([dT,0]),1))/6.0
	Rij = transpose(tresps*((1+dblarr(nt))#exptimes))#Bij ; Matrix mapping coefficients to data
	Dij = diag_matrix([0,1/dT]+[1/dT,0]) - shift(diag_matrix([0,1/dT]),-1) - $
		shift(diag_matrix([1/dT,0]),1)
	regmat = Dij*n_elements(exptimes)/(drv_con^2*(logt[nt-1]-logt[0]))
	rvec = total(Rij,2)

	dems=fltarr(nx,ny,nt)
	chi2=fltarr(nx,ny)-1.
	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			err = reform(errors[i,j,*])
			dat0 = reform(data[i,j,*]) > 0.0
			s = alog(total((rvec)*((dat0 > 1.0e-2)/err^2))/total((rvec/err)^2)/(1+dblarr(nt)))
			for k=0,kmax-1 do begin
				dat = (dat0-Rij#((1-s)*exp(s)))/err ; Correct data by f(s)-s*f'(s)...
				mmat = Rij*((1.0/err)#exp(s)) ; Weight mapping by 1/err and f'(s)...
				amat = transpose(mmat)#mmat+regmat
				la_choldc,amat,status=stat
				if(stat eq 0) then begin
					c2p = mean((dat0-Rij#(exp(s)))^2.0/err^2)
					deltas = la_cholsol(amat,transpose(mmat)#dat)-s
					deltas *= (max(abs(deltas)) < 0.5/steps[0])/max(abs(deltas))
					ds = 1-2*(c2p lt chi2_th) ; Direction sign; is chi^2 too large or too small?
					c20 = mean((dat0-Rij#(exp(s+deltas*ds*steps[0])))^2.0/err^2)
					c21 = mean((dat0-Rij#(exp(s+deltas*ds*steps[1])))^2.0/err^2)
					interp_step = ((steps[0]*(c21-chi2_th)+steps[1]*(chi2_th-c20))/(c21-c20))
					s += deltas*ds*((interp_step > steps[0]) < steps[1])
					chi2[i,j] = mean((dat0-Rij#(exp(s)))^2.0/err^2)
				endif else break
				if((ds*(c2p-c20)/steps[0] lt tol)*(k gt kcon) or abs(chi2[i,j]-chi2_th) lt tol) then break
			endfor
			dems[i,j,*] = exp(s)
		endfor
	endfor
	return,dems

end
