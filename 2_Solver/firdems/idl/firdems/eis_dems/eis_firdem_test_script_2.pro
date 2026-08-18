; This script computes and plots the response of the firdem code to a set of log-normal
; distributions, with specified width and total em, at temperatures ranging from 10^5.5 to
; 10^7.0 Kelvin. This version uses synthetic EIS data
.compile ../estimate_quantile.pro
.compile ../dist_median.pro
.compile ../firdem.pro
.compile ../firdem_get_pixel.pro
.compile eis_firdem_wrapper.pro
.compile eis_firdem_test_2.pro
logtmin = 5.5
logtmax = 7.0
tmin = 10^logtmin
tmax = 10^logtmax
nt = 100
nclvl = 8

nx = 25
ny = 40
emtot = 5.0e28
width = 0.3

tempsin = tmin*(tmax/tmin)^(findgen(nt)/(nt-1))
demarr = eis_firdem_test_2(tempsin(0), tempsout, dem_in, tr_struct=tr_struct, nx=nx, ny=ny, emtot=emtot,sigma_logn=width*0.5,niter=niter,logtmin=logtmin-0.5,logtmax=logtmax+0.5)

itout_min = value_locate(tempsout,tmin)
itout_max = value_locate(tempsout,tmax)

nt_out = n_elements(tempsout)

demarr = dblarr(nt,nt_out)
chi2arr = dblarr(nt,nx,ny)
times = dblarr(nt)
for i=0,nt-1 do begin &$
	demarr(i,*) = eis_firdem_test_2(tempsin(i),tempsout,dem_in, tr_struct=tr_struct, chi2s=chi2s, nx=nx, ny=ny, times=time, emtot=emtot, sigma_logn=width*0.5,niter=niter, logtmin=logtmin-0.5, logtmax=logtmax+0.5) &$
	times(i) = time &$
	chi2arr(i,*,*) = chi2s &$
	print,i, nt &$
endfor

logtavgs = dblarr(nt)
logtout = alog10(tempsout)
for i=0,nt-1 do logtavgs(i) = dist_median(logtout,demarr(i,*)>0)

nt_out = n_elements(tempsout)

logtin = alog10(tempsin)
logtout = alog10(tempsout)

ia = lindgen(nt_out)
ishifts = lonarr(nt_out)
for i=0,nt-1 do ishifts(i) = ia(where(abs(tempsout-tempsin(i)) eq min(abs(tempsout-tempsin(i)))))

negdemarr = dblarr(nt,nt_out)
posdemarr = dblarr(nt,nt_out)
clrdemarr = dblarr(nt,nt_out,3)
ipos = where(demarr gt 0,count)
if count gt 0 then posdemarr(ipos) = demarr(ipos)
clrdemarr(*,*,0) = posdemarr
negs = where(demarr lt 0,count)
if count gt 0 then negdemarr(negs) = - demarr(negs)
clrdemarr(*,*,2) = negdemarr

xscale = alog10(tempsin(1))-alog10(tempsin(0))
yscale = alog10(tempsout(1))-alog10(tempsout(0))
x0 = alog10(tempsin(0))
y0 = alog10(tempsout(0))
yrange = [alog10(min(tempsout)),alog10(max(tempsout))]

yt = 'Extracted Log!D10!N(T) (Kelvin)'
xt = 'Injected Log!D10!N(T) (Kelvin)'

set_plot,'ps'
emtot_str=string(format='(%"%7.1e")',emtot)
width_str=string(format='(%"%g")',width)
device,file="plots/eis_gresp_width"+width_str+"_"+emtot_str+"EM.eps",/encapsulated
!p.font=0
device,/inches,xsize=5.6,ysize=1.6,/times
!p.multi=[0,3,1]

cmax = emtot/sqrt(4*acos(0)*(0.5*width)^2)
;cmax = max(smooth(max(posdemarr,dimension=2),round(0.25*nt)))

plot_image,posdemarr, min=0, max=cmax, scale=[xscale,yscale], origin=[x0,y0], xtitle = xt, ytitle = yt, title='EIS DEM Gaussian Response', yrange=yrange,/nosquare

clevels = min(demarr)+(cmax-min(demarr))*findgen(nclvl)/(nclvl-1)
contour, demarr,alog10(tempsin),alog10(tempsout),levels=clevels,/overplot,c_colors=64

oplot,alog10(tempsin),alog10(tempsin),linestyle=2,color=255,thick=2

oplot,logtin,logtavgs,linestyle=0,color=0,thick=6
oplot,logtin,logtavgs,linestyle=0,color=255,thick=2

legend,'Extracted EMWMT',linestyle=0,color=0,textcolors=255,charsize=0.45,box=0,thick=6
legend,'Extracted EMWMT',linestyle=0,color=255,textcolors=255,charsize=0.45,box=0,thick=2

title = 'EM=' + emtot_str +' cm!E-5!N' + ', Width='+width_str + ' Dex'
qlvls = [0.95,0.9,0.75,0.5]
qlvllabels = ['95th percentile','90th percentile','75th percentile','50th percentile']
nqlvls = n_elements(qlvls)
qlvlcurves = dblarr(nt, nqlvls)

for i=0,nt-1 do for j=0,nqlvls-1 do qlvlcurves(i,j) = estimate_quantile(chi2arr(i,*,*),qlvls(j))
plot,alog10(tempsin),qlvlcurves(*,0),linestyle=0,title = title,xtitle=xt,ytitle="Chi squared",yrange=[0,20]

for i=0,nqlvls-1 do oplot,alog10(tempsin),qlvlcurves(*,i),linestyle=i
legend,qlvllabels,linestyle=indgen(nqlvls),charsize=0.55,box=0

plot,alog10(tempsin),times/nx/ny,title="Average Time per DEM",xtitle=xt,ytitle="Time (Seconds)"

device,/close
set_plot,'x'
