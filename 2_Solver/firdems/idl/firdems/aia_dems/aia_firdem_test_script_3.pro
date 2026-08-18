; This script demonstrates how to create synthetic AIA data from a random test DEM, reconstruct
; the DEM using the firdem code, and plot the results. Repeatedly running the script can
; provide invaluable insight into how well input DEMs can generally be recovered.
.compile ../estimate_quantile.pro
.compile ../firdem_get_pixel.pro
.compile ../firdem.pro
.compile aia_firdem_wrapper.pro

; Generate AIA Temperature Response functions:
 if(n_elements(tr) eq 0) then tr=aia_get_response(/temp,/dn,/evenorm)

npeaks=5
logt = tr.logte
nt = n_elements(tr.logte)
pi = 2*acos(0)
; Generate some Gaussian peaks:
tcen_min = 6.0 ; Minimum log_10(T) for Gaussian peaks
tcen_max = 7.0 ; Maximum log_10(T) for Gaussian peaks
tsig_min = 0.025 ; Minimum standard deviation (in log_10(T)) for Gaussian peaks
tsig_max = 0.150 ; Maximum standard deviation (in log_10(T)) for Gaussian peaks
etot_min = 0.5e28 ; Minimum total emission for Gaussian peaks
etot_max = 1.5e28 ; Maximum total emission for Gaussian peaks
cens = tcen_min+(tcen_max-tcen_min)*randomu(seed,npeaks)
sigs = tsig_min+(tsig_max-tsig_min)*randomu(seed,npeaks)
amps = etot_min+(etot_max-etot_min)*randomu(seed,npeaks)
dem = dblarr(nt)
for i=0,npeaks-1 do dem += amps(i)*exp(-(logt-cens(i))^2/(2*sigs(i)^2))/sqrt(2*pi*sigs(i)^2)

; Produce AIA data numbers. The subroutine aia_firdem_synthesize_data in aia_firdem_wrapper.pro
; will do the same thing:
channels_good = where(tr.channels ne 'A304') ; This channel is unsuitable for DEM analysis
channels = tr.channels(channels_good)
nchan = n_elements(channels)
tresps = tr.all(*,channels_good)
dns0 = dblarr(1,1,nchan)
for i=0,nchan-1 do dns0(*,*,i) = int_tabulated(logt,dem*tresps(*,i))

G = aia_firdem_dnperphot(channels)
rdn = aia_firdem_readnoise(channels)
phots0 = dns0/G > 0; Convert to photons
nx=50 ; The total number of random noise realizations to produce is nx*ny.
ny=50
dns=dblarr(nx,ny,nchan)
; Generate data numbers with noise:
for i=0,nchan-1 do dns(*,*,i) = G(i)*randomu(seed,nx,ny,poisson=phots0(i))+rdn(i)*randomn(seed)

exptimes = 1.0+dblarr(nchan)
; Compute the DEM:
aia_firdem_wrapper, channels, exptimes, dns, dem_struc, tr_struct=tr, telapsed=t
; Plot The DEMs:
plot,logt,dem,title=string(t/nx/ny)+" Seconds per DEM"
for i=0,4 do begin &$
	for j=0,4 do begin &$
		firdem_get_pixel,dem_struc,i,j,t,dem_out &$
		oplot,alog10(t),dem_out,linestyle=1 &$
	endfor &$
endfor
legend,["Input","Output"], linestyle=[0,1]
