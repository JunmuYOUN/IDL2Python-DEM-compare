4/23/98

This is the README.TXT document for a package of routines which computes
the Velocity Differential Emission Measure (VDEM) from spectral line
profiles observed by the Yohkoh/BCS or SoHO/SUMER instruments.  Please
direct any questions about these routines to Dr. Elizabeth Newton
(elizabeth.newton@msfc.nasa.gov).

The Velocity Differential Emission Measure (VDEM) can be inverted or
deconvolved from observed spectral lines to recover the line-of-sight
velocity distribution of the emitting ion.   VDEM has units of
(photons s-1 (km s-1)-1).  The moments of VDEM are related to the
plasma's mass, momentum, and enthalpy flux.   Details on the derivation
and various uses of VDEM may be found in:

Newton, Emslie, Mariska 1995, ApJ, 447,915
Newton, Emslie, Mariska 1996, ApJ, 459, 804
Newton 1997, ApJ, 484, 455

There are 6 IDL .pro routines and one ASCII datafile associated with the
VDEM package:  vdem.pro, prep_vdem_bcs.pro, prep_vdem_sumer.pro,
gcvcompute.pro, gcvfunction.pro, v5_invert.pro, and excitcoef.dat.  
In order to utilize the data file, the user must set a vdem 
environment variable by adding the following line to the 
idl_startup.pro:
		set_env, 'VDEM=/where_vdem_files_are/".

A brief description of each routine follows below.  More detailed
documentation also appears at the beginning of each routine, in the
standard SolarSoft format.  These routines permit the user to input
observed spectra and obtain VDEMs for the line profiles.  No
post-processing (mass, momentum, or enthalpy) calculations are provided
in these routines however.

It should also be noted that the VDEM routines have really only
been validated for the Ca XIX and SXV lines observed by the BCS
or for SUMER lines.  The Fe lines observed by BCS present a particular
challenge for VDEM analysis since so many satellite lines are evidently
blended into the resonance line and there is fairly low resolution.
For this reason, it is not recommended to employ the VDEM package for
analysis of the Fe  lines.

Routines:

(1)  PREP_VDEM_BCS:
This routine prepares BCS data for use in VDEM.PRO.  In this routine,
spectra are extracted from Yohkoh BCS spectrally calibrated data files.
User-provided background is subtracted, and the user may define the
wavelength range over which the line will be inverted.  The background
model and wavelength window are important for successfully computing VDEM.

Prior to running this routine, the user must have generated Yohkoh/BCS
index and data arrays using yodat.pro (to select the flare, times,
channel, integration interval, etc.) as well as mk_bsc.pro (to perform
instrument calibrations and corrections).  In order to determine the
other inputs required by PREP_VDEM_BCS.PRO, the user should examine the
extracted spectra using plot_bsc.pro.   Specifically, the user must
determine the bin numbers of the first and last wavelength bins which
will define the wavelength range over which the line will be inverted.
Since the BCS wavelength is not absolutely calibrated, the user must
also determine the amount of correction to apply to the wavelength scale
to account for the pointing offset which arises from the flare's location
on the Sun.  Generally, to do this the user may examine the spectra
occurring late in the flare to determine the "rest wavelength" appropriate
for the flare and subtract this from the ion's laboratory rest wavelength
to determine the offset correction.   Finally, the user must have also
determined the background (constant or vector) to subtract from the spectra,
utilizing whatever method one prefers.

Once PREP_VDEM_BCS.PRO has been run, the data is ready for input into VDEM.PRO.


(2)  PREP_VDEM_SUMER
This routine prepares SUMER data for use in VDEM.PRO.  User-defined pixels
are averaged and calibrated.  Then a background may be computed for the
spectra and subtracted.  Finally SUMER's specific intensity data are
converted to flux data for processing by VDEM.PRO.

Prior to running PREP_VDEM_SUMER, the user must have generated SUMER
index and data arrays using the rd_sumer.pro routine.  It is also suggested
that the user calculate the wavelength array for the data set utilizing
observations of reference lines.

Once PREP_VDEM_SUMER.PRO has been run, the data is ready for input into VDEM.PRO.


(3) VDEM.PRO
This routine is the main program for inverting spectral line profiles
(BCS or SUMER) to obtain the VDEM.  Input data must be prepared using
either PREP_VDEM_BCS or PREP_VDEM_SUMER.  The user must also input the
appropriate atomic information for the emitting ion, such as the mass,
laboratory rest wavelength, and the average temperature at which the
ion is formed.

As a default, VDEM.PRO performs the spectral line inversions utilizing the
GCV technique to calculate the optimal smoothing parameter.  However, the
user may provide this smoothing weight as an optional  input.  VDEM.PRO
also utilizes a default velocity bin size which is the instrument's
observational limit, i.e.,
dvel= speed of light * wavelength binwidth/rest wavelength.

For the deconvolution, VDEM.PRO constructs a gaussian kernel function
with a default width equal to the thermal Doppler velocity of the ion.
However, the user may provide a nonthermal velocity as an optional input
in order to account for any nonthermal broadening, and, in the case of
SUMER, the instrumental broadening, that may be consistently
present in the line.  In the special case of the CaXIX line profiles
observed by Yohkoh/BCS, for which the user sets the keyword /CAXIX,
VDEM.PRO constructs a kernel which contains not only the main resonance
line feature, but also a structure which accounts for emission by the
d13 satellite feature.  This satellite feature appears on the long
wavelength side of the CaXIX resonance line and results in persistent
line asymmetry.  Hence the VDEM inversion must account for this feature
in order to arrive at a velocity distribution that does not overestimate
plasma downflows.

This routine outputs a structure containing the VDEM, errors in the VDEM,
and the velocity array.  VDEM.PRO may optionally provide a flux computed
by forward convolving the VDEM, so that this flux may be compared to the
observed flux.

Users are cautioned that the deconvolved VDEM's behavior at its
endpoints is poorly constrained by the data, and hence its endpoint
values will exhibit artifacts which typically arise in inverse problems.


(4) GCVFUNCTION.PRO
This function computes the Generalized Cross Validation (GCV) function
which is minimized by GCVCOMPUTE.PRO in order to find the optimal smoothing
parameter for the VDEM inversion.  Details on the derivation of GCV may
be found in Golub, Heath, & Wahba.  Newton et al. 1996 may be consulted
for its application to VDEM.


(5) GCVCOMPUTE.PRO
This routine minimizes the GCV function to compute the optimal smoothing
parameter for the VDEM inversion.  GCVCOMPUTE.PRO utilizes the
single-dimension bracketing and Brent's minimization schemes which are
published in Numerical Recipes.


(6) v5_INVERT.PRO
This routine determines which inversion procedure is employed, 
depending on the user's IDL version.

(7) EXCITCOEF.DAT
This ASCII file contains a (4 columns x 11 rows) array which
has the excitation coefficients for the Ca XIX 'w' line, 'd13' satellite
line, and the dielectronic recombination features (which contribute to
emission in the same wavelength range as the 'w' line).  This file
is utilized only if Yohkoh/BCS CaXIX spectra are being inverted and the /CAXIX
switch is employed.  The data come from Bely-Dubau et al. Dec. 1982,
RAS, Monthly Notices, v. 201, p. 1155-1169, "Dielectronic Satellite
Spectra for highly-charged He-like ions. VII. Calcium Spectra."  The
first column comprises temperature; the second, 'w' line coefficients;
the third, 'd13' coefficients; and the fourth, dielectronic
recombination values.
