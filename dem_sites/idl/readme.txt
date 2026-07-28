This package implements the SITES DEM inversion method described in [1]. The package also includes the Grid-SITES method, that can greatly increase the efficiency of SITES (or other DEM inversions) described in [2].

[1] https://ui.adsabs.harvard.edu/abs/2019SoPh..294..135M/abstract
[2] https://ui.adsabs.harvard.edu/abs/2019SoPh..294..136P/abstract

If you use these methods in your work, please cite the papers.

QUICK-START
Perhaps the quickest way of using these codes for AIA is to use the example calls provided in dem_example_sites.pro and dem_example_gridsites.pro. Users will need to change some values in these codes such as dates and directory structure. Note that you must be in your AIA/SDO SSW environment for these examples to work.

SITES
The core inversion function is dem_sites.pro. This will invert a single multi-channel measurement (i.e. one pixel). See dem_example_sites.pro for an example of using SITES for AIA data, particularly for looping through multiple pixels to create a DEM map of larger regions.

GRID-SITES
The Grid-SITES function is dem_gridsites.pro. This will invert multi-channel measurements for multiple pixels. See dem_example_gridsites.pro for an example of using Grid-SITES for AIA data. Grid-SITES uses the dem_sites.pro function to invert the DEMs, but is far more efficient, particularly for larger images or time series. 

ERRORS, PREPARATION OF DATA, AND TEMPERATURE RESPONSES
SITES requires the estimated errors of both the measurement and the response functions. The procedure dem_openaiafiles_sites.pro provides the measurement data cube in multiple channels, plus the errors. It can also combine multiple consecutive observations to improve the signal to noise (and other convenient tasks, for example extracting a  region). See the dem_example_sites procedure for an example on using dem_openaiafiles_sites. 

The dem_aiainterpolresponse_sites.pro procedure takes the structure returned by aia_get_response (standard AIA SSW calibration/analysis software) and interpolates the response curves to the user-defined temperature bins. It also estimates the uncertainty in the responses. 

OTHER INSTRUMENTS
In principle, if one has multi-wavelength measurements, appropriate temperature response curves, and measurement/response uncertainties, SITES (and Grid-SITES) should work. Note that SITES has been described exclusively in the context of AIA in [1], and I would be very interested to see results applied to other instruments. See contact below.

INVERSION OPTIONS
There are only 3 user-set parameters that effect the core inversion routine. These are:
1) Number of temperature bins. I usually set this at around 45.
2) Range of temperatures. I usually set a range of around 0.5 to 8MK, on a regular logarithmic temperature scale. I have not tested the code on flare-like temperatures. See [1] for a discussion of low <1MK temperatures in the context of AIA (basically never go below 0.6MK with AIA only).
3) Width of smoothing kernel. The default setting has been found through trial and error - a compromise between stable solutions and over-smoothing the results. Users can provide their own smoothing kernel if they wish to experiment with this parameter.

CONTACT
Huw Morgan, Aberystwyth University, hmorgan@aber.ac.uk

