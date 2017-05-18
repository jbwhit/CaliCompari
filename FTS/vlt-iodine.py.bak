#!/usr/bin/env python
# vlt-iodine
import pyfits as pf
import numpy as np

# Download and gunzip from:
# http://www.eso.org/sci/facilities/paranal/instruments/uves/tools/uves_i2_70_wn_cor.fits.gz

hdu = pf.open("uves_i2_70_wn_cor.fits.gz")

# CRPIX1  =                   1. / Reference pixel                                
# CRVAL1  =        14999.9976299 / Coordinate at reference pixel                  
# CDELT1  =         0.0066704159 / Coordinate increment per pixel                 
# CTYPE1  = 'LINEAR          '   / Units of coordinate                            
# BUNIT   = '                '   / Units of data values            

starting_pixel = hdu[0].header['CRPIX1']
starting_wavenumber = hdu[0].header['CRVAL1']
delta = hdu[0].header['CDELT1']
    
fluxes = hdu[0].data
x = np.array(range(len(fluxes)))
    
wavelengths = 1.0 / (np.ones_like(x) * starting_wavenumber + x * delta)
# reverse order
wavelength = wavelengths[::-1] * 1.0e8 # make into Angstroms
flux = fluxes[::-1]

np.savetxt("vlt.2013-12-07.txt.gz", (wavelength, flux)) # natively compresses
