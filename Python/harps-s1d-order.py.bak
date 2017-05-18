#!/usr/bin/env python
# encoding: utf-8
"""
harps.py

Turn order merged HARPS output into a split by order output.

Created by Jonathan Whitmore on 2011-10-04.
Copyright (c) 2011. All rights reserved.
"""

import sys
import getopt
import numpy as np
import pyfits as pf

help_message = '''
To change HARPS exposures into a format that 
CaliCompari can understand.
'''


class Usage(Exception):
  def __init__(self, msg):
    self.msg = msg


def main(argv=None):
  if argv is None:
    argv = sys.argv
  try:
    try:
      opts, args = getopt.getopt(argv[1:], "hi:o:v", ["help", "inputfile=", "output="])
    except getopt.error, msg:
      raise Usage(msg)
    
    # option processing
    for option, value in opts:
      if option == "-v":
        verbose = True
      if option in ("-h", "--help"):
        raise Usage(help_message)
      if option in ("-o", "--output"):
        outputFile = value
      if option in ("-i", "--inputfile"):
        inputfile = value
  
    hdulist = pf.open(inputfile) # Read in input file
    # startingPixel = int(hdulist[0].header['CRPIX1'])
    # print "Starting pixel: ", startingPixel
    startingWavelength = hdulist[0].header['CRVAL1']
    print "Starting wavelength: ", startingWavelength
    stepSize = hdulist[0].header['CDELT1']
    print "Step size: ", stepSize
  
    fluxArray = hdulist[0].data
    wavelengthArray = np.arange(startingWavelength, len(fluxArray)*stepSize + startingWavelength, stepSize)
    
    if len(wavelengthArray) == len(fluxArray):
      print "Wavelength and flux arrays are the same length"
    else:
      print "Problems: wavelength and flux arrays are unequal lengths!"
      if len(wavelengthArray) > len(fluxArray):
        repeat = len(wavelengthArray) - len(fluxArray)
        wavelist = list(wavelengthArray)
        for x in range(repeat):
          wavelist.pop(-1)
        wavelengthArray = np.array(wavelist)
      print len(wavelengthArray)
      print len(fluxArray)
  
    Overlap = 50.0 # Angstroms
    orders = 20 # number of orders to break it into. 
    slicelength = len(wavelengthArray)/orders
    wav = []
    flx = []
    err = []
    outputFile = str("jw_" + inputfile)
    for i in range(orders): 
      temp1 = np.array([wavelengthArray[slicelength * i: slicelength*(i+1) + 5000], \
        fluxArray[slicelength * i: slicelength*(i+1) + 5000],\
        np.sqrt(np.abs(fluxArray[slicelength * i: slicelength*(i+1) + 5000]))])
      pf.append(outputFile, temp1)
      
  except Usage, err:
    print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
    print >> sys.stderr, "\t for help use --help"
    return 2


if __name__ == "__main__":
  sys.exit(main())