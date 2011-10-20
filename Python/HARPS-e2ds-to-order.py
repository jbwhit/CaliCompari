#!/usr/bin/env python
# encoding: utf-8
"""
HARPS-e2ds-to-order.py

Created by Jonathan Whitmore on 2011-10-14.
Copyright (c) 2011. All rights reserved.
"""

import sys
import os
import argparse
import pyfits as pf
import numpy as np

help_message = '''
Takes reduced HARPS***e2ds_A.fits data and reads the header to output a fits file that has the wavelength per pixel per order.
'''


class Usage(Exception):
  def __init__(self, msg):
    self.msg = msg


def main(argv=None):
  parser = argparse.ArgumentParser(description='Process input file.')
  parser.add_argument('f', type=str, help='input a filename')
  args = parser.parse_args()
  inputfile = vars(args)['f']
  outputFITSfile = str('order_' + inputfile)
  hdu = pf.open(inputfile)
  print len(hdu[0].data)
  polynomialOrder = hdu[0].header['HIERARCH ESO DRS CAL TH DEG LL']
  orders = hdu[0].header['HIERARCH ESO DRS CAL LOC NBO']
  coefficients = (polynomialOrder + 1) * orders # number of coefficients per order
  allCoefficients = []
  for x in xrange(coefficients):
    allCoefficients.append(hdu[0].header['HIERARCH ESO DRS CAL TH COEFF LL' + str(x)])
  A = {}
  for y in xrange(0,len(allCoefficients),polynomialOrder+1):
    A[y/(polynomialOrder+1)] = allCoefficients[y:y+polynomialOrder+1]
  for order in range(orders):
    print "order: ", order
    wavelength = []
    for pixel in xrange(len(hdu[0].data[order])):
      temp = 0.0
      for x in range(polynomialOrder):
        temp += allCoefficients[order + 1 + x] * pixel ** x
        # temp += A[order][x] * pixel ** x
      wavelength.append(temp)
    pf.append(outputFITSfile, np.array([np.array(wavelength), hdu[0].data[order], np.sqrt(np.abs(hdu[0].data[order]))]))
  

if __name__ == "__main__":
  sys.exit(main())
