#!/usr/bin/env python
# encoding: utf-8
"""
FTS-continuum.py

Created by Jonathan Whitmore on 2011-09-30.
Copyright (c) 2011. All rights reserved.
"""

import sys
import os
import pylab as pl
import numpy as np
from scipy.interpolate import splrep, splev
import heapq

def smoothContinuum(binsize, nlarge, plot=False):
  """docstring for smoothContinuum"""
  iow, iof, inten, kinten = np.loadtxt('sao2010.solref.txt',unpack=True, comments='C')
  siow = []
  siof = []
  for i in xrange(0, len(iof), binsize):
    siow.append(np.average(iow[i:i + binsize]))
    siof.append(np.average(heapq.nlargest(nlarge, iof[i:i + binsize])))
  
  tck = splrep(siow, siof)
  print "Multiplying wavelength by 10 to convert to Angstroms."
  np.savetxt('sao2010.contin.ascii', zip(iow * 10.0, iof/splev(iow,tck)))
  if plot == True:
    plot(iow, iof/splev(iow,tck))
  pass

def main():
  smoothContinuum(200, 1)
  # iow, iof = np.loadtxt('KPNO2010.ascii',unpack=True)
  # binsize = 5000 # FTS pixels
  # miow = []
  # miof = []
  # for i in xrange(0, len(iof), binsize):
  #     miow.append(iow[i+np.argmax(iof[i:i + binsize])])
  #     miof.append(np.max(iof[i:i + binsize]))
  # 
  # tck = splrep(miow, miof)
  # np.savetxt('KPNO2010_contin.ascii', zip(iow, iof/splev(iow,tck)))
  pass


if __name__ == '__main__':
  main()


