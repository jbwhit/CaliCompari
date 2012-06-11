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

def twoPassSmooth(iow, iof, binsize1, binsize2, nlarge):
  """Pass the data through twice..."""
  fiow = []
  fiof = []
  for i in xrange(0, len(iof), binsize1):
    fiow.append(np.average(iow[i:i + binsize1]))
    fiof.append(np.average(heapq.nlargest(nlarge, iof[i:i + binsize1])))
  
  tck = splrep(fiow, fiof)
  iof = iof/splev(iow,tck)
  iow = iow * 10.0
  # pl.plot(fiow, fiof)
  # pl.plot(iow, iof, color="black")
  siow = []
  siof = []
  for i in xrange(0, len(iof), binsize2):
    siow.append(np.average(iow[i:i + binsize2]))
    siof.append(np.average(heapq.nlargest(nlarge, iof[i:i + binsize2])))
  
  print siof
  tck = splrep(siow, siof)
  iof = iof/splev(iow, tck)
  # print tck 
  # pl.plot(iow, splev(iow, tck), color="blue")
  # pl.plot(iow, iof, color="red")
  # pl.plot(siow, siof)
  # pl.show()
  return iow, iof

def main():
  # smoothContinuum(200, 1)
  iow, iof, inten, kinten = np.loadtxt('sao2010.solref.txt',unpack=True, comments='C')
  iow, iof = twoPassSmooth(iow, iof, 300, 300, 5)
  np.savetxt('sao2010.contin.ascii', zip(iow, iof))
  # pl.plot(tiow, tiof)
  # pl.savefig('plot.pdf')
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


