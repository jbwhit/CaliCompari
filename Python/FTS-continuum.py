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

def main():
  iow, iof = np.loadtxt('KPNO2010.ascii',unpack=True)
  binsize = 5000 # FTS pixels
  miow = []
  miof = []
  for i in xrange(0, len(iof), binsize):
      miow.append(iow[i+np.argmax(iof[i:i + binsize])])
      miof.append(np.max(iof[i:i + binsize]))

  tck = splrep(miow, miof)
  np.savetxt('KPNO2010_contin.ascii', zip(iow, iof/splev(iow,tck)))
  pass


if __name__ == '__main__':
  main()


