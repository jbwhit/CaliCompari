#!/usr/bin/env python
# encoding: utf-8
"""
vpfitinterface.py

Created by Jonathan Whitmore on 2011-07-13.
Copyright (c) 2011 All rights reserved.

Provides for a python interface for setting up and interfacing with vpfit.

"""

import sys
import shutil
import os
import glob
import csv
import numpy as np
import scipy as sp
from subprocess import *
import time
import datetime
from optparse import OptionParser
import argparse
from ConfigParser import RawConfigParser
import random as ra
import tempfile
import itertools
import fileinput

FWHM = 6.5

rootName='test'

def main(rootName=rootName):
  # root name, ex: a1
  # Fit setup file, ex: a1.setup
  # vpfit commands file, ex: a1.vpcommands
  # Move fort.26 to have same name as input, ex: a1.results
  pass

def options():
  """docstring for options
  I - interactive setup and fit
  F - run from an input file
  D - display profiles from input file
  """
  if len(glob.glob(rootName + '.setup')) > 0:
    pass
  else:
    print "No basis setup file."
  pass

def commands():
  """docstring for commands"""
  commandFileName = rootName + '.vpcommands'
  commandFileHandle = open(commandFileName, 'w')
  print >>commandFileHandle, 'f ad\n\n\n', commandFileName, '\n'
  commandFileHandle.close()
  pass

# check if zero/continuum adjust, and if so, continue run

if __name__ == '__main__':
  main()