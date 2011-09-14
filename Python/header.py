#!/usr/bin/env python
# encoding: utf-8
"""
header.py

Created by Jonathan Whitmore on 2011-07-12.
Copyright (c) 2011. All rights reserved.
"""

import sys
import os
import csv

def main():
  """
  header.py grabs the header information from a file and stores it into a dictionary.
  """
  # Need to run: 
  # grep -v 'COMMENT' headerFile | awk '{if ($2 >0) print} ' > cleanedFile
  inFile = 'test.2092'
  header = {} # create/clear header dictionary
  for row in csv.reader(open(inFile),delimiter='='): # note
    header[row[0].strip()]= row[1].strip() # removes unnecessary white space
    # should check if float, bool, string
    # strip off "double quotes"
  pass

if __name__ == '__main__':
  main()

