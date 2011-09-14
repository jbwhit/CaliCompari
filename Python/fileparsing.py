#!/usr/bin/env python
# encoding: utf-8
"""
template-fitting.py

Created by Jonathan Whitmore on 2011-02-15
"""
version = 0.0

from config_calibration import *

OutputCalibrationDir = "Ascii_Output/Calibration/"
OutputResolutionDir = "Ascii_Output/Resolution/"

allFiles = glob.glob(OutputCalibrationDir + "*")

print "There are", len(allFiles), "files to analyze."

namingSections = len(allFiles[0].split("/")[-1].split(".")[:-1])

print "Each file has", namingSections, "naming sections."

print allFiles[0].split("/")[-1].split(".")[:-1]

# for rootString in allFiles:
#   print glob.glob(rootString.split("/")[-1].split(".") + "*")

# len(CalibrationFile.split("/")[-1].split("."))
# CalibrationFile.split("/")[-1].split(".")[:-1]
# 
# marc = CalibrationFile.split("/")[-1].split(".")[:-1]
# glob.glob("Ascii_Output/Calibration/" + marc[0] + "." + marc[1] + "*")


def unique(items):
    found = set([])
    keep = []
    for item in items:
        if item not in found:
            found.add(item)
            keep.append(item)

    return keep

print unique([1, 1, 2, 'a', 'a', 3])

parseList = [rootString.split("/")[-1].split(".")[:-1] for rootString in allFiles]


class MultiDict(dict):
  def __getitem__(self, key):
    if key in self:
      return dict.__getitem__(self, key)
    result = []
    for complete_key in sorted(self.keys()):
      if complete_key.startswith(key):
          result.extend(self[complete_key])
    return result

try1 = MultiDict()
avwav = MultiDict()
avcal = MultiDict()
averr = MultiDict()
avpix = MultiDict()
avres = MultiDict()
avres_err = MultiDict()

for files in allFiles:  
  reader = csv.reader(open(files),delimiter=' ')
  try1[files] = []
  for row in reader:
    try1[files].append(float(row[0]))

for files in allFiles:
  reader = csv.reader(open(files), delimiter=' ')
  avwav[files] = []
  avcal[files] = []
  averr[files] = []
  avres[files] = []
  avres_err[files] = []
  avpix[files] = []
  for row in reader:
    avwav[files].append(float(row[0]))
    avcal[files].append(float(row[1]))
    averr[files].append(float(row[2]))
    avres[files].append(float(row[3]))
    avres_err[files].append(float(row[4]))
    if len(row) > 5:
      avpix[files].append(float(row[5]))


objectList = list(set([x[0] for x in parseList])) # probably should scrape the beginning for just object
objectList.sort()
arcList = list(set([x[1] for x in parseList]))
arcList.sort()
exposureList = list(set([x[2] for x in parseList]))
exposureList.sort()
chipList = list(set([x[3] for x in parseList]))
chipList.sort()

def unholymess():
  """docstring for unholymess"""
  arcthing = {}
  jeff = avwav.keys()
  jeff.sort()
  for thing in arcList:
      print thing
      arcthing[thing] = []  
  for thing in jeff:
      arcthing[thing.split('.')[1]].append(thing)
  for counter, listthing in enumerate(arcList):
    for files in arcthing[listthing]:
        pl.plot(avwav[files],avcal[files])
    pl.title(listthing.title())
    pl.savefig('PDF_Output/' + listthing + ".pdf")
    pl.close()
  pass