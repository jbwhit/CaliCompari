#!/usr/bin/env python
# encoding: utf-8
"""
config_calibration.py

Created by Jonathan Whitmore on 2010-06-14.
"""

from load_calibration import *


# ==============================
# = Configuration file parsing =
# ==============================

filename='config'
config = RawConfigParser()
config.read(filename)
error_lower = float(config.get('Parameters','error_lower'))
sn_lower = float(config.get('Parameters','sn_lower'))
# TODO Make this s_upper optional
s_upper = float(config.get('Parameters','s_upper'))
edgebuffer = int(config.get('Parameters','edgebuffer'))
telescope = config.get('Information','telescope')
astro_object = config.get('Information','astro_object')
spectrograph = config.get('Information','spectrograph')

bin_size = float(config.get('Parameters','bin_size'))  # 300.0 # in km/s
step = float(config.get('Parameters','step_size')) #  50.0 # km/s
gauss_elements = float(config.get('Parameters','gauss_elements'))
sigma = float(config.get('Parameters','sigma'))
i2exposures = config.get('Information','Exposures').strip().splitlines()
chips = config.get('Information','Multi_Chips').strip().splitlines()
destroy = config.get('Parameters','Wavelength Remove').strip().splitlines()
begin_kill_array = []
end_kill_array = []
for i in range(len(destroy)):
    t1, t2 = destroy[i].split()
    begin_kill_array.append(float(t1))
    end_kill_array.append(float(t2))

begin_kill_array = np.array(begin_kill_array)
end_kill_array = np.array(end_kill_array)

if spectrograph=="UVES": 
    print "Instrument: ", spectrograph
    FTSFile = config.get('Information','path_to_UVES_FTS')
if spectrograph=="HIRES":
    print "Instrument: ", spectrograph
    FTSFile = config.get('Information','path_to_HIRES_FTS')

# ==================================
# = Print input values to terminal =
# ==================================

print "Chips: ", chips
print "Telescope: ", telescope, " Object: ", astro_object
print "Error Lower Bound: ", error_lower
print "S/N Lower Bound: ", sn_lower
print "Edgebuffer: ", edgebuffer
print "FTS File: ", FTSFile
print "Bin size: ", bin_size, "km/s", "Step size: ", step, "km/s"
print "# of Elements in gaussian array: ", gauss_elements, "Initial sigma: ", sigma

# =================
# = Parse Options =
# =================

core_test = 0

parser = OptionParser()
parser.add_option("-e","--exposure", 
                  action="store", type="string", dest="exponame", 
                  help="name exposures", metavar="EXPO")
parser.add_option("-c","--core", 
                action="store", type="string", dest="corename", 
                help="Dual Core Advantage", metavar="CORE")

(options, args) = parser.parse_args()

# Optional use of dual core stuff
try: core_test += float(options.corename) 
except: pass

# Run with flag -c 1    or --core 1
core1 = i2exposures[:len(i2exposures)/2]
core2 = i2exposures[len(i2exposures)/2:]

if core_test == 1:
    exposures_to_analyze = core1
    print "Exposures: ", exposures_to_analyze
if core_test == 2:
    exposures_to_analyze = core2
    print "Exposures: ", exposures_to_analyze
if core_test != 1 and core_test !=2: 
    exposures_to_analyze = core1 + core2
    print "Exposures: ", exposures_to_analyze
