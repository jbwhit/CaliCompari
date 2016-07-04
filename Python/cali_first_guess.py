#!/usr/bin/env python

# Usage: cali_first_guess [shift] [slope] [sigma] [multiple] [offset]
# Units are:              [Ang.]  [???]   [km/s?] [None]     [Norm.]
# Reasonable values:      0.001   -0.002  3.0     1.37       0.002

import sys
import json

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

first_guesses = {}
first_guesses.update({'shift':float(sys.argv[1]),
                      'fix_shift':False,
                      'limit_shift':(-1.5, 1.5),
                      'error_shift':0.03})
first_guesses.update({'slope':float(sys.argv[2]),
                      'fix_slope':False,
                      'limit_slope':(-2.0, 2.0),
                      'error_slope':0.04})
first_guesses.update({'sigma':float(sys.argv[3]),
                      'fix_sigma':False,
                      'limit_sigma':(1.0, 10.0),
                      'error_sigma':0.2})
first_guesses.update({'multiple':float(sys.argv[4]),
                      'fix_multiple':False,
                      'limit_multiple':(0.1, 20.0),
                      'error_multiple':0.03})
first_guesses.update({'offset':float(sys.argv[5]),
                      'fix_offset':False,
                      'limit_offset':(-2.0, 2.0),
                      'error_offset':0.03})
first_guesses.update({'minuit':0, 'fix_minuit':True})

with open("first_guesses.json", 'w') as file_handle:
    json.dump(first_guesses, file_handle, indent=2, sort_keys=True)
