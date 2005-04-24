#!/usr/bin/env python

"""dvar.py   calculate variance of test versus target distances
Usage: dvar.py [options] target.dst test.dst

both target.dst and test.dst can be simple lists of distances or list of
coordinates in plain or atomeye format.  Number of distances in target.dst
can be larger than in test.dst.

Options:
  -n, --noscale   do not scale test distances to target distance list
"""

__id__ = "$Id$"


