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

def posToDist(r):
    """calculate sorted list of distances from a given coordinate matrix
    return builtin list"""

    from math import sqrt
    d = []
    for i in range(len(r)):
        for j in range(i+1,len(r)):
            dij2 = sum([ pow(xi-xj,2) for (xi,xj) in zip(r[i], r[j]) ])
            d.append(sqrt(dij2))
    d.sort()
    return d

def readPlainPos(file):
    """read plain coordinates from open file
    return list of coordinates"""
    
    rlist = []
    file.seek(0)
    for l in file:
        if l[0] == '#':  continue
        rlist.append( [ float(w) for w in l.split()[0:3] ] )
    rarray = rlist
    return rarray

if __name__ == "__main__":
    file = open("/u24/juhas/programs/BGA/solids/icosahedron.xyz")
    r = readPlainPos(file)
    file.close()
    d = posToDist(r)
    for v in d:
        print v
