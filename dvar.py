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

def coordinatesToDist(r):
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

def readXYZ(file):
    """read plain coordinates from open file
    return list of coordinates"""
    rlist = []
    file.seek(0)
    for l in file:
        if l[0] == '#':
            continue
        rlist.append( [ float(w) for w in l.split()[0:3] ] )
    return rlist

def readAtomEye(file):
    """read coordinates from open atomeye file
    return list of coordinates"""
    metrics = [  [ None, None, None ]  for i in range(3)  ]
    file.seek(0)
    import re
    for l in file:
        l = l.rstrip()
        if re.match('#', l):
            continue
        elif re.match('Number of particles', l):
            Natoms = float( l.split()[-1] )
        elif re.match('H0\(\d,\d\)', l):
            i = int(l[3]) - 1
            j = int(l[5]) - 1
            m = re.compile('H0\(\d,\d\) *= *(\S*?) *A?$').match(l)
            metrics[i][j] = float(m.group(1))
        elif l == "" :
            break
    # here we should be at start of coordinates
    rlist = []
    for l in file:
        l = l.rstrip()
        w = l.split()
        if len(w) == 1:
            try:
                atomMass = float(w[0])
            except ValueError:
                atomType = w[0]
        else:
            # first three numbers are coordinates
            rhkl = [ float(w[i]) for i in range(3) ]
            r = [0.0, 0.0, 0.0]
            for i in range(3):
                for j in range(3):
                    r[i] += metrics[i][j]*rhkl[j]
            rlist.append(r)
    # all done here
    return rlist

def readDistanceList(file):
    """read list of distances from a plain file
    return sorted list of distances"""
    dlist = []
    file.seek(0)
    for l in file:
        if l[0] == '#':
            continue
        dlist.append( float(l.split()[0]) )
    dlist.sort()
    return dlist

def getDistancesFrom(filename):
    """open and detect format of specified distance or coordinate file
    return sorted list of distances"""
    file = open(filename)
    lhead = file.readlines(2000)
    lhead = [ l.rstrip() for l in lhead[0:-1] ]
    file.seek(0)
    import re
    if re.search('^Number of particles', '\n'.join(lhead), re.M):
        dlist = coordinatesToDist(readAtomEye(file))
    else:
        for l in lhead:
            numbercount = 0
            for w in l.split():
                try:
                    x = float(w)
                    numbercount += 1
                except ValueError:
                    pass
            if numbercount:
                break
        if numbercount == 1:
            dlist = readDistanceList(file)
        else:
            dlist = coordinatesToDist(readXYZ(file))
    file.close()
    return dlist

def dvar(dTarget, dTest, rescale = True):
    """calculate variance of two sorted distance lists
    return float"""
    dTest.sort()
    if len(dTarget) < len(dTest):
        raise RuntimeError, "test distance list too long"
    elif len(dTarget) == len(dTest):
        dtgt1 = dTarget
    else:
        dtgt1 = []
        dtgt0 = list(dTarget)
        for t in dTest:
            dtg = [abs(t-d) for d in dtgt0]
            i = dtg.index(min(dtg))
            dtgt1.append( dtgt0.pop(i) )
    # here we need to compare dTest and dtgt1
    # default scale is 1.0
    scale = 1.0
    # if required update scale
    if rescale:
        sxy = sxx = 0.0
        for i in range(len(dTest)):
            sxy += dTest[i]*dtgt1[i]
            sxx += dTest[i]*dTest[i]
        if sxx > 0.0:
            scale = sxy/sxx
    # now we can calculate dvar
    v = 0.0
    for i in range(len(dTest)):
        v += pow(scale*dTest[i] - dtgt1[i], 2)
    v /= len(dTest)
    return v

def usage():
    import os.path
    print __doc__.replace("dvar.py", os.path.basename(sys.argv[0]))

import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "n", ["noscale"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    if len(args) < 2:
        usage()
        sys.exit()
    elif len(args) > 2:
        print >> sys.stderr, "too many file arguments:"
        print >> sys.stderr, "\n    ".join(args)
        sys.exit(2)
    # process options
    rescale = True
    for o, a in opts:
        if o in ("-n", "--noscale"):
            rescale = False
    dTarget = getDistancesFrom(args[0])
    dTest = getDistancesFrom(args[1])
    try:
        dv = dvar(dTarget, dTest, rescale)
        print dv
    except RuntimeError, s:
        print >> sys.stderr, s
        sys.exit(1)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
