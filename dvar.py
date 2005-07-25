#!/usr/bin/env python

"""dvar.py   calculate variance of test versus target distances
Usage: dvar.py [options] target.dst test1.dst [test2.dst]...

both target.dst and test.dst can be simple lists of distances or list of
coordinates in plain or atomeye format.  Number of distances in target.dst
can be larger than in test.dst.

Options:
  -m, --multiple  allow arbitrary multiplicity for test distances
  -r, --rescale   apply least squares scaling of test to target distances
  -h, --help      display this message
  -V, --version   show script version
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

def parseXYZ(lines):
    """parse list of line strings
    return list of coordinates"""
    rlist = []
    for l in lines:
        rlist.append( [ float(w) for w in l.split()[0:3] ] )
    return rlist

def parseAtomEye(lines):
    """parse a list of lines obtained from atomeye file
    return list of coordinates"""
    metrics = [  [ None, None, None ]  for i in range(3)  ]
    import re
    lstart = 0
    for l in lines[lstart:]:
        lstart += 1
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
    for l in lines[lstart:]:
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

def getDistancesFrom(filename):
    """open and detect format of specified distance or coordinate file
    return sorted list of distances"""
    file = open(filename)
    lines = file.readlines()
    import re
    if [ l for l in lines[:1000] if re.match('Number of particles', l) ]:
        dlist = coordinatesToDist(parseAtomEye(lines))
    else:
        last = len(lines)
        for first in range(last):
            try:
                xl = [ float(w) for w in lines[first].split() ]
            except ValueError:
                xl = []
            if len(xl):
                break
        # check if we have dumb XYZ file
        if len(xl) == 1 and xl[0] == round(xl[0]) and first+1 < last:
            try:
                xl1 = [ float(w) for w in lines[first+1].split() ]
            except:
                raise RuntimeError, "invalid data format at line %i" % \
                        (first+2, )
            # here the first line is number of records
            if len(xl1) > 1:
                first += 1
                last = first + int(xl[0])
                xl = xl1
        if len(xl) == 1:
            dlist = [ float(l) for l in lines[first:last] ]
            dlist.sort()
        elif len(xl) in (2, 3):
            dlist = coordinatesToDist(parseXYZ(lines[first:last]))
        else:
            raise RuntimeError, "invalid data format in " + filename
    return dlist

def findNearest(x, c):
    """find nearest item in the sorted list x to constant c
    returns index to the nearest value"""
    lo = 0; hi = len(x)-1
    while x[lo] < c < x[hi] and lo < hi-1:
        mi = (lo+hi)/2
        if x[mi] < c:
            lo = mi
        else:
            hi = mi
    if c < x[lo] and lo > 0 and c-x[lo-1] < x[lo]-c:
        idx = lo - 1
    elif c > x[hi] and hi+1 < len(x) and c-x[hi] < x[hi+1]-c:
        idx = hi+1
    elif c-x[lo] <= x[hi]-c:
        idx = lo
    else:
        idx = hi
    return idx

def dvar(dTarget, dTest, rescale = False, multiple = False):
    """calculate variance of two sorted distance lists
    return float"""
    dTest.sort()
    if multiple:
        # get a unique list of distances
        dtgt0 = dict(zip(dTarget, [1]*len(dTarget)))
        dtgt0 = dtgt0.keys()
        dtgt0.sort()
        dtgt1 = [ dtgt0[findNearest(dtgt0,t)] for t in dTest ]
    elif len(dTarget) < len(dTest):
        raise RuntimeError, "test distance list too long"
    elif len(dTarget) == len(dTest):
        dtgt1 = dTarget
    else:
        dtgt1 = []
        dtgt0 = list(dTarget)
        for t in dTest:
            i = findNearest(dtgt0, t)
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
    if len(dTest):
        v /= len(dTest)
    return v

def usage(style = None):
    import os.path
    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("dvar.py", myname)
    if style == 'brief':
        msg = msg.split("\n")[1] + "\n" + \
                "Try `%s --help' for more information." % myname
    print msg
    return

import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "mrhV", \
                ["multiple", "rescale", "help", "version"])
    except getopt.GetoptError, errmsg:
        print >> sys.stderr, errmsg
        sys.exit(2)
    # process options
    multiple, rescale = (False, False)
    for o, a in opts:
        if o in ("-m", "--multiple"):
            multiple = True
        if o in ("-r", "--rescale"):
            rescale = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-V", "--version"):
            print __id__
            sys.exit()
    if len(args) < 2:
        usage('brief')
        sys.exit()
    try:
        dTarget = getDistancesFrom(args[0])
    except IOError, (errno, errmsg):
        print >> sys.stderr, "%s: %s" % (args[0], errmsg)
        sys.exit(1)
    maxflen = max([ len(a) for a in args[1:] ])
    for f in args[1:]:
        try:
            v = dvar(dTarget, getDistancesFrom(f), rescale, multiple)
            if len(args) > 2:
                print f.ljust(maxflen+1), v
            else:
                print "%.8g" % v
        except RuntimeError, errmsg:
            print >> sys.stderr, "%s: %s" % (f, errmsg)
            sys.exit(1)
        except IOError, (errno, errmsg):
            print >> sys.stderr, "%s: %s" % (f, errmsg)
            sys.exit(1)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
