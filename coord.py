#!/usr/bin/env python

"""coord.py   show atom coordination numbers for given structure
Usage: coord.py [options] structfile1 [structfile2]...

structfile can be in plain or atomeye format.

Options:
  -c, --cutoff=NUM  change bond cutoff from atomeye default of 1.73823
  -l, --list      list neighbors of every atom
  -h, --help      display this message
  -v, --version   show script version
"""

__id__ = "$Id$"

def findNeighbors(r, cutoff):
    """for every position find neighbors closer than cutoff
    return list of neighbor lists"""
    from math import sqrt
    nbl = [ list([]) for i in range(len(r)) ]
    for i in range(len(r)):
        for j in range(i+1,len(r)):
            dij2 = sum([ pow(xi-xj,2) for (xi,xj) in zip(r[i], r[j]) ])
            dij = sqrt(dij2)
            if dij < cutoff:
                nbl[i].append(j)
                nbl[j].append(i)
    return nbl

def getCoordinatesFrom(filename):
    """open and detect format of specified coordinate file
    return list of coordinates"""
    import dvar
    file = open(filename)
    lhead = file.readlines(2000)
    lhead = [ l.rstrip() for l in lhead[0:-1] ]
    file.seek(0)
    import re
    if re.search('^Number of particles', '\n'.join(lhead), re.M):
        r = dvar.readAtomEye(file)
    else:
        headlines = 0
        for l in lhead:
            numbercount = 0
            for w in l.split():
                try:
                    x = float(w)
                    numbercount += 1
                except ValueError:
                    numbercount = 0
                    break
            if numbercount:
                break
            headlines += 1
        if numbercount == 3:
            r = dvar.readXYZ(file, headlines)
        else:
            raise RuntimeError, "invalid data format in " + filename
    file.close()
    return r

def usage(style = None):
    import os.path
    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("coord.py", myname)
    if style == 'brief':
        msg = msg.split("\n")[1] + "\n" + \
                "Try `%s --help' for more information." % myname
    print msg
    return

import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "c:lhv", \
                ["cutoff", "list", "help", "version"])
    except getopt.GetoptError, errmsg:
        print >> sys.stderr, errmsg
        sys.exit(2)
    # process options
    cutoff = 1.73823
    listnb = False
    for o, a in opts:
        if o in ("-c", "--cutoff"):
            cutoff = float(a)
        elif o in ("-l", "--list"):
            listnb = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-v", "--version"):
            print __id__
            sys.exit()
    if len(args) < 1:
        usage('brief')
        sys.exit()
    maxflen = max([ len(a) for a in args ])
    for f in args:
        try:
            nbl = findNeighbors(getCoordinatesFrom(f), cutoff)
            if listnb:
                if len(args) > 2:
                    print f
                for i in range(len(nbl)):
                    print "%-6i%s" % ( i, " ".join([str(n) for n in nbl[i]]) )
            else:
                if len(args) > 2:
                    print f.ljust(maxflen+1),
                zmult = {}
                for z in [len(l) for l in nbl]:
                    zmult[z] = zmult.get(z, 0) + 1
                zks = zmult.keys()
                zks.sort()
                print " ".join([ "%ix%i" % (zmult[z], z) for z in zks ])
        except RuntimeError, errmsg:
            print >> sys.stderr, "%s: %s" % (f, errmsg)
            sys.exit(1)
        except IOError, (errno, errmsg):
            print >> sys.stderr, "%s: %s" % (f, errmsg)
            sys.exit(1)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
