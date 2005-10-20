#!/usr/bin/env python

"""uniqdst.py   list unique distances and their multiplicities
Usage: uniqdst.py [options] distanceList.dst

Options:
  -r, --resolution  [0.01] resolution for equal distances
"""

__id__ = "$Id$"

def readDistancesFrom(filename):
    """readDistancesFrom(filename)  open and read specified distance file

    return sorted list of distances"""
    dlist = [ float(line) for line in open(filename) ]
    return dlist

def findDuplicates(dlist, resolution):
    """findDuplicates(dlist, resolution)

    find duplicate distances in sorted distance list
    returns list of list of duplicates"""
    # fill with dummy first value
    dupslist = [ [-2.0*resolution] ]
    for d in dlist:
        if (d - dupslist[-1][0]) > resolution:
            dupslist.append([d])
        else:
            dupslist[-1].append(d)
    # remove dummy value
    del dupslist[0]
    return dupslist

def usage(style = None):
    import os.path
    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("uniqdst.py", myname)
    if style == 'brief':
        msg = msg.split("\n")[1] + "\n" + \
                "Try `%s --help' for more information." % myname
    print msg
    return

import sys
import getopt

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "r:hV", \
                ["resolution", "help", "version"])
    except getopt.GetoptError, errmsg:
        print >> sys.stderr, errmsg
        sys.exit(2)
    # process options
    resolution = 0.01
    for o, a in opts:
        if o in ("-r", "--resolution"):
            resolution = float(a)
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-V", "--version"):
            print __id__
            sys.exit()
    if len(args) > 1:
        print "too many distance files"
        usage('brief')
        sys.exit(1)
    try:
        dlist = readDistancesFrom(args[0])
    except IOError, (errno, errmsg):
        print >> sys.stderr, "%s: %s" % (args[0], errmsg)
        sys.exit(1)
    dupslist = findDuplicates(dlist, resolution)
    for dups in dupslist:
        print sum(dups)/len(dups), len(dups)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
