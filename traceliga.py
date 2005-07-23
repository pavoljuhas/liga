#!/usr/bin/env python

"""traceliga.py   trace development of specified team in liga system
Usage: traceliga.py liga.log [season level place]

liga.log must be created using verbose run of gizaliga, 
"season level place" define the last traced team, by default the last
team in liga.log

Options:
  -f, --frameselection   print trace in gizaliga parfile format
  -h, --help      display this message
  -V, --version   show script version
"""

__id__ = "$Id$"

class TeamIdType:
    'class TeamIdType stores season, level and place of liga team'
    def __init__(self, season, level, place):
        self.season = season
        self.level = level
        self.place = place
        return
    def __eq__(self, other):
        return self.season == other.season and \
                self.level == other.level  and \
                self.place == other.place
    # definition of __eq__ requires to define __hash__ as well
    def __hash__(self):
        return hash( (self.season, self.level, self.place) )
    def __str__(self):
        return "%i %i %i" % (self.season, self.level, self.place)

def readLog(filename):
    """read verbose liga info from specified log file
    return reverse ordered list of verbose message lines"""
    file = open(filename)
    import re
    verbosepattern = re.compile("""^\d+(?: PUSH .* TO | SWAP .* WITH )""")
    vmsgs = []
    for l in file:
        if not verbosepattern.match(l): continue
        vmsgs.append(l.rstrip())
    vmsgs.reverse()
    file.close()
    return vmsgs

def findTrace(tracedTeam, vmsgs):
    '''find trace of the specified [default=last] team from verbose messages
    return ordered list of TeamIdType instances'''
    from copy import copy
    if not tracedTeam:
        f = vmsgs[0].split()
        tracedTeam = TeamIdType( int(f[0]), int(f[5]), int(f[6]) )
    trace = [tracedTeam]
    for l in vmsgs:
        f = l.split()
        season = int(f[0])
        if season > trace[-1].season: continue
        action = f[1]
        (l1, p1, l2, p2) = [int(f[i]) for i in (2, 3, 5, 6)]
        lastlp = (trace[-1].level, trace[-1].place)
        if action == 'SWAP':
            if lastlp == (l1, p1):
                trace.append( TeamIdType(season, l2, p2) )
            if lastlp == (l2, p2):
                trace.append( TeamIdType(season, l1, p1) )
        elif action == 'PUSH' and lastlp == (l2, p2):
            trace.append( TeamIdType(season, l1, p1) )
    trace.reverse()
    # erase all up to and including the last occurence of the 0th level
    idx0 = [ i for i in range(len(trace)) if trace[i].level == 0 ]
    if idx0: del trace[:idx0[-1]+1]
    # remove duplicate entries
    selected = {}
    uniqued_trace = []
    for t in trace:
        if not t in selected:
            uniqued_trace.append(t)
            selected[t] = True
    return uniqued_trace

import sys
import getopt

def usage(style = None):
    """print usage information, usage('brief') gives shorter info"""
    import os.path
    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("traceliga.py", myname)
    if style == 'brief':
        msg = msg.split("\n")[1] + "\n" + \
                "Try `%s --help' for more information." % myname
    print msg
    return

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "fhV", \
                ["frameselection", "help", "version"])
    except getopt.GetoptError, errmsg:
        print >> sys.stderr, errmsg
        sys.exit(2)
    # process options
    tracefmt = "plain"
    for o, a in opts:
        if o in ("-f", "--frameselection"):
            tracefmt = "frameselection"
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-V", "--version"):
            print __id__
            sys.exit()
    if len(args) < 1:
        usage('brief')
        sys.exit()
    elif len(args) not in (1, 4):
        print >> sys.stderr, "incorrect number of arguments"
        sys.exit(2)
    logfile = args[0]
    tracedTeam = None
    if len(args) == 4:
        try:
            slp = [ int(a) for a in args[1:5] ]
            tracedTeam = TeamIdType(*slp)
        except ValueError, errmsg:
            print >> sys.stderr, errmsg
            sys.exit(2)
    try:
        vmsgs = readLog(logfile)
    except IOError, (errno, errmsg):
        print >> sys.stderr, "%s: %s" % (logfile, errmsg)
        sys.exit(2)
    if len(vmsgs) == 0:
        sys.exit()
    trace = findTrace(tracedTeam, vmsgs)
    # print it out
    if tracefmt == 'plain':
        for t in trace:
            print t
    elif tracefmt == 'frameselection':
        lines = [ 'frameselection=' ] + \
                [ str(t) for t in trace ]
        print ' \\\n    '.join(lines)
    return

if __name__ == "__main__":
    main(sys.argv[1:])
