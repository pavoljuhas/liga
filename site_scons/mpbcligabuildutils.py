#!/usr/bin/env python

'''Utility functions used by scons and sphinx for version extraction.
'''

import os
import re

MYDIR = os.path.dirname(os.path.abspath(__file__))

def gitinfo():
    'Extract dictionary of version data from git records.'
    from subprocess import Popen, PIPE
    global _cached_gitinfo
    if _cached_gitinfo is not None:
        return _cached_gitinfo
    nullfile = open(os.devnull, 'w')
    kw = dict(stdout=PIPE, stderr=nullfile, cwd=MYDIR,
              universal_newlines=True)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    if proc.wait():
        _cached_gitinfo = {}
        return gitinfo()
    proc = Popen(['git', 'log', '-1', '--format=%H%n%ci%n%at%n%an'], **kw)
    glog = proc.stdout.read()
    words = desc.strip().split('-')
    vtag = words[0].lstrip('v')
    vnum = int((words[1:2] or ['0'])[0])
    rv = {}
    rv['commit'], rv['date'], rv['timestamp'], rv['author'] = [
            s.strip() for s in glog.strip().split('\n')]
    rv['timestamp'] = int(rv['timestamp'])
    mx = re.match(r'^(\d+)\.(\d+)(?:\.(\d+))?((?:[ab]|rc)\d*)?', vtag)
    rv['major'] = int(mx.group(1))
    rv['minor'] = int(mx.group(2))
    rv['micro'] = int(mx.group(3) or 0)
    rv['prerelease'] = mx.group(4)
    rv['patchnumber'] = vnum
    rv['version'] = vtag
    if vnum:
        rv['version'] += '.post' + str(vnum)
    _cached_gitinfo = rv
    return gitinfo()

_cached_gitinfo = None
