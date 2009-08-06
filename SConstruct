Environment().EnsureSConsVersion(0, 98)

import os

# Customizable compile variables
vars = Variables('sconsvars.py')

vars.Add('tests', 'Custom list of unit test sources', None)
vars.Add(EnumVariable('build',
    'compiler settings', 'debug',
    allowed_values=('debug', 'fast')))
vars.Add(BoolVariable('profile',
    'build with profiling information', False))
vars.Add(PathVariable('prefix',
    'installation prefix directory', '/usr/local'))
vars.Add(PathVariable('libdir',
    'object code library directory [prefix/lib]', None))
vars.Add(PathVariable('bindir',
    'installation directory for executable [prefix/bin]', None))

# Default build environment
env = DefaultEnvironment(variables=vars, ENV={
    'PATH' : os.environ['PATH'],
    'PYTHONPATH' : os.environ['PYTHONPATH'],
    })
Help(vars.GenerateHelpText(env))

env.ParseConfig("python-config --includes")
env.ParseConfig("python-config --ldflags")
env.ParseConfig("gsl-config --cflags --libs")

env.AppendUnique(LIBS='libboost_python')

if env['build'] == 'debug':
    env.Append(CCFLAGS='-g')
    exesuffix = ''
elif env['build'] == 'fast':
    env.AppendUnique(CCFLAGS=['-O3', '-ffast-math'])
    env.AppendUnique(CPPDEFINES='NDEBUG')
    exesuffix = '-fast'

if env['profile']:
    env.AppendUnique(CCFLAGS='-pg')
    env.AppendUnique(LINKFLAGS='-pg')

env.Append(CCFLAGS='-Wall')
env.PrependUnique(CPPPATH=['#/'])

# make env available to subsidiary SConscripts
Export('env')

# Define lists for storing library source and include files.
def isLibSource(f):
    f = str(f)
    rv = f[:1].isupper() and not f.startswith('Test') and f != 'Version.cpp'
    return rv

env['binaries'] = []
env['lib_sources'] = filter(isLibSource, Glob('*.cpp'))
# This SConscript updates Version.cpp and adds it to lib_sources
SConscript('SConscript.version')
env['lib_objects'] = map(env.Object, env['lib_sources'])

# Top Level Targets ----------------------------------------------------------

# mpbcliga -- application
mpbcliga = env.Program('mpbcliga' + exesuffix,
        ['mpbcliga.cpp'] + env['lib_objects'])
Alias('mpbcliga', mpbcliga)
env['binaries'] += mpbcliga
env['mpbcliga'] = mpbcliga

# This SConscript defines all test targets
SConscript('SConscript.tests')

# Installation targets.

prefix = env['prefix']

# install-bin
bindir = env.get('bindir', os.path.join(env['prefix'], 'bin'))
get_target_path = lambda f : os.path.join(bindir, f.path)
bin_targets = map(get_target_path, env['binaries'])

Alias('install-bin', InstallAs(bin_targets, env['binaries']))

# install
Alias('install', ['install-bin'])

# vim: ft=python
