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
vars.Add(PathVariable('bindir',
    'installation directory for executable [prefix/bin]', None))

# Default build environment
env = DefaultEnvironment(variables=vars, ENV={
    'PATH' : os.environ['PATH'],
    'PYTHONPATH' : os.environ['PYTHONPATH'],
    'CPATH' : os.environ['CPATH'],
    'LIBRARY_PATH' : os.environ['LIBRARY_PATH'],
    'LD_LIBRARY_PATH' : os.environ['LD_LIBRARY_PATH'],
    })
Help(vars.GenerateHelpText(env))

builddir = Dir('build/%s' % env['build'])

Export('env')

env.SConscript('src/SConscript.main', variant_dir=builddir, duplicate=0)
env.Clean(env['binaries'], builddir)
