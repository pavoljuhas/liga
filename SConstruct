import os

# copy system environment variables related to compilation
DefaultEnvironment(ENV={
        'PATH' : os.environ['PATH'],
        'PYTHONPATH' : os.environ.get('PYTHONPATH', ''),
        'CPATH' : os.environ.get('CPATH', ''),
        'LIBRARY_PATH' : os.environ.get('LIBRARY_PATH', ''),
        'LD_LIBRARY_PATH' : os.environ.get('LD_LIBRARY_PATH', ''),
    }
)

# Create construction environment
env = DefaultEnvironment().Clone()

# Variables syntax was changed in 0.98
env.EnsureSConsVersion(0, 98)

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
vars.Update(env)
vars.Add(PathVariable('bindir',
    'installation directory for executable [prefix/bin]',
    env['prefix'] + '/bin'))
vars.Update(env)
env.Help(vars.GenerateHelpText(env))

builddir = env.Dir('build/%s' % env['build'])

Export('env')

env.SConscript('src/SConscript.main', variant_dir=builddir, duplicate=0)

# vim: ft=python
