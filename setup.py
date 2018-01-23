#!/usr/bin/env python

package_name = "MMTK"

import os
import sys
import types
import ctypes
import ctypes.util
import distutils.sysconfig
import numpy.distutils.misc_util

from setuptools import setup, find_packages
from Cython.Build import cythonize
from Scientific import N

from distutils.core import Command, Extension
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.command.install_data import install_data
from distutils import dir_util
from distutils.filelist import FileList, translate_pattern

from glob import glob
from stat import ST_MODE, S_ISREG, S_ISDIR, S_ISLNK
from sysconfig import get_paths
from six import PY3

from wheel.paths import get_install_paths

from numpy.distutils.system_info import default_include_dirs, default_lib_dirs
from numpy.distutils.misc_util import get_numpy_include_dirs
from sysconfig import get_paths

sysconfig = distutils.sysconfig.get_config_vars()

# FIXME: remove redundnat code
src_ext = 'pyx'
use_cython = True


def get_wheel_header_paths():
    return [os.path.dirname(get_install_paths('dummy')['headers'])]


py3c_include_dirs = get_wheel_header_paths() + \
    [get_paths()['include']] + \
    [get_paths()['platinclude']]
numpy_include_dirs = get_numpy_include_dirs()
netcdf_include_dirs = default_include_dirs


include_dirs = ['Include'] + py3c_include_dirs + numpy_include_dirs + \
    netcdf_include_dirs

# num_package will always be numpy
num_package = "NumPy"
compile_args = ["-DNUMPY=1"]

headers = glob(os.path.join("Include", "MMTK", "*.h"))
headers.extend(glob(os.path.join("Include", "MMTK", "*.px[di]")))

paths = [
    os.path.join('MMTK', 'ForceFields'),
    os.path.join('MMTK', 'ForceFields', 'Amber'),
    os.path.join('MMTK', 'Database', 'Atoms'),
    os.path.join('MMTK', 'Database', 'Groups'),
    os.path.join('MMTK', 'Database', 'Molecules'),
    os.path.join('MMTK', 'Database', 'Complexes'),
    os.path.join('MMTK', 'Database', 'Proteins'),
    os.path.join('MMTK', 'Database', 'PDB'),
    os.path.join('MMTK', 'Tools', 'TrajectoryViewer')
]

data_files = []

for dir in paths:
    files = []
    for f in glob(os.path.join(dir, '*')):
        if f[-3:] != '.py' and f[-4:-1] != '.py' and os.path.isfile(f):
            files.append(f)
    data_files.append((dir, files))


#################################################################
# Check various compiler/library properties

libraries = []
if 'LIBM' in sysconfig and sysconfig['LIBM'] != '':
    libraries.append('m')

macros = []
try:
    from Scientific.MPI import world
except ImportError:
    world = None
if world is not None:
    if PY3:
        is_type = type(world) == type
    else:
        is_type = type(world) == types.InstanceType
    if is_type:
        world = None
if world is not None:
    macros.append(('WITH_MPI', None))

if sys.platform != 'win32':
    if hasattr(ctypes.CDLL(ctypes.util.find_library('m')), 'erfc'):
        macros.append(('LIBM_HAS_ERFC', None))

if sys.platform != 'win32':
    if ctypes.sizeof(ctypes.c_long) == 8:
        macros.append(('_LONG64_', None))

if sys.version_info[0] == 2 and sys.version_info[1] >= 2:
    macros.append(('EXTENDED_TYPES', None))

#################################################################
# System-specific optimization options

low_opt = []
if sys.platform != 'win32' and 'gcc' in sysconfig['CC']:
    low_opt = ['-O0']
low_opt.append('-g')

high_opt = []
if sys.platform[:5] == 'linux' and 'gcc' in sysconfig['CC']:
    high_opt = [
        '-O3', '-ffast-math', '-fomit-frame-pointer', '-fkeep-inline-functions'
    ]
if sys.platform == 'darwin' and 'gcc' in sysconfig['CC']:
    high_opt = [
        '-O3', '-ffast-math', '-fomit-frame-pointer', '-fkeep-inline-functions'
    ]
if sys.platform == 'aix4':
    high_opt = ['-O4']
if sys.platform == 'odf1V4':
    high_opt = ['-O2', '-fp_reorder', '-ansi_alias', '-ansi_args']

high_opt.append('-g')

#################################################################
# Extensions

# ext_pkg = 'MMTK'
ext_pkg = ''

extensions = [
    Extension(
        '%s.MMTK_DCD' % ext_pkg, ['Src/MMTK_DCD.c', 'Src/ReadDCD.c'],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_deformation' % ext_pkg, ['Src/MMTK_deformation.c'],
        extra_compile_args=compile_args + high_opt,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_dynamics' % ext_pkg, ['Src/MMTK_dynamics.c'],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_minimization' % ext_pkg, ['Src/MMTK_minimization.c'],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_surface' % ext_pkg, ['Src/MMTK_surface.c'],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_trajectory' % ext_pkg, ['Src/MMTK_trajectory.c'],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_universe' % ext_pkg, ['Src/MMTK_universe.c'],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_forcefield' % ext_pkg, [
            'Src/MMTK_forcefield.c', 'Src/bonded.c', 'Src/nonbonded.c',
            'Src/ewald.c', 'Src/sparsefc.c'
        ],
        extra_compile_args=compile_args + high_opt,
        include_dirs=include_dirs + ['Src'],
        define_macros=[('SERIAL', None), ('VIRIAL', None),
                       ('MACROSCOPIC', None)] + macros,
        libraries=libraries),
    Extension(
        '%s.MMTK_energy_term' % ext_pkg, ['Src/MMTK_energy_term.%s' % src_ext],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_restraints' % ext_pkg, ['Src/MMTK_restraints.%s' % src_ext],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_trajectory_action' % ext_pkg,
        ['Src/MMTK_trajectory_action.%s' % src_ext],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_trajectory_generator' % ext_pkg,
        ['Src/MMTK_trajectory_generator.%s' % src_ext],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
    Extension(
        '%s.MMTK_state_accessor' % ext_pkg,
        ['Src/MMTK_state_accessor.%s' % src_ext],
        extra_compile_args=compile_args,
        include_dirs=include_dirs,
        libraries=libraries,
        define_macros=macros),
]

if use_cython:
    extensions = cythonize(extensions, include_path=include_dirs)

#################################################################

setup(
    name=package_name,
    version="0.1",
    description="Molecular Modelling Toolkit",
    long_description="""
The Molecular Modelling Toolkit (MMTK) is an Open Source program
library for molecular simulation applications. It provides the most
common methods in molecular simulations (molecular dynamics, energy
minimization, normal mode analysis) and several force fields used for
biomolecules (Amber 94, Amber 99, several elastic network
models). MMTK also serves as a code basis that can be easily extended
and modified to deal with non-standard situations in molecular
simulations.
""",
    author="Konrad Hinsen",
    author_email="hinsen@cnrs-orleans.fr",
    url="http://dirac.cnrs-orleans.fr/MMTK/",
    license="CeCILL-C",
    packages=find_packages(),
    headers=headers,
    ext_modules=extensions,
    data_files=data_files,
    scripts=['tviewer'],
    command_options={'build_sphinx': {
        'source_dir': ('setup.py', 'Doc')
    }},
    use_2to3=True,
    install_requires=[
        'ScientificPython',
        'six',
        'scipy',
    ]
)
