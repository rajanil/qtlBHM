
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import sys

# setup bed parser
ext_modules = [Extension("bhmodel", ["bhmodel.pyx"])]
ext_modules = cythonize(ext_modules)

setup(
    name = 'Bayesian hierarchical model',
    author = 'Anil Raj',
    version = '1.0',
    cmdclass = {'build_ext': build_ext},
    include_dirs=[numpy.get_include(), '.'],
    ext_modules = ext_modules
)
