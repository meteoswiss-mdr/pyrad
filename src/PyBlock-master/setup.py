"""
PyBlock
"""

import os
import sys
from setuptools import setup, find_packages
from distutils.sysconfig import get_python_lib

# - Pull the header into a variable
doclines = __doc__.split("\n")

VERSION = '1.3'

DATADIR = os.sep.join([os.path.dirname(__file__), 'data'])

# - Set variables for setup
PACKAGES = ['pyblock']
package_dir = {'': 'pyblock'}

# - Run setup
setup(
      name='pyblock',
      version=VERSION,
      url='http://github.com/tjlang',
      author='Timothy Lang',
      author_email='timothy.j.lang@nasa.gov',
      description=doclines[0],
      packages=PACKAGES,
      package_data={'pyblock': ['data/*.txt']},
      include_package_data=True,
      classifiers=["""
        Development Status :: Beta,
        Programming Language :: Python",
        Topic :: Scientific/Engineering
        Topic :: Scientific/Engineering :: Atmospheric Science
        Operating System :: Unix
        Operating System :: POSIX :: Linux
        Operating System :: MacOS
        """],
      long_description="""Blockage Correction Tools""",
      install_requires=['numpy'],
      )
