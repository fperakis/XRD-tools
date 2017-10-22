#!/usr/bin/env python
#coding: utf8

"""
Setup script for XRD-tools.
"""

from glob import glob

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='xrd-tools',
      version='0.0.1',
      author="Fivos Perakis",
      author_email="f.perakis@fysik.su.se",
      url='https://github.com/fperakis/XRD-tools',
      description='Collection of tools to analyse x-ray diffraction data',
      packages=["xrdtools"],
      package_dir={'xrdtools': 'src'},
      #install_requires=['scikit-beam', 'scikit-image']
      #scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')],
      #test_suite="test"
      )
