#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys

if sys.version_info < (2, 6, 0):
    sys.stderr.write("Requires Python 2.6 or newer.\n")
    sys.exit(-1)

from distutils.core import setup

def get_version():
    VERSIONFILE="../VERSION"
    return open(VERSIONFILE, "r").read().strip()

setup(
    name='GenomeTools',
    version=get_version(),
    description= 'Python bindings for GenomeTools',
    author='Sascha Steinbiss',
    author_email='steinbiss@zbh.uni-hamburg.de',
    url= 'http://www.genometools.org ',
    packages=['gt', 'gt.core',
             'gt.annotationsketch', 'gt.extended']
)
