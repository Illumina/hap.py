#!/usr/bin/env python
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt
#
# ONLY USE THIS FILE IF YOU WANT TO INSTALL HAP.PY GLOBALLY
# Otherwise, read README.md and have a look at install.py

import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).closeread()


setup(name="hap.py",
      version="0.1.3",
      author="Peter Krusche",
      author_email="pkrusche@illumina.com",
      description="Tools for generating and comparing haplotypes from VCF files.",
      license="BSD",
      keywords="VCF comparison",
      url="http://github.com/sequencing/hap.py",
      packages=['an_example_pypi_project', 'tests'],
      long_description=read('README.md'),
      )
