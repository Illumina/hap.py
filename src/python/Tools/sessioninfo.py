# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# 31/05/2017
#
# Collect session and run information
#

import os
import sys
import time
import platform
import copy

import Tools


def sessionInfo():
    """ Return a dictionary with session and run information
    """

    version = "%s" % Tools.version

    result = {'name': os.path.basename(sys.argv[0]),
              'timestamp': time.strftime("%a %b %d %X %Y"),
              'version': version,
              'runInfo': [{"key": "commandline", "value": " ".join(sys.argv)}],
              'uname': " / ".join(platform.uname()),
              'dist': " / ".join(platform.dist()),
              'mac_ver': " / ".join([platform.mac_ver()[0], platform.mac_ver()[2]]),
              'python_implementation': platform.python_implementation(),
              'python_version': platform.python_version(),
              'metadata': {
                  "required": {
                      "id": "haplotypes",
                      'version': version,
                      "module": "%s" % os.path.basename(sys.argv[0]),
                      "description": "%s generated this JSON file via command line %s" % (
                          sys.argv[0], " ".join(sys.argv))}},
              'environment': {str(k): str(os.environ[k]) for k in os.environ.keys()}}

    result["python_prefix"] = sys.prefix
    if hasattr(sys, 'real_prefix'):
        result["python_virtualenv"] = True
        result["python_real_prefix"] = sys.real_prefix

    try:
        import psutil
        result["cpus"] = psutil.cpu_count()
        result["logical_cpus"] = psutil.cpu_count(True)
        result["cpu_freq"] = psutil.cpu_freq()
        result["memory"] = dict(psutil.virtual_memory().__dict__)
    except:
        pass

    try:
        import pip
        pip_packages = []
        for i in pip.get_installed_distributions(local_only=True):
            pip_packages.append(str(i))

        result["pip_packages"] = pip_packages
    except:
        pass

    return result

