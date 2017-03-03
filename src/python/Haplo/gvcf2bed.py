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
# 10/02/2017
#
# Run gvcf2bed

import os
import tempfile
import subprocess
import copy
import json
import logging
import Tools
import pipes


def gvcf2bed(vcf,
             ref,
             regions=None,
             scratch_prefix=None):
    """ Run gvcf2bed and return temporary region bed file in temp folder """

    tf = tempfile.NamedTemporaryFile(dir=scratch_prefix, suffix=".bed")
    tf.close()
    cmdline = "gvcf2bed %s -r %s -o %s" % (pipes.quote(vcf), pipes.quote(ref), tf.name)
    if regions:
        cmdline += " -T %s" % pipes.quote(regions)
    logging.info("Running gvcf2bed: '%s'" % cmdline)
    subprocess.check_call(cmdline, shell=True)
    return tf.name

