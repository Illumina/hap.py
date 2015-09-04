#!/illumina/development/haplocompare/hc-virtualenv/bin/python
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
# 3/9/2014
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import sys
import os
import logging
import subprocess
import gzip
import tempfile
import time
import shutil

import Tools.bcftools

def runVCFEval(vcf1, vcf2, target, args):
    """ Run VCFEval and convert it's output to something quantify
        understands
    """
    starttime = time.time()

    vtf = tempfile.NamedTemporaryFile(dir=args.scratch_prefix,
                                      prefix="vcfeval.result",
                                      suffix=".dir")
    vtf.close()

    try:
        runme = "%s vcfeval -b %s -c %s -t %s -o %s -T %i --baseline-tp" % (
            args.engine_vcfeval,
            vcf1.replace(" ", "\\ "),
            vcf2.replace(" ", "\\ "),
            args.engine_vcfeval_template.replace(" ", "\\ "),
            vtf.name,
            args.threads)
        logging.info(runme)
        po = subprocess.Popen(runme,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        o, e = po.communicate()

        po.wait()

        rc = po.returncode

        if rc != 0:
            raise Exception("Error running rtg tools / vcfeval. Return code was %i, output: %s / %s \n" % (rc, o, e))
        else:
            logging.info("vcfeval output: %s / %s" % (o, e))

        runme = "postvcfeval -o %s %s" % (target, vtf.name)
        logging.info(runme)
        po = subprocess.Popen(runme,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        o, e = po.communicate()

        po.wait()

        rc = po.returncode

        if rc != 0:
            raise Exception("Error running postvcfeval. Return code was %i, output: %s / %s \n" % (rc, o, e))
        else:
            logging.info("postvcfeval output: %s / %s" % (o, e))
    finally:
        # remove temp path
        shutil.rmtree(vtf.name)

    elapsed = time.time() - starttime
    logging.info("vcfeval for %s vs. %s -- time taken %.2f" % (vcf1, vcf2, elapsed))
