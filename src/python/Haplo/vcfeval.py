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

import os
import logging
import subprocess
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

    del_vcf1 = False
    del_vcf2 = False

    try:
        if not os.path.exists(vcf1 + ".tbi"):
            vcf1 = Tools.bcftools.makeIndex(vcf1)
            del_vcf1 = True
        if not os.path.exists(vcf2 + ".tbi"):
            vcf2 = Tools.bcftools.makeIndex(vcf2)
            del_vcf2 = True

        runme = "%s vcfeval -b %s -c %s -t %s -o %s -T %i --baseline-tp" % (
            args.engine_vcfeval,
            vcf1.replace(" ", "\\ "),
            vcf2.replace(" ", "\\ "),
            args.engine_vcfeval_template.replace(" ", "\\ "),
            vtf.name,
            args.threads)
        logging.info(runme)
        po = subprocess.Popen(runme,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        o, e = po.communicate()

        po.wait()

        rc = po.returncode

        if rc != 0:
            raise Exception("Error running rtg tools / vcfeval. Return code was %i, output: %s / %s \n" % (rc, o, e))
        else:
            logging.info("vcfeval output: \n%s\n / \n%s\n" % (o, e))

        runme = "postvcfeval %s %s -r %s" % (vtf.name, target, args.ref)
        logging.info(runme)
        po = subprocess.Popen(runme,
                              shell=True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        o, e = po.communicate()

        po.wait()

        rc = po.returncode

        if rc != 0:
            raise Exception("Error running postvcfeval. Return code was %i, output: %s / %s \n" % (rc, o, e))
        else:
            logging.info("postvcfeval output: %s / %s" % (o, e))

        Tools.bcftools.runBcftools("index", "-t", target)
    finally:
        if del_vcf1:
            try:
                os.unlink(vcf1)
                os.unlink(vcf1 + ".tbi")
            except:
                pass
        if del_vcf2:
            try:
                os.unlink(vcf2)
                os.unlink(vcf2 + ".tbi")
            except:
                pass
        # remove temp path
        try:
            shutil.rmtree(vtf.name)
        except:
            pass

    elapsed = time.time() - starttime
    logging.info("vcfeval for %s vs. %s -- time taken %.2f" % (vcf1, vcf2, elapsed))

    return [target, target + ".tbi"]
