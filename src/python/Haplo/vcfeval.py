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
import logging
import Tools.bcftools

import Haplo.version


def findVCFEval():
    """ Return default version of rtgtools if hap.py was built with
        rtgtools included.
    """
    if Haplo.version.has_vcfeval:
        base = os.path.join(os.path.dirname(__file__),
                            "..",   # python27
                            "..",   # lib
                            "..",   # hap.py-base
                            "libexec",
                            "rtg-tools-install")
        # prefer wrapper when it's there
        bfile = os.path.join(base, "rtg-wrapper.sh")
        bfile2 = os.path.join(base, "rtg")
        if os.path.isfile(bfile) and os.access(bfile, os.X_OK):
            return bfile
        elif os.path.isfile(bfile2) and os.access(bfile2, os.X_OK):
            return bfile2
        else:
            logging.warn(("Could not find our included version of rtg-tools at %s. " % base) +
                         "To use vcfeval for comparison, you might have to specify "
                         "its location on the command line.")
            return "rtg"
    else:
        # default: return
        return "rtg"


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
    del_sdf = False

    try:
        if not os.path.exists(vcf1 + ".tbi"):
            vcf1 = Tools.bcftools.makeIndex(vcf1)
            del_vcf1 = True
        if not os.path.exists(vcf2 + ".tbi"):
            vcf2 = Tools.bcftools.makeIndex(vcf2)
            del_vcf2 = True

        if not args.engine_vcfeval_template or not os.path.exists(args.engine_vcfeval_template):
            logging.warn("Creating template for vcfeval. "
                         "You can speed this up in the future by doing this in advance.")
            del_sdf = True
            stf = tempfile.NamedTemporaryFile(dir=args.scratch_prefix,
                                              prefix="vcfeval.sdf",
                                              suffix=".dir")
            stf.close()
            args.engine_vcfeval_template = stf.name
            runme = "%s format -o %s %s" % (
                args.engine_vcfeval,
                args.engine_vcfeval_template.replace(" ", "\\ "),
                args.ref.replace(" ", "\\ "))
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
            elif o.strip() or e.strip():
                logging.info("vcfeval output: \n%s\n / \n%s\n" % (o, e))

        runme = "%s vcfeval -b %s -c %s -t %s -o %s -T %i -m ga4gh --ref-overlap" % (
            args.engine_vcfeval,
            vcf1.replace(" ", "\\ "),
            vcf2.replace(" ", "\\ "),
            args.engine_vcfeval_template.replace(" ", "\\ "),
            vtf.name,
            args.threads)

        if not args.pass_only:
            runme += " --all-records"

        if args.roc:
            runme += " -f %s" % args.roc

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
        elif o.strip() or e.strip():
            logging.info("vcfeval output: \n%s\n / \n%s\n" % (o, e))

        # in GA4GH mode, this is what vcfeval should output
        shutil.copy(os.path.join(vtf.name, "output.vcf.gz"), target)
        shutil.copy(os.path.join(vtf.name, "output.vcf.gz.tbi"), target + ".tbi")
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
        if del_sdf:
            try:
                shutil.rmtree(args.engine_vcfeval_template)
            except:
                pass

    elapsed = time.time() - starttime
    logging.info("vcfeval for %s vs. %s -- time taken %.2f" % (vcf1, vcf2, elapsed))

    if os.path.exists(target) and os.path.exists(target + ".tbi"):
        return [target, target + ".tbi"]
    else:
        return None
