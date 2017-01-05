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
import shutil
import logging
import traceback

import Haplo.version

from Tools.bcftools import runBcftools
from Tools import LoggingWriter

def runSCmp(vcf1, vcf2, target, args):
    """ Runs scmp, which outputs a file quantify can produce counts on

    """

    tf1 = tempfile.NamedTemporaryFile(suffix=".bcf")
    tf2 = tempfile.NamedTemporaryFile(suffix=".bcf")
    tf3 = tempfile.NamedTemporaryFile(suffix=".bcf")
    for f in [tf1, tf2, tf3]:
        f.close()
    try:
        # change GTs so we can compare them
        alleles_cmd = "alleles '%s' -o '%s'" % (vcf1, tf1.name)
        subprocess.check_call(alleles_cmd, shell=True)

        alleles_cmd = "alleles '%s' -o '%s'" % (vcf2, tf2.name)
        subprocess.check_call(alleles_cmd, shell=True)

        runBcftools("index", tf1.name)
        runBcftools("index", tf2.name)

        runBcftools("merge", "--force-samples", tf1.name, tf2.name, "-o", tf3.name, "-O", "b")
        runBcftools("index", tf3.name)

        scmp_cmd = "scmp '%s' -r '%s' --threads %i -o '%s'" % (tf3.name, args.ref, args.threads, target)
        if args.roc:
            scmp_cmd += " --q '%s'" % args.roc
        subprocess.check_call(scmp_cmd, shell=True)
        runBcftools("index", target)
    except Exception as e:
        logging.error("Exception when running scmp: %s" % str(e))
        logging.error('-'*60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-'*60)
        raise
    except BaseException as e:
        logging.error("Exception when running scmp: %s" % str(e))
        logging.error('-'*60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-'*60)
        raise
    finally:
        for f in [tf1, tf2, tf3]:
            try:
                os.unlink(f.name)
            except:
                pass

    return [target, target + ".csi"]

