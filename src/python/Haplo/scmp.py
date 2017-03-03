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
import pipes

import Haplo.version

from Tools.bcftools import runBcftools
from Tools import LoggingWriter

def runSCmp(vcf1, vcf2, target, args):
    """ Runs scmp, which outputs a file quantify can produce counts on
    vcf1 and vcf2 must be indexed and only contain a single sample column.
    """

    try:
        # change GTs so we can compare them
        vargs = ["merge", "--force-samples", vcf1, vcf2,
                 "|", "scmp", "-", "-r", args.ref,
                 "--threads", str(args.threads),
                 "-o", target]

        if args.roc:
            vargs += ["--q", args.roc]

        runBcftools(*vargs)

        if target.endswith(".vcf.gz"):
            runBcftools("index", "-t", target)
            return [target, target + ".tbi"]
        else:
            runBcftools("index", target)
            return [target, target + ".csi"]
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


