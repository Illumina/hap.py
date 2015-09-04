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
# 3/9/2015
#
# Diploid VCF File Comparison
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import os
import logging
import tempfile
import time
import subprocess


def xcmpWrapper(location_str, args):
    """ Haplotype block comparison wrapper function
    """
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     dir=args.scratch_prefix,
                                     prefix="result.%s" % location_str,
                                     suffix=".vcf.gz")
    tf.close()

    if args.write_bed:
        tf2 = tempfile.NamedTemporaryFile(delete=False,
                                          dir=args.scratch_prefix,
                                          prefix="result.blocks.%s" % location_str,
                                          suffix=".bed")
        tf2.close()
        bname = "-e %s" % tf2.name
    else:
        bname = ""

    to_run = "xcmp %s %s -l %s -o %s %s -r %s -f %i --apply-filters-truth %i -n %i -V %i --leftshift %i --expand-hapblocks %i " \
             "--window %i --compare-raw %i --no-hapcmp %i --roc-vals %i" % \
             (args.vcf1.replace(" ", "\\ "),
              args.vcf2.replace(" ", "\\ "),
              location_str,
              tf.name,
              bname,
              args.ref,
              0 if args.usefiltered else 1,
              0 if args.usefiltered_truth else 1,
              args.max_enum,
              1 if args.int_preprocessing else 0,
              1 if args.int_preprocessing_ls else 0,
              args.hb_expand,
              args.window,
              1 if args.int_match_raw else 0,
              1 if args.no_hc else 0,
              1 if args.roc else 0
              )

    # regions / targets already have been taken care of in blocksplit / preprocessing

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      dir=args.scratch_prefix,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      dir=args.scratch_prefix,
                                      prefix="stdout",
                                      suffix=".log")

    try:
        logging.info("Running '%s'" % to_run)
        subprocess.check_call(to_run, shell=True, stdout=tfo, stderr=tfe)
    finally:
        tfo.close()
        tfe.close()
        with open(tfo.name) as f:
            for l in f:
                logging.info(l.replace("\n", ""))
        os.unlink(tfo.name)
        with open(tfe.name) as f:
            for l in f:
                logging.warn(l.replace("\n", ""))
        os.unlink(tfe.name)

    elapsed = time.time() - starttime
    logging.info("xcmp for chunk %s -- time taken %.2f" % (location_str, elapsed))

    if bname == "":
        bname = None
    else:
        bname = bname[3:]

    return tf.name, bname
