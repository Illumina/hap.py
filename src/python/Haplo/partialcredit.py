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
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import os
import logging
import subprocess
import tempfile
import time
import itertools
import multiprocessing

from Tools.parallel import runParallel
from Tools.bcftools import runBcftools


def preprocessWrapper(file_and_location, args):
    starttime = time.time()
    filename, location_str = file_and_location
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     dir=args.scratch_prefix,
                                     prefix="input.%s" % location_str,
                                     suffix=".prep.vcf.gz")
    tf.close()

    to_run = "preprocess %s:* %s-o %s -V %i -L %i -r %s" % \
             (filename.replace(" ", "\\ "),
              ("-l %s " % location_str) if location_str else "",
              tf.name,
              1 if args.int_preprocessing else 0,
              1 if args.int_preprocessing_ls else 0,
              args.ref)

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
    logging.info("preprocess for %s -- time taken %.2f" % (location_str, elapsed))
    runBcftools("index", "-t", tf.name)
    return tf.name


def blocksplitWrapper(location_str, bargs):
    """ Blocksplit for partial credit preprocessing """
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     prefix="result.%s" % location_str,
                                     suffix=".chunks.bed")
    result = None
    try:
        tf.close()

        to_run = "blocksplit %s -l %s -o %s --window %i --nblocks %i -f 0" % \
                 (bargs["vcf"],
                  location_str,
                  tf.name,
                  bargs["dist"],
                  bargs["pieces"])

        tfe = tempfile.NamedTemporaryFile(delete=False,
                                          prefix="stderr",
                                          suffix=".log")
        tfo = tempfile.NamedTemporaryFile(delete=False,
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

        r = []
        with open(tf.name) as f:
            for l in f:
                ll = l.strip().split("\t", 3)
                if len(ll) < 3:
                    continue
                xchr = ll[0]
                start = int(ll[1]) + 1
                end = int(ll[2])
                r.append("%s:%i-%i" % (xchr, start, end))
        result = r

        elapsed = time.time() - starttime
    finally:
        logging.info("blocksplit for %s -- time taken %.2f" % (location_str, elapsed))
        os.unlink(tf.name)
    return result

def partialCredit(vcfname, outputname, args):
    """ Partial-credit-process a VCF file according to our args """

    locations = [""]
    if args.threads > 1:
        logging.info("Partial credit processing uses %i parallel processes." % args.threads)
        pool = multiprocessing.Pool(int(args.threads))

        # use blocksplit to subdivide input
        res = runParallel(pool,
                          blocksplitWrapper,
                          args.locations,
                          {"vcf": vcfname, "dist": 20000, "pieces": min(40, args.threads*4)})

        if None in res:
            raise Exception("One of the blocksplit processes failed.")

        locations = itertools.chain.from_iterable(res)
    else:
        pool = None

    res = []
    try:
        res = runParallel(pool,
                          preprocessWrapper,
                          itertools.izip(itertools.repeat(vcfname), locations),
                          args)

        if None in res:
            raise Exception("One of the preprocess jobs failed")

        cmd = ["concat", "-a", "-o", outputname, "-O", "z"] + res
        runBcftools(*cmd)
        runBcftools("index", "-t", outputname)
    finally:
        for r in res:
            try:
                os.unlink(r)
                os.unlink(r + ".tbi")
            except:
                pass
