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
import pipes

from Tools.parallel import runParallel, getPool
from Tools.bcftools import runBcftools, concatenateParts
from Tools.vcfextract import extractHeadersJSON


def preprocessWrapper(file_and_location, args):
    starttime = time.time()
    filename, location_str = file_and_location
    if args["bcf"]:
        int_suffix = "bcf"
    else:
        int_suffix = "vcf.gz"

    tf = tempfile.NamedTemporaryFile(delete=False,
                                     prefix="input.%s" % location_str,
                                     suffix=".prep." + int_suffix)
    tf.close()

    to_run = "preprocess %s:* %s-o %s -V %i -L %i -r %s" % \
             (pipes.quote(filename),
              ("-l %s " % pipes.quote(location_str)) if location_str else "",
              tf.name,
              args["decompose"],
              args["leftshift"],
              pipes.quote(args["reference"]))

    if args["haploid_x"]:
        to_run += " --haploid-x 1"

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stdout",
                                      suffix=".log")
    finished = False
    try:
        logging.info("Running '%s'" % to_run)
        subprocess.check_call(to_run, shell=True, stdout=tfo, stderr=tfe)
        finished = True
    finally:
        if finished:
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
        else:
            logging.error("Preprocess command %s failed. Outputs are here %s / %s" % (to_run, tfo.name, tfe.name))
            with open(tfo.name) as f:
                for l in f:
                    logging.error(l.replace("\n", ""))
            with open(tfe.name) as f:
                for l in f:
                    logging.error(l.replace("\n", ""))

    elapsed = time.time() - starttime
    logging.info("preprocess for %s -- time taken %.2f" % (location_str, elapsed))
    runBcftools("index", tf.name)
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
                 (pipes.quote(bargs["vcf"]),
                  pipes.quote(location_str),
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

    finally:
        elapsed = time.time() - starttime
        logging.info("blocksplit for %s -- time taken %.2f" % (location_str, elapsed))
        os.unlink(tf.name)
    return result

def partialCredit(vcfname,
                  outputname,
                  reference,
                  locations,
                  threads=1,
                  window=10000,
                  leftshift=True,
                  decompose=True,
                  haploid_x=False):
    """ Partial-credit-process a VCF file according to our args """

    pool = getPool(int(threads))
    if threads > 1:
        logging.info("Partial credit processing uses %i parallel processes." % threads)

        if not locations:
            h = extractHeadersJSON(vcfname)
            if not h["tabix"]["chromosomes"]:
                logging.warn("Empty input or not tabix indexed")
                if outputname.endswith(".bcf"):
                    runBcftools("view", "-O", "b", "-o", outputname, vcfname)
                    runBcftools("index", outputname)
                else:
                    runBcftools("view", "-O", "z", "-o", outputname, vcfname)
                    runBcftools("index", "-t", outputname)
                # just return the same file
                return
            locations = h["tabix"]["chromosomes"]
        elif type(locations) is str or type(locations) is unicode:
            locations = locations.split(",")

        # use blocksplit to subdivide input
        res = runParallel(pool,
                          blocksplitWrapper,
                          locations,
                          {"vcf": vcfname,
                           "dist": window,
                           "pieces": min(40, threads*4)})

        if None in res:
            raise Exception("One of the blocksplit processes failed.")

        locations = list(itertools.chain.from_iterable(res))
        if not len(locations):
            logging.warn("Blocksplit returned no blocks. This can happen when "
                         "an input contains no valid variants.")
            locations = [""]
    else:
        locations = [""]

    res = []
    try:
        res = runParallel(pool,
                          preprocessWrapper,
                          itertools.izip(itertools.repeat(vcfname), locations),
                          {"reference": reference,
                           "decompose": decompose,
                           "leftshift": leftshift,
                           "haploid_x": haploid_x,
                           "bcf": outputname.endswith(".bcf")})

        if None in res:
            raise Exception("One of the preprocess jobs failed")
        if not res:
            raise Exception("No blocks were processed. List of locations: %s" % str(list(locations)))

        concatenateParts(outputname, *res)
        if outputname.endswith(".vcf.gz"):
            runBcftools("index", "-f", "-t", outputname)
        else:  # use bcf
            runBcftools("index", "-f", outputname)
    finally:
        for r in res:
            try:
                os.unlink(r)
            except:
                pass
            try:
                os.unlink(r + ".tbi")
            except:
                pass
            try:
                os.unlink(r + ".csi")
            except:
                pass
