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
# 01/04/2015
#
# Process raw counts coming out of quantify

import os
import tempfile
import subprocess
import copy
import json
import logging
import Tools
import pipes

from Tools.bcftools import runBcftools


def _locations_tmp_bed_file(locations):
    """ turn a list of locations into a bed file """
    if type(locations) is str:
        locations = locations.split(",")
    if type(locations) is not list:
        raise Exception("Invalid list of locations (must be str or list): %s" % str(locations))

    llocations = []

    for l in locations:
        xchr, _, _pos = l.partition(":")
        start, _, end = _pos.partition("-")
        if not xchr:
            raise Exception("Invalid chromosome name in %s" % str(l))
        try:
            start = int(start)
        except:
            start = 0

        try:
            end = int(end)
        except:
            end = 2 ** 31 - 1

        llocations.append((xchr, start, end))

    locations = sorted(llocations)

    tf = tempfile.NamedTemporaryFile(delete=False)
    for xchr, start, end in locations:
        print >> tf, "%s\t%i\t%i" % (xchr, start - 1, end)
    tf.close()

    return tf.name


def run_quantify(filename,
                 output_file=None, write_vcf=False, regions=None,
                 reference=Tools.defaultReference(),
                 locations=None, threads=1,
                 output_vtc=False,
                 output_rocs=False,
                 qtype=None,
                 roc_file=None,
                 roc_val=None,
                 roc_header=None,
                 roc_filter=None,
                 roc_delta=None,
                 roc_regions=None,
                 clean_info=True,
                 strat_fixchr=False):
    """Run quantify and return parsed JSON

    :param filename: the VCF file name
    :param output_file: output file name (if None, will use a temp file)
    :param write_vcf: write annotated VCF (give filename)
    :type write_vcf: str
    :param regions: dictionary of stratification region names and file names
    :param reference: reference fasta path
    :param locations: a location to use
    :param output_vtc: enable / disable the VTC field
    :param output_rocs: enable / disable output of ROCs by QQ level
    :param roc_file: filename for a TSV file with ROC observations
    :param roc_val: field to use for ROC QQ
    :param roc_header: name of ROC value for tables
    :param roc_filter: ROC filtering settings
    :param roc_delta: ROC minimum spacing between levels
    :param roc_regions: List of regions to output full ROCs for
    :param clean_info: remove unused INFO fields
    :param strat_fixchr: fix chr naming in stratification regions
    :returns: parsed counts JSON
    """

    if not output_file:
        output_file = tempfile.NamedTemporaryFile().name

    run_str = "quantify %s -o %s" % (
            pipes.quote(filename),
            pipes.quote(output_file))
    run_str += " -r %s" % pipes.quote(reference)
    run_str += " --threads %i" % threads

    if output_vtc:
        run_str += " --output-vtc 1"
    else:
        run_str += " --output-vtc 0"

    if output_rocs:
        run_str += " --output-rocs 1"
    else:
        run_str += " --output-rocs 0"

    if qtype:
        run_str += " --type %s" % qtype

    if roc_file:
        run_str += " --output-roc %s" % pipes.quote(roc_file)

    if roc_val:
        run_str += " --qq %s" % pipes.quote(roc_val)
        if roc_header != roc_val:
            # for xcmp, we extract the QQ value into the IQQ INFO field
            # we pass the original name along here
            run_str += " --qq-header %s" % pipes.quote(roc_header)

    if roc_filter:
        run_str += " --roc-filter '%s'" % pipes.quote(roc_filter)

    if roc_delta:
        run_str += " --roc-delta %f" % roc_delta

    if clean_info:
        run_str += " --clean-info 1"
    else:
        run_str += " --clean-info 0"

    if strat_fixchr:
        run_str += " --fix-chr-regions 1"
    else:
        run_str += " --fix-chr-regions 0"

    if write_vcf:
        if not write_vcf.endswith(".vcf.gz") and not write_vcf.endswith(".bcf"):
            write_vcf += ".vcf.gz"
        run_str += " -v %s" % pipes.quote(write_vcf)

    if regions:
        for k, v in regions.iteritems():
            run_str += " -R '%s:%s'" % (k, v)

    if roc_regions:
        for r in roc_regions:
            run_str += " --roc-regions '%s'" % r

    location_file = None
    if locations:
        location_file = _locations_tmp_bed_file(locations)
        run_str += " --only '%s'" % location_file

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stdout",
                                      suffix=".log")

    logging.info("Running '%s'" % run_str)

    try:
        subprocess.check_call(run_str, shell=True, stdout=tfo, stderr=tfe)
    except:
        tfo.close()
        tfe.close()
        with open(tfo.name) as f:
            for l in f:
                logging.error("[stdout] " + l.replace("\n", ""))
        os.unlink(tfo.name)
        with open(tfe.name) as f:
            for l in f:
                logging.error("[stderr] " + l.replace("\n", ""))
        os.unlink(tfe.name)
        if location_file:
            os.unlink(location_file)
        raise

    tfo.close()
    tfe.close()
    with open(tfo.name) as f:
        for l in f:
            logging.info("[stdout] " + l.replace("\n", ""))
    os.unlink(tfo.name)
    with open(tfe.name) as f:
        for l in f:
            logging.info("[stderr] " + l.replace("\n", ""))
    os.unlink(tfe.name)
    if location_file:
        os.unlink(location_file)

    if write_vcf and write_vcf.endswith(".bcf"):
        runBcftools("index", write_vcf)
    elif write_vcf:
        to_run = "tabix -p vcf %s" % pipes.quote(write_vcf)
        logging.info("Running '%s'" % to_run)
        subprocess.check_call(to_run, shell=True)


