#!/usr/bin/env python
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
# 17/2/2015
#
# Somatic VCF feature extraction
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import sys
import os
import argparse
import logging
import traceback
import tempfile
import shutil
import pandas

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools.bcftools import runBcftools, preprocessVCF
from Tools.bamstats import bamStats

import Somatic


def main():
    parser = argparse.ArgumentParser("Somatic VCF Feature Extraction")

    parser.add_argument("input", help="Input VCF file")

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name. Output will be in CSV format")

    parser.add_argument("-l", "--location", dest="location", default="",
                        help="Location for bcftools view (e.g. chr1)")

    parser.add_argument("-R", "--restrict-regions", dest="regions_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (sparse) regions (using -R in bcftools).")

    parser.add_argument("-T", "--target-regions", dest="targets_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (dense) regions (using -T in bcftools).")

    parser.add_argument("-P", "--include-nonpass", dest="inc_nonpass", action="store_true", default=False,
                        help="Use to include failing variants in comparison.")

    parser.add_argument("--feature-table", dest="features", default="generic",
                        help="Select a feature table to output. Options are: %s" % str(Somatic.FeatureSet.sets.keys()))

    parser.add_argument("--feature-label", dest="label", default=None,
                        help="We will output a lable column, this value will go in there -- default is "
                             "the input filename.")

    parser.add_argument("--bam", dest="bams", default=[], action="append",
                        help="pass one or more BAM files for feature table extraction")

    parser.add_argument("-r", "--reference", dest="ref", default=Tools.defaultReference(),
                        help="Specify a reference file for normalization.")

    parser.add_argument("--normalize", dest="normalize", default=False, action="store_true",
                        help="Enable running of bcftools norm on the input file.")

    parser.add_argument("--fix-chr", dest="fixchr", default=False, action="store_true",
                        help="Replace numeric chromosome names in the query by chr*-type names")

    args = parser.parse_args()

    scratch = tempfile.mkdtemp()

    try:
        logging.info("Scratch path is %s" % scratch)

        if not args.label:
            args.label = os.path.basename(args.input)

        bams = []
        md = None
        for x in args.bams:
            bams.append(bamStats(x))

        if bams:
            bres = pandas.concat(bams).groupby("CHROM").mean()
            md = {}
            for x in bres.index:
                logging.info("Mean coverage on %s is %f" % (x, bres.loc[x]["COVERAGE"]))
                md[x] = float(bres.loc[x]["COVERAGE"])*3.0

        nqpath = os.path.join(scratch, "normalized_query.vcf.gz")

        logging.info("Preprocessing input...")
        preprocessVCF(args.input, nqpath, args.location,
                      not args.inc_nonpass,  # pass_only
                      args.fixchr,  # chrprefix
                      args.normalize,  # norm,
                      args.regions_bedfile,
                      args.targets_bedfile,
                      args.ref)

        runBcftools("index", nqpath)

        logging.info("Extracting features...")
        fset = Somatic.FeatureSet.make(args.features)
        fset.setChrDepths(md)
        featuretable = fset.collect(nqpath, args.label)

        if not args.output.endswith(".csv"):
            args.output += ".csv"
        logging.info("Saving feature table %s..." % args.output)
        featuretable.to_csv(args.output)

    finally:
        logging.info("Deleting scratch folder %s " % scratch)
        shutil.rmtree(scratch)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
