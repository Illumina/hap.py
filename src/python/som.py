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
# 9/9/2014
#
# Somatic VCF Comparison and feature extraction
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
import numpy as np
np.seterr(all='ignore')
import gzip
import json
from collections import Counter

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools.bcftools import runBcftools, parseStats, preprocessVCF
from Tools.bamstats import bamStats
from Tools.bedintervaltree import BedIntervalTree
from Tools.roc import ROC
from Tools.ci import binomialCI
from Tools.metric import makeMetricsObject, dataframeToMetricsTable
from Tools.fastasize import fastaContigLengths, calculateLength
import Somatic


def parse_args():
    parser = argparse.ArgumentParser("Somatic Comparison")

    parser.add_argument("truth", help="Truth VCF file")
    parser.add_argument("query", help="Query VCF file")

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file prefix for statistics and feature table (when selected)")

    parser.add_argument("-l", "--location", dest="location", default="",
                        help="Location for bcftools view (e.g. chr1)")

    parser.add_argument("-R", "--restrict-regions", dest="regions_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (sparse) regions (using -R in bcftools).")

    parser.add_argument("-T", "--target-regions", dest="targets_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (dense) regions (using -T in bcftools).")

    parser.add_argument("-f", "--false-positives", dest="FP",
                        help="False-positive region bed file to distinguish UNK from FP")

    parser.add_argument("-a", "--ambiguous", dest="ambi", action='append',
                        help="Ambiguous region bed file(s) to distinguish from FP (e.g. variant only observed "
                             "in some replicates)")

    parser.add_argument("--ambi-fp", dest="ambi_fp", action='store_true', default=False,
                        help="Use FP calls from ambiguous region files also.")

    parser.add_argument("--no-ambi-fp", dest="ambi_fp", action='store_false',
                        help="Do not use FP calls from ambiguous region files also.")

    parser.add_argument("--count-unk", dest="count_unk", action='store_true', default=False,
                        help="Assume the truth set covers the whole genome and only count FPs in regions "
                             "specified by the truth VCF or ambiguous/false-positive regions.")

    parser.add_argument("--no-count-unk", dest="count_unk", action='store_false',
                        help="Do not use FP calls from ambiguous region files also.")

    parser.add_argument("-e", "--explain_ambiguous", dest="explain_ambiguous", required=False,
                        default=False, action="store_true",
                        help="print a table giving the number of ambiguous events per category")

    parser.add_argument("-r", "--reference", dest="ref", default=Tools.defaultReference(),
                        help="Specify a reference file.")

    parser.add_argument("--scratch-prefix", dest="scratch_prefix",
                        default=None,
                        help="Filename prefix for scratch report output.")

    parser.add_argument("--keep-scratch", dest="delete_scratch",
                        default=True, action="store_false",
                        help="Filename prefix for scratch report output.")

    parser.add_argument("--continue", dest="cont", default=False, action="store_true",
                        help="Continue from scratch space (i.e. use VCFs in there if they already exist).")

    parser.add_argument("-P", "--include-nonpass", dest="inc_nonpass", action="store_true", default=False,
                        help="Use to include failing variants in comparison.")

    parser.add_argument("--feature-table", dest="features", default=False, choices=Somatic.FeatureSet.sets.keys(),
                        help="Select a feature table to output.")

    parser.add_argument("--happy-stats", dest="happy_stats", default=False, action="store_true",
                        help="Generate summary.csv.")

    parser.add_argument("--bam", dest="bams", default=[], action="append",
                        help="pass one or more BAM files for feature table extraction")

    parser.add_argument("--normalize-truth", dest="normalize_truth", default=False, action="store_true",
                        help="Enable running of bcftools norm on the truth file.")

    parser.add_argument("--normalize-query", dest="normalize_query", default=False, action="store_true",
                        help="Enable running of bcftools norm on the query file.")

    parser.add_argument("-N", "--normalize-all", dest="normalize_all", default=False, action="store_true",
                        help="Enable running of bcftools norm on both truth and query file.")

    parser.add_argument("--fixchr-truth", dest="fixchr_truth", action="store_true", default=True,
                        help="Add chr prefix to truth file (default: true).")

    parser.add_argument("--fixchr-query", dest="fixchr_query", action="store_true", default=True,
                        help="Add chr prefix to query file (default: true).")

    parser.add_argument("--fix-chr-truth", dest="fixchr_truth", action="store_true", default=None,
                        help="Same as --fixchr-truth.")

    parser.add_argument("--fix-chr-query", dest="fixchr_query", action="store_true", default=None,
                        help="Same as --fixchr-query.")

    parser.add_argument("--no-fixchr-truth", dest="fixchr_truth", action="store_false", default=False,
                        help="Disable chr replacement for truth (default: false).")

    parser.add_argument("--no-fixchr-query", dest="fixchr_query", action="store_false", default=False,
                        help="Add chr prefix to query file (default: false).")

    parser.add_argument("--no-order-check", dest="disable_order_check", default=False, action="store_true",
                        help="Disable checking the order of TP features (dev feature).")

    parser.add_argument("--roc", dest="roc", default=None, choices=ROC.list(),
                        help="Create a ROC-style table. This is caller specific "
                             " - this will override the --feature-table switch!")

    parser.add_argument("--bin-afs", dest="af_strat", default=None, action="store_true",
                        help="Stratify into different AF buckets. This needs to have features available "
                             "for getting the AF both in truth and query variants.")
    parser.add_argument("--af-binsize", dest="af_strat_binsize", default=0.2,
                        help="Bin size for AF binning (should be < 1). Multiple bin sizes can be specified using a comma, "
                             "e.g. 0.1,0.2,0.5,0.2 will split at 0.1, 0.3, 0.8 and 1.0.")
    parser.add_argument("--af-truth", dest="af_strat_truth", default="I.T_ALT_RATE",
                        help="Feature name to use for retrieving AF for truth variants (TP and FN)")
    parser.add_argument("--af-query", dest="af_strat_query", default="T_AF",
                        help="Feature name to use for retrieving AF for query variants (FP/UNK/AMBI)")

    parser.add_argument("-FN", "--count-filtered-fn", dest="count_filtered_fn", action="store_true",
                        help="Count filtered vs. absent FN numbers. This requires the -P switch (to use all "
                             "variants) and either the --feature-table or --roc switch.")

    parser.add_argument("--fp-region-size", dest="fpr_size",
                        help="How to obtain the normalisation constant for FP rate. By default, this will use the FP region bed file size when using"
                             " --count-unk and the size of all reference contigs that overlap with the location specified in -l otherwise."
                             " This can be overridden with: 1) a number of nucleotides, or 2) \"auto\" to use the lengths of all contigs that have calls."
                             " The resulting value is used as fp.region.size.")

    parser.add_argument("--ci-level", dest="ci_level", default=0.95, type = float,
                        help="Confidence level for precision/recall confidence intervals (default: 0.95)")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    args = parser.parse_args()

    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

    try:
        if type(args.af_strat_binsize) is str:
            args.af_strat_binsize = map(float, args.af_strat_binsize.split(","))
        else:
            args.af_strat_binsize = map(float, [args.af_strat_binsize])

        if not args.af_strat_binsize:
            raise Exception("Bin size list is empty")
    except:
        logging.error("Failed to parse stratification bin size: %s" % str(args.af_strat_binsize))
        exit(1)

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        level=loglevel)

    if args.normalize_all:
        args.normalize_truth = True
        args.normalize_query = True

    if args.roc:
        args.roc = ROC.make(args.roc)
        args.features = args.roc.ftname
        if not args.inc_nonpass:
            logging.warn("When creating ROCs without the -P switch, the ROC data points will only "
                         "include filtered variants (i.e. they will normally end at the caller's "
                         "quality threshold).")

    if not (args.ci_level > 0.0 and args.ci_level < 1.0):
        raise Exception("Confidence interval level must be > 0.0 and < 1.0.")

    if args.af_strat and not args.features:
        raise Exception("To stratify by AFs, a feature table must be selected -- use this switch together "
                        "with --feature-table or --roc")

    if args.count_filtered_fn and (not args.inc_nonpass or not args.features):
        raise Exception("Counting filtered / unfiltered FNs only works when a feature table is selected, "
                        "and when using unfiltered variants. Specify -P --feature-table <...> or use "
                        "--roc to select a ROC type.")

    if args.happy_stats and (not args.inc_nonpass or not args.features):
        raise Exception("Obtaining hap.py stats only works when a feature table is selected, "
                        "and when using unfiltered variants. Specify -P --feature-table <...> or use "
                        "--roc to select a ROC type.")

    if args.scratch_prefix:
        scratch = os.path.abspath(args.scratch_prefix)
        args.delete_scratch = False
        Tools.mkdir_p(scratch)
    else:
        scratch = tempfile.mkdtemp()

    return args


def resolve_vtype(args):
    fvtype = args.features.split(".")[-1]
    if fvtype == "snv":
        vtype = "SNP"
    elif fvtype == "indel":
        vtype = "INDEL"
    else:
        logging.warning("Could not parse variant type from feature table arg, setting to NA")
        vtype = "NA"
    return vtype


def summary_from_featuretable(f, args):

    # define vtype to match those in hap.py
    vtype = resolve_vtype(args)

    # define output data frame and temporary dict
    happy_cols = ["Type", "Filter", "TRUTH.TOTAL", "TRUTH.TP",
        "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK", "FP.gt", "METRIC.Recall",
        "METRIC.Precision", "METRIC.Frac_NA", "METRIC.F1_Score",
        "TRUTH.TOTAL.TiTv_ratio", "QUERY.TOTAL.TiTv_ratio",
        "TRUTH.TOTAL.het_hom_ratio", "QUERY.TOTAL.het_hom_ratio"]
    happy_summary = pandas.DataFrame(columns = happy_cols)
    r = dict()
    r["Type"] = vtype
    decimals = 4

    # truth numbers do not depend on filter
    truth = f[(f["REF.truth"].notnull()) & (f["REF.truth"] != "") & (f["REF.truth"] != ".")]
    r["TRUTH.TOTAL"] = truth.shape[0]
    r["TRUTH.TP"] = truth[(truth["tag"] == "TP")].shape[0]
    r["TRUTH.FN"] = truth[(truth["tag"] == "FN")].shape[0]
    r["TRUTH.TOTAL.het_hom_ratio"] = "NA"
    r["TRUTH.TOTAL.TiTv_ratio"] = "NA"

    # query numbers can be from PASS / ALL variants
    for vfilter in ["PASS", "ALL"]:
        r["Filter"] = vfilter

        query = f[(f["REF"].notnull()) & (f["REF"] != "") & (f["REF"] != ".") & \
                  (f["ALT"].notnull()) & (f["ALT"] != "") & (f["ALT"] != ".")]
        if vfilter == "PASS":
            query = query[(query["FILTER"].isnull()) | (query["FILTER"] == "PASS") | (query["FILTER"] == "") | (query["FILTER"] == ".")]

        # query counts across all AFs
        r["QUERY.FP"] = query[(query["tag"] == "FP")].shape[0]
        r["QUERY.UNK"] = query[(query["tag"] == "UNK")].shape[0] + query[(query["tag"] == "AMBI")].shape[0]
        r["QUERY.TOTAL"] = r["TRUTH.TP"] + r["QUERY.FP"] + r["QUERY.UNK"]
        r["FP.gt"] = "NA"
        r["QUERY.TOTAL.het_hom_ratio"] = "NA"
        r["QUERY.TOTAL.TiTv_ratio"] = "NA"

        # performance metrics across all AFs
        r["METRIC.Recall"] = np.round(np.true_divide(r["TRUTH.TP"], r["TRUTH.TOTAL"]), decimals)
        r["METRIC.Precision"] = np.round(np.true_divide(r["TRUTH.TP"], (r["TRUTH.TP"] + r["QUERY.FP"])), decimals)
        r["METRIC.Frac_NA"] = np.round(np.true_divide(r["QUERY.UNK"], r["QUERY.TOTAL"]), decimals)
        r["METRIC.F1_Score"] = np.round(np.true_divide(2 * r["METRIC.Precision"] * r["METRIC.Recall"], (r["METRIC.Precision"] + r["METRIC.Recall"])), decimals)

        happy_summary = happy_summary.append(pandas.DataFrame([r]))

    # reorder columns and return
    happy_summary = happy_summary[happy_cols]

    # append breaks col types - fix them
    int_cols = ["TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK"]
    for col in int_cols:
        happy_summary[col] = happy_summary[col].astype("int")

    # return
    return happy_summary


def extended_from_featuretable(f, args):
    # define vtype to match those in hap.py
    vtype = resolve_vtype(args)

    # define output data frame and temporary dict
    happy_cols = ["Type", "Subtype", "Subset", "Filter",
        "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
        "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK", "FP.gt", "METRIC.Recall",
        "METRIC.Precision", "METRIC.Frac_NA", "METRIC.F1_Score",
        "TRUTH.TOTAL.TiTv_ratio", "QUERY.TOTAL.TiTv_ratio",
        "TRUTH.TOTAL.het_hom_ratio", "QUERY.TOTAL.het_hom_ratio"]
    happy_extended = pandas.DataFrame(columns = happy_cols)
    r = dict()
    decimals = 4

    # set up AF iterators
    start = 0.0
    end = 1.0
    current_binsize = args.af_strat_binsize[0]
    next_binsize = 0
    af_t_feature = args.af_strat_truth
    af_q_feature = args.af_strat_query

    # iterate over AF bins
    while start < 1.0:
        # static columns
        r["Subtype"] = "*"
        r["Type"] = vtype

        # include 1 in last interval
        end = start + current_binsize
        r["Subset"] = "[%.2f,%.2f)" % (start, end)
        if end >= 1:
            end = 1.00000001
            r["Subset"] = "[%.2f,1.00]" % start

        # truth numbers do not depend on filter
        truth = f[(f["REF.truth"].notnull()) & (f["REF.truth"] != "") & (f["REF.truth"] != ".") & \
                  (f[af_t_feature] >= start) & (f[af_t_feature] < end)]
        r["TRUTH.TOTAL"] = truth.shape[0]
        r["TRUTH.TP"] = truth[(truth["tag"] == "TP")].shape[0]
        r["TRUTH.FN"] = truth[(truth["tag"] == "FN")].shape[0]
        r["TRUTH.TOTAL.het_hom_ratio"] = "NA"
        r["TRUTH.TOTAL.TiTv_ratio"] = "NA"

        # query numbers can be from PASS / ALL variants
        for vfilter in ["PASS", "ALL"]:
            r["Filter"] = vfilter

            query = f[(f["REF"].notnull()) & (f["REF"] != "") & (f["REF"] != ".") & \
                      (f["ALT"].notnull()) & (f["ALT"] != "") & (f["ALT"] != ".") & \
                      (f[af_q_feature] >= start) & (f[af_q_feature] < end)]
            if vfilter == "PASS":
                query = query[(query["FILTER"].isnull()) | (query["FILTER"] == "PASS") | (query["FILTER"] == "") | (query["FILTER"] == ".")]

            # query counts across all AFs
            r["QUERY.FP"] = query[(query["tag"] == "FP")].shape[0]
            r["QUERY.UNK"] = query[(query["tag"] == "UNK")].shape[0] + query[(query["tag"] == "AMBI")].shape[0]
            r["QUERY.TOTAL"] = r["TRUTH.TP"] + r["QUERY.FP"] + r["QUERY.UNK"]
            r["FP.gt"] = "NA"
            r["QUERY.TOTAL.het_hom_ratio"] = "NA"
            r["QUERY.TOTAL.TiTv_ratio"] = "NA"

            # performance metrics across all AFs
            r["METRIC.Recall"] = np.round(np.true_divide(r["TRUTH.TP"], r["TRUTH.TOTAL"]), decimals)
            r["METRIC.Precision"] = np.round(np.true_divide(r["TRUTH.TP"], (r["TRUTH.TP"] + r["QUERY.FP"])), decimals)
            r["METRIC.Frac_NA"] = np.round(np.true_divide(r["QUERY.UNK"], r["QUERY.TOTAL"]), decimals)
            r["METRIC.F1_Score"] = np.round(np.true_divide(2 * r["METRIC.Precision"] * r["METRIC.Recall"], (r["METRIC.Precision"] + r["METRIC.Recall"])), decimals)
            happy_extended = happy_extended.append(pandas.DataFrame([r]))

        # re-configure counters
        start = end
        next_binsize += 1
        if next_binsize >= len(args.af_strat_binsize):
            next_binsize = 0
        current_binsize = args.af_strat_binsize[next_binsize]

    # reorder columns and return
    happy_extended = happy_extended[happy_cols]

    # append breaks col types - fix them
    int_cols = ["TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK"]
    for col in int_cols:
        happy_extended[col] = happy_extended[col].astype("int")

    # return
    return happy_extended


# noinspection PyBroadException
def main():

    args = parse_args()

    if args.scratch_prefix:
        scratch = os.path.abspath(args.scratch_prefix)
        args.delete_scratch = False
        Tools.mkdir_p(scratch)
    else:
        scratch = tempfile.mkdtemp()

    logging.info("Scratch path is %s" % scratch)
    try:
        bams = []
        md = None
        for x in args.bams:
            bams.append(bamStats(x))

        if bams:
            bres = pandas.concat(bams).groupby("CHROM").mean()

            md = {}

            for x in bres.index:
                logging.info("Mean coverage on %s is %f" % (x, bres.loc[x]["COVERAGE"]))
                md[x] = float(bres.loc[x]["COVERAGE"]) * 3.0

        logging.info("Normalizing/reading inputs")

        ntpath = os.path.join(scratch, "normalized_truth.vcf.gz")

        if not (args.cont and os.path.exists(ntpath)):
            preprocessVCF(args.truth, ntpath, args.location,
                          True,  # pass_only
                          args.fixchr_truth,  # chrprefix
                          args.normalize_truth,  # norm,
                          args.regions_bedfile,
                          args.targets_bedfile,
                          args.ref)
        else:
            logging.info("Continuing from %s" % ntpath)

        if not (args.cont and os.path.exists(ntpath + ".csi")):
            runBcftools("index", ntpath)

        nqpath = os.path.join(scratch, "normalized_query.vcf.gz")

        if not (args.cont and os.path.exists(nqpath)):
            preprocessVCF(args.query, nqpath, args.location,
                          not args.inc_nonpass,  # pass_only
                          args.fixchr_query,  # chrprefix
                          args.normalize_query,  # norm,
                          args.regions_bedfile,
                          args.targets_bedfile,
                          args.ref)
        else:
            logging.info("Continuing from %s" % nqpath)

        if not (args.cont and os.path.exists(nqpath + ".csi")):
            runBcftools("index", nqpath)

        logging.info("Intersecting")

        tpfn_files = all([os.path.exists(os.path.join(scratch, "tpfn", "0000.vcf.gz")),
                          os.path.exists(os.path.join(scratch, "tpfn", "0001.vcf.gz")),
                          os.path.exists(os.path.join(scratch, "tpfn", "0002.vcf.gz"))])

        tpfn_r_files = all([os.path.exists(os.path.join(scratch, "tpfn", "0000.vcf.gz")),
                            os.path.exists(os.path.join(scratch, "tpfn", "0001.vcf.gz")),
                            os.path.exists(os.path.join(scratch, "tpfn", "0002.vcf.gz"))])

        if not (args.cont and tpfn_files):
            runBcftools("isec", ntpath, nqpath, "-p", os.path.join(scratch, "tpfn"), "-O", "z")
        else:
            logging.info("Continuing from %s" % os.path.join(scratch, "tpfn"))

        if args.features and not (args.cont and tpfn_r_files):
            # only need to do this for getting the feature table
            runBcftools("isec", nqpath, ntpath, "-p", os.path.join(scratch, "tpfn_r"), "-O", "z")

        logging.info("Getting FPs / Ambi / Unk")

        fppath = os.path.join(scratch, "fp.vcf.gz")
        unkpath = os.path.join(scratch, "unk.vcf.gz")
        ambipath = os.path.join(scratch, "ambi.vcf.gz")

        # get header to print to unk and ambi VCFs
        rununiquepath = os.path.join(scratch, "tpfn", "0001.vcf.gz")
        header = runBcftools("view", rununiquepath, "--header-only")

        fp = Tools.BGZipFile(fppath, True)
        fp.write(header)

        unk = Tools.BGZipFile(unkpath, True)
        unk.write(header)

        ambi = Tools.BGZipFile(ambipath, True)
        ambi.write(header)

        ambiClasses = Counter()
        ambiReasons = Counter()

        fpclasses = BedIntervalTree()
        if args.ambi:
            # can have multiple ambiguous BED files
            for aBED in args.ambi:
                # auto-label from first value after chr start end
                # new ambi files have the label in position 4
                # old ones will look weird here.
                fpclasses.addFromBed(aBED, lambda xe: xe[4], args.fixchr_truth)

        if args.FP:
            fpclasses.addFromBed(args.FP, "FP", args.fixchr_truth)

        # split VCF into FP, UNK and AMBI
        toProcess = gzip.open(rununiquepath, "rb")
        for entry in toProcess:
            if entry[0] == '#':
                continue

            fields = entry.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            stop = int(fields[1]) + len(fields[3])

            overlap = fpclasses.intersect(chrom, start, stop)

            is_fp = False
            is_ambi = False

            classes_this_pos = set()

            for o in overlap:
                reason = o.value[0]
                if reason == "fp" and args.ambi_fp:
                    reason = "FP"
                elif reason == "fp":
                    reason = "ambi-fp"
                elif reason == "unk":
                    reason = "ambi-unk"

                classes_this_pos.add(reason)
                try:
                    ambiReasons["%s: rep. count %s" % (reason, o.value[1])] += 1
                except IndexError:
                    ambiReasons["%s: rep. count *" % reason] += 1
                for x in o.value[3:]:
                    ambiReasons["%s: %s" % (reason, x)] += 1
                if reason == "FP":
                    is_fp = True
                else:
                    is_ambi = True

            for reason in classes_this_pos:
                ambiClasses[reason] += 1

            if is_fp:
                fp.write(entry)
            elif is_ambi:
                ambi.write(entry)
            elif not args.count_unk:
                # when we don't have FP regions, unk stuff becomes FP
                fp.write(entry)
            else:
                unk.write(entry)

        toProcess.close()

        # since 0001.vcf.gz should already be sorted, we can just convert to bgzipped vcf
        # and create index
        fp.close()
        ambi.close()
        unk.close()

        runBcftools("index", "--tbi", fppath)
        runBcftools("index", "--tbi", unkpath)
        runBcftools("index", "--tbi", ambipath)

        logging.info("Counting variants...")

        truthcounts = parseStats(runBcftools("stats", ntpath), "total.truth")
        querycounts = parseStats(runBcftools("stats", nqpath), "total.query")

        tpcounts = parseStats(runBcftools("stats", os.path.join(scratch, "tpfn", "0002.vcf.gz")), "tp")
        fncounts = parseStats(runBcftools("stats", os.path.join(scratch, "tpfn", "0000.vcf.gz")), "fn")
        fpcounts = parseStats(runBcftools("stats", fppath), "fp")
        ambicounts = parseStats(runBcftools("stats", ambipath), "ambi")
        unkcounts = parseStats(runBcftools("stats", unkpath), "unk")

        res = pandas.merge(truthcounts, querycounts, on="type")
        res = pandas.merge(res, tpcounts, on="type")
        res = pandas.merge(res, fpcounts, on="type")
        res = pandas.merge(res, fncounts, on="type")
        res = pandas.merge(res, unkcounts, on="type")
        res = pandas.merge(res, ambicounts, on="type")

        # no explicit guarantee that total.query is equal to unk + ambi + fp + tp
        # testSum = res["fp"] + res["tp"] + res["unk"] + res["ambi"]

        # filter and relabel
        res = res[res["type"] != "samples"]
        res = res[res["type"] != "multiallelic SNP sites"]
        res = res[res["type"] != "multiallelic sites"]
        res.loc[res["type"] == "SNPs", "type"] = "SNVs"

        metrics_output = makeMetricsObject("som.py.comparison")

        if args.ambi and args.explain_ambiguous:
            ac = list(ambiClasses.iteritems())
            if ac:
                ambie = pandas.DataFrame(ac, columns=["class", "count"])
                ambie.sort_values(["class"], inplace=True)
                pandas.set_option("display.max_rows", 1000)
                pandas.set_option("display.max_columns", 1000)
                pandas.set_option("display.width", 1000)
                pandas.set_option("display.height", 1100)
                logging.info("FP/ambiguity classes with info (multiple classes can "
                             "overlap):\n" + ambie.to_string(index=False))
                # in default mode, print result summary to stdout
                if not args.quiet and not args.verbose:
                    print "FP/ambiguity classes with info (multiple classes can " \
                          "overlap):\n" + ambie.to_string(index=False)
                ambie.to_csv(args.output + ".ambiclasses.csv")
                metrics_output["metrics"].append(dataframeToMetricsTable("ambiclasses", ambie))
            else:
                logging.info("No ambiguous variants.")

            ar = list(ambiReasons.iteritems())
            if ar:
                ambie = pandas.DataFrame(ar, columns=["reason", "count"])
                ambie.sort_values(["reason"], inplace=True)
                pandas.set_option("display.max_rows", 1000)
                pandas.set_option("display.max_columns", 1000)
                pandas.set_option("display.width", 1000)
                pandas.set_option("display.height", 1100)
                logging.info("Reasons for defining as ambiguous (multiple reasons can overlap):\n" + ambie.to_string(
                    formatters={'reason': '{{:<{}s}}'.format(ambie['reason'].str.len().max()).format}, index=False))
                # in default mode, print result summary to stdout
                if not args.quiet and not args.verbose:
                    print "Reasons for defining as ambiguous (multiple reasons can overlap):\n" + ambie.to_string(
                        formatters={'reason': '{{:<{}s}}'.format(ambie['reason'].str.len().max()).format}, index=False)
                ambie.to_csv(args.output + ".ambireasons.csv")
                metrics_output["metrics"].append(dataframeToMetricsTable("ambireasons", ambie))
            else:
                logging.info("No ambiguous variants.")

        if args.features:
            logging.info("Extracting features...")
            fset = Somatic.FeatureSet.make(args.features)
            fset.setChrDepths(md)

            logging.info("Collecting TP info (1)...")
            tps = fset.collect(os.path.join(scratch, "tpfn", "0002.vcf.gz"), "TP")

            # TP_r is a hint for fset, they are both TPs
            logging.info("Collecting TP info (2)...")
            tps2 = fset.collect(os.path.join(scratch, "tpfn_r", "0002.vcf.gz"), "TP_r")

            # this is slow because it tries to sort
            # ... which we don't need to do since tps1 and tps2 have the same ordering

            logging.info("Sorting...")
            tps.sort_values(["CHROM", "POS"], inplace=True)
            tps2.sort_values(["CHROM", "POS"], inplace=True)
            tps = tps.reset_index(drop=True)
            tps2 = tps2.reset_index(drop=True)

            logging.info("Merging TP info...")
            columns_tps = list(tps)
            columns_tps2 = list(tps2)

            len1 = tps.shape[0]
            len2 = tps2.shape[0]

            if len1 != len2:
                raise Exception("Cannot read TP features, lists have different lengths : %i != %i" % (len1, len2))

            if not args.disable_order_check:
                logging.info("Checking order %i / %i" % (len1, len2))

                for x in xrange(0, len1):
                    for a in ["CHROM", "POS"]:
                        if tps.loc[x][a] != tps2.loc[x][a]:
                            raise Exception("Cannot merge TP features, inputs are out of order at %s / %s" % (
                                str(tps[x:x + 1]), str(tps2[x:x + 1])))

            logging.info("Merging...")

            cdata = {
                "CHROM": tps["CHROM"],
                "POS": tps["POS"],
                "tag": tps["tag"]
            }

            tpc = pandas.DataFrame(cdata, columns=["CHROM", "POS", "tag"])

            all_columns = list(set(columns_tps + columns_tps2))
            for a in all_columns:
                if a in columns_tps and a not in columns_tps2:
                    tpc[a] = tps[a]
                elif a not in columns_tps and a in columns_tps2:
                    tpc[a] = tps2[a]
                elif a not in ["CHROM", "POS", "tag"]:
                    tpc[a] = tps2[a]
                    tpc[a + ".truth"] = tps[a]

            logging.info("Collecting FP info...")
            fps = fset.collect(fppath, "FP")
            ambs = fset.collect(ambipath, "AMBI")

            logging.info("Collecting FN info...")
            fns = fset.collect(os.path.join(scratch, "tpfn", "0000.vcf.gz"), "FN")

            renamed = {}
            tp_cols = list(tpc)
            for col in list(fns):
                if col + ".truth" in tp_cols:
                    renamed[col] = col + ".truth"
            fns.rename(columns=renamed, inplace=True)

            featurelist = [tpc, fps, fns, ambs]

            if unkpath is not None:
                logging.info("Collecting UNK info...")
                unk = fset.collect(unkpath, "UNK")
                featurelist.append(unk)

            logging.info("Making feature table...")
            featuretable = pandas.concat(featurelist)

            # reorder to make more legible
            first_columns = ["CHROM", "POS", "tag"]
            # noinspection PyTypeChecker
            all_columns = list(featuretable)

            if "REF" in all_columns:
                first_columns.append("REF")

            if "REF.truth" in all_columns:
                first_columns.append("REF.truth")

            if "ALT" in all_columns:
                first_columns.append("ALT")

            if "ALT.truth" in all_columns:
                first_columns.append("ALT.truth")

            ordered_columns = first_columns + sorted([x for x in all_columns if x not in first_columns])
            featuretable = featuretable[ordered_columns]
            # make sure positions are integers
            featuretable["POS"] = featuretable["POS"].astype(int)

            logging.info("Saving feature table...")
            featuretable.to_csv(args.output + ".features.csv", float_format='%.8f')

            if args.roc is not None:
                roc_table = args.roc.from_table(featuretable)
                roc_table.to_csv(args.output + ".roc.csv", float_format='%.8f')

            featuretable["FILTER"].fillna("", inplace=True)
            featuretable.ix[featuretable["REF"].str.len() < 1, "absent"] = True
            featuretable.ix[featuretable["tag"] == "FN", "REF"] = featuretable.ix[featuretable["tag"] == "FN",
                                                                                  "REF.truth"]
            featuretable.ix[featuretable["tag"] == "FN", "ALT"] = featuretable.ix[featuretable["tag"] == "FN",
                                                                                  "ALT.truth"]
            af_t_feature = args.af_strat_truth
            af_q_feature = args.af_strat_query
            for vtype in ["records", "SNVs", "indels"]:
                featuretable["vtype"] = resolve_vtype(args)
                featuretable_this_type = featuretable

                if args.count_filtered_fn:
                    res.ix[res["type"] == vtype, "fp.filtered"] = featuretable_this_type[
                        (featuretable_this_type["tag"] == "FP") & (featuretable_this_type["FILTER"] != "")].shape[0]
                    res.ix[res["type"] == vtype, "tp.filtered"] = featuretable_this_type[
                        (featuretable_this_type["tag"] == "TP") & (featuretable_this_type["FILTER"] != "")].shape[0]
                    res.ix[res["type"] == vtype, "unk.filtered"] = featuretable_this_type[
                        (featuretable_this_type["tag"] == "UNK") & (featuretable_this_type["FILTER"] != "")].shape[0]
                    res.ix[res["type"] == vtype, "ambi.filtered"] = featuretable_this_type[
                        (featuretable_this_type["tag"] == "AMBI") & (featuretable_this_type["FILTER"] != "")].shape[0]

                if args.af_strat:
                    start = 0.0
                    end = 1.0
                    current_binsize = args.af_strat_binsize[0]
                    next_binsize = 0
                    while start < 1.0:
                        # include 1 in last interval
                        end = start + current_binsize
                        if end >= 1:
                            end = 1.00000001
                        if start >= end:
                             break
                        n_tp = featuretable_this_type[(featuretable_this_type["tag"] == "TP") &
                                                      (featuretable_this_type[af_t_feature] >= start) &
                                                      (featuretable_this_type[af_t_feature] < end)]
                        n_fn = featuretable_this_type[(featuretable_this_type["tag"] == "FN") &
                                                      (featuretable_this_type[af_t_feature] >= start) &
                                                      (featuretable_this_type[af_t_feature] < end)]
                        n_fp = featuretable_this_type[(featuretable_this_type["tag"] == "FP") &
                                                      (featuretable_this_type[af_q_feature] >= start) &
                                                      (featuretable_this_type[af_q_feature] < end)]
                        n_ambi = featuretable_this_type[(featuretable_this_type["tag"] == "AMBI") &
                                                        (featuretable_this_type[af_q_feature] >= start) &
                                                        (featuretable_this_type[af_q_feature] < end)]
                        n_unk = featuretable_this_type[(featuretable_this_type["tag"] == "UNK") &
                                                       (featuretable_this_type[af_q_feature] >= start) &
                                                       (featuretable_this_type[af_q_feature] < end)]

                        r = {"type": "%s.%f-%f" % (vtype, start, end),
                             "total.truth": n_tp.shape[0] + n_fn.shape[0],
                             "total.query": n_tp.shape[0] + n_fp.shape[0] + n_ambi.shape[0] + n_unk.shape[0],
                             "tp": n_tp.shape[0],
                             "fp": n_fp.shape[0],
                             "fn": n_fn.shape[0],
                             "unk": n_unk.shape[0],
                             "ambi": n_ambi.shape[0]}

                        if args.count_filtered_fn:
                            r["fp.filtered"] = n_fp[n_fp["FILTER"] != ""].shape[0]
                            r["tp.filtered"] = n_tp[n_tp["FILTER"] != ""].shape[0]
                            r["unk.filtered"] = n_unk[n_unk["FILTER"] != ""].shape[0]
                            r["ambi.filtered"] = n_ambi[n_ambi["FILTER"] != ""].shape[0]

                        res = pandas.concat([res, pandas.DataFrame([r])])

                        if args.roc is not None and (n_tp.shape[0] + n_fn.shape[0] + n_fp.shape[0]) > 0:
                            roc_table_strat = args.roc.from_table(pandas.concat([n_tp, n_fp, n_fn]))
                            rtname = "%s.%s.%f-%f.roc.csv" % (args.output, vtype, start, end)
                            roc_table_strat.to_csv(rtname, float_format='%.8f')
                        start = end
                        next_binsize += 1
                        if next_binsize >= len(args.af_strat_binsize):
                            next_binsize = 0
                        current_binsize = args.af_strat_binsize[next_binsize]

        if not args.af_strat:
            res = res[(res["total.truth"] > 0)]

        # summary metrics with confidence intervals
        ci_alpha = 1.0 - args.ci_level

        recall = binomialCI(res["tp"], res["tp"]+res["fn"], ci_alpha)
        precision = binomialCI(res["tp"], res["tp"]+res["fp"], ci_alpha)
        res["recall"], res["recall_lower"], res["recall_upper"] = recall
        res["recall2"] = res["tp"] / (res["total.truth"])
        res["precision"], res["precision_lower"], res["precision_upper"] = precision
        res["na"] = res["unk"] / (res["total.query"])
        res["ambiguous"] = res["ambi"] / res["total.query"]

        any_fp = fpclasses.countbases(label="FP")

        fp_region_count = 0
        auto_size = True
        if args.fpr_size:
            try:
                fp_region_count = int(args.fpr_size)
                auto_size = False
            except:
                pass
        if auto_size:
            if any_fp:
                if args.location:
                    chrom, _, rest = args.location.partition(":")
                    if rest:
                        start, _, end = rest.partition("_")
                        if start:
                            start = int(start)
                        if end:
                            end = int(end)
                    else:
                        fp_region_count += fpclasses.countbases(chrom, label="FP")
                else:
                    fp_region_count = any_fp
            else:
                cs = fastaContigLengths(args.ref)
                if args.location:
                    fp_region_count = calculateLength(cs, args.location)
                else:
                    # use all locations we saw calls on
                    h1 = Tools.vcfextract.extractHeadersJSON(ntpath)
                    h1_chrs = h1["tabix"]["chromosomes"]

                    if not h1_chrs:
                        logging.warn("No contigs in truth file")
                        h1_chrs = []

                    if len(h1_chrs) > 0:
                        qlocations = " ".join(h1_chrs)
                        fp_region_count = calculateLength(cs, qlocations)
                    else:
                        fp_region_count = 0

        res["fp.region.size"] = fp_region_count
        res["fp.rate"] = 1e6 * res["fp"] / res["fp.region.size"]

        if args.count_filtered_fn:
            res["recall.filtered"] = (res["tp"] - res["tp.filtered"]) / (res["tp"] + res["fn"])

            res["precision.filtered"] = (res["tp"] - res["tp.filtered"]) / (res["tp"] - res["tp.filtered"] +
                                                                            res["fp"] - res["fp.filtered"])

            res["fp.rate.filtered"] = 1e6 * (res["fp"] - res["fp.filtered"]) / res["fp.region.size"]

            res["na.filtered"] = (res["unk"] - res["unk.filtered"]) / (res["total.query"])
            res["ambiguous.filtered"] = (res["ambi"] - res["ambi.filtered"]) / res["total.query"]

        # HAP-162 remove inf values
        res.replace([np.inf, -np.inf], 0)

        metrics_output["metrics"].append(dataframeToMetricsTable("result", res))
        vstring = "som.py-%s" % Tools.version

        logging.info("\n" + res.to_string())
        # in default mode, print result summary to stdout
        if not args.quiet and not args.verbose:
            print "\n" + res.to_string()

        res["sompyversion"] = vstring

        vstring = " ".join(sys.argv)
        res["sompycmd"] = vstring

        # save results
        res.to_csv(args.output + ".stats.csv")

        with open(args.output + ".metrics.json", "w") as fp:
            json.dump(metrics_output, fp)

        if args.happy_stats:
            # parse saved feature table as the one in memory has been updated
            featuretable = pandas.read_csv(args.output + ".features.csv", low_memory = False, dtype = {"FILTER": str})

            # hap.py summary.csv
            summary = summary_from_featuretable(featuretable, args)
            summary.to_csv(args.output + ".summary.csv")

            #  hap.py extended.csv
            if args.af_strat:
                extended = extended_from_featuretable(featuretable, args)
                extended.to_csv(args.output + ".extended.csv", index=False, na_rep="NA")

    finally:
        if args.delete_scratch:
            shutil.rmtree(scratch)
        else:
            logging.info("Scratch kept at %s" % scratch)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
