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
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt
#
# 9/9/2014
#
# Diploid VCF File Comparison
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
import subprocess
import multiprocessing
import gzip
import tempfile
import time
import pandas
import json

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools import vcfextract
from Tools.bcftools import preprocessVCF, bedOverlapCheck
from Tools.parallel import runParallel
from Tools.metric import makeMetricsObject, dataframeToMetricsTable

import Haplo.quantify


def blocksplitWrapper(location_str, args):
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     dir=args.scratch_prefix,
                                     prefix="result.%s" % location_str,
                                     suffix=".chunks.bed")
    tf.close()

    to_run = "blocksplit %s %s -l %s -o %s --window %i --nblocks %i -f %i" % \
                (args.vcf1,
                 args.vcf2,
                 location_str,
                 tf.name,
                 args.window,
                 args.pieces,
                 0 if args.usefiltered else 1)

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
    logging.info("blocksplit for %s -- time taken %.2f" % (location_str, elapsed))
    return tf.name


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

    to_run = "xcmp %s %s -l %s -o %s %s -r %s -f %i -n %i -V %i --expand-hapblocks %i --window %i" % \
                (args.vcf1,
                 args.vcf2,
                 location_str,
                 tf.name,
                 bname,
                 args.ref,
                 0 if args.usefiltered else 1,
                 args.max_enum,
                 1 if args.int_preprocessing else 0,
                 args.hb_expand,
                 args.window)

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


def main():
    parser = argparse.ArgumentParser("Haplotype Comparison")

    # input
    parser.add_argument('--location', '-l', dest='locations', required=False, default=None,
                        help='Add a location to the compare list (when not given, will use chr1-22, chrX, chrY).')

    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Show version number and exit.")

    parser.add_argument("-P", "--include-nonpass", dest="usefiltered", action="store_true", default=False,
                        help="Use to include failing variants in comparison.")

    parser.add_argument("-R", "--restrict-regions", dest="regions_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (sparse) regions (using -R in bcftools).")

    parser.add_argument("-T", "--target-regions", dest="targets_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (dense) regions (using -T in bcftools).")

    parser.add_argument("-f", "--false-positives", dest="fp_bedfile",
                        default=None, type=str,
                        help="False positive / confident call regions (.bed or .bed.gz).")

    parser.add_argument("-r", "--reference", dest="ref", default=None, help="Specify a reference file.")

    # output
    parser.add_argument("-o", "--report-prefix", dest="reports_prefix",
                        default=None,
                        help="Filename prefix for report output.")

    parser.add_argument("-V", "--write-vcf", dest="write_vcf",
                        default=False, action="store_true",
                        help="Write an annotated VCF.")

    parser.add_argument("-B", "--write-bed", dest="write_bed",
                        default=False, action="store_true",
                        help="Write a bed file with the haplotype blocks that were used.")

    parser.add_argument("-X", "--write-counts", dest="write_counts",
                        default=False, action="store_true",
                        help="Write advanced counts and metrics.")

    parser.add_argument("--scratch-prefix", dest="scratch_prefix",
                        default=None,
                        help="Directory for scratch files.")

    parser.add_argument("--keep-scratch", dest="delete_scratch",
                        default=True, action="store_false",
                        help="Filename prefix for scratch report output.")

    # detailed control of comparison
    parser.add_argument("--preprocess-truth", dest="preprocessing_truth", action="store_true", default=False,
                        help="Preprocess truth file using bcftools.")

    parser.add_argument("--external-preprocessing", dest="preprocessing", action="store_true", default=False,
                        help="Perform VCF preprocessing using bcftools.")

    parser.add_argument("--bcftools-norm", dest="preprocessing_norm", action="store_true", default=False,
                        help="Enable preprocessing through bcftools norm -c x -D (requires external preprocessing to be switched on).")

    parser.add_argument("--fixchr-truth", dest="fixchr_truth", action="store_true", default=None,
                        help="Add chr prefix to truth file (default: auto).")

    parser.add_argument("--fixchr-query", dest="fixchr_query", action="store_true", default=None,
                        help="Add chr prefix to query file (default: auto).")

    parser.add_argument("--no-fixchr-truth", dest="fixchr_truth", action="store_false",
                        help="Disable chr replacement for truth (default: auto).")

    parser.add_argument("--no-fixchr-query", dest="fixchr_query", action="store_false",
                        help="Add chr prefix to query file (default: auto).")

    parser.add_argument("--no-internal-preprocessing", dest="int_preprocessing", action="store_false", default=True,
                        help="Switch off xcmp's internal VCF leftshift preprocessing.")

    parser.add_argument("--no-auto-index", dest="auto_index", action="store_false", default=True,
                        help="Disable automatic index creation for input files. "
                             "The index is only necessary at this stage if we want to auto-detect locations. "
                             "When used with -l, and when it is known that there are variants at all given locations "
                             "this is not needed and can be switched off to save time.")

    parser.add_argument("-w", "--window-size", dest="window",
                        default=30, type=int,
                        help="Window size for haplotype block finder.")

    parser.add_argument("--enumeration-threshold", dest="max_enum",
                        default=2048, type=int,
                        help="Enumeration threshold / maximum number of sequences to enumerate per block.")

    parser.add_argument("-e", "--expand-hapblocks", dest="hb_expand",
                        default=30, type=int,
                        help="Expand haplotype blocks by this many basepairs left and right.")

    parser.add_argument("--threads", dest="threads",
                        default=multiprocessing.cpu_count(), type=int,
                        help="Number of threads to use.")

    if Tools.has_sge:
        parser.add_argument("--force-interactive", dest="force_interactive",
                            default=False, action="store_true",
                            help="Force running interactively (i.e. when JOB_ID is not in the environment)")

    parser.add_argument("_vcfs", help="Two VCF files.", default=[], nargs="*")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    args, _ = parser.parse_known_args()

    if not Tools.has_sge:
        args.force_interactive = True

    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

    ## reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        level=loglevel)

    if len(sys.argv) < 2:
        parser.print_help()
        exit(0)

    if args.version:
        if Tools.has_muscle:
            print "Hap.py %s-muscle" % Tools.version
        else:
            print "Hap.py %s-no-muscle" % Tools.version
        exit(0)

    # sanity-check regions bed file (HAP-57)
    if args.regions_bedfile:
        logging.info("Checking input regions.")
        if bedOverlapCheck(args.regions_bedfile):
            raise Exception("The regions bed file (specified using -R) has overlaps, this will not work with xcmp."
                            " You can either use -T, or run the file through bedtools merge")
        args.preprocessing_truth = True
        args.preprocessing = True

    if args.targets_bedfile:
        args.preprocessing_truth = True
        args.preprocessing = True

    if args.fp_bedfile and not os.path.exists(args.fp_bedfile):
        raise Exception("FP/confident call region bed file does not exist.")

    tempfiles = []

    try:
        if not args.force_interactive and not "JOB_ID" in os.environ:
            parser.print_help()
            raise Exception("Please qsub me so I get approximately 1 GB of RAM per thread.")

        if not args.ref:
            args.ref = Tools.defaultReference()

        if not os.path.exists(args.ref):
            raise Exception("Please specify a valid reference path using -r.")

        if not args.reports_prefix:
            raise Exception("Please specify an output prefix using -o ")

        # noinspection PyProtectedMember
        if not args._vcfs or len(args._vcfs) != 2:
            raise Exception("Please specify exactly two input VCFs.")

        # noinspection PyProtectedMember
        args.vcf1 = args._vcfs[0]
        # noinspection PyProtectedMember
        args.vcf2 = args._vcfs[1]

        if not os.path.exists(args.vcf1):
            raise Exception("Input file %s does not exist." % args.vcf1)
        if not os.path.exists(args.vcf2):
            raise Exception("Input file %s does not exist." % args.vcf2)

        logging.info("Comparing %s and %s" % (args.vcf1, args.vcf2))

        h1 = vcfextract.extractHeadersJSON(args.vcf1)
        if args.auto_index and not h1["tabix"]:
            logging.info("Creating indexed version of %s -- consider creating an index beforehand to save time here." % args.vcf1)
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="truth.ix",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            tempfiles.append(vtf.name + ".tbi")
            args.vcf1 = Tools.bcftools.makeIndex(args.vcf1, vtf.name)
            h1 = vcfextract.extractHeadersJSON(args.vcf1)

        h2 = vcfextract.extractHeadersJSON(args.vcf2)
        if args.auto_index and not h2["tabix"]:
            logging.info("Creating indexed version of %s -- consider creating an index beforehand to save time here." % args.vcf2)
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="query.ix",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            tempfiles.append(vtf.name + ".tbi")
            args.vcf2 = Tools.bcftools.makeIndex(args.vcf2, vtf.name)
            h2 = vcfextract.extractHeadersJSON(args.vcf2)

        if args.locations is None or len(args.locations) == 0:
            # all chromosomes
            args.locations = ["chr" + x for x in map(str, range(1, 23))]

        if type(args.locations) is not list and args.locations is not None:
            # noinspection PyUnresolvedReferences
            args.locations = args.locations.split(",")

        if not h1["tabix"]:
            args.preprocessing_truth = True
            logging.warn("Truth file is not Tabix indexed. Switching on pre-processing + chr name conversion.")
            if args.fixchr_truth is None:
                args.fixchr_truth = True
        elif args.fixchr_truth is None:
            # autodetect chr naming
            count_with_fix = len([__ for __ in h1["tabix"]["chromosomes"]
                                     if ("chr%s" % str(_)) in args.locations])
            count_no_fix = len([__ for __ in h1["tabix"]["chromosomes"] if str(_) in args.locations])
            logging.info("Truth: Number of chromosome names matching with / without renaming : %i / %i " % (
                count_with_fix, count_no_fix))
            if count_with_fix > count_no_fix:
                args.fixchr_truth = True
                logging.info("Will fix chromosome names (truth).")
            else:
                logging.info("Will not fix chromosome names (truth).")
                args.fixchr_truth = False

        if not h2["tabix"]:
            args.preprocessing = True
            logging.warn("Query file is not Tabix indexed. Switching on pre-processing + chr name conversion.")
            # don't overwrite setting, but if it's None, replace with True to be sure
            if args.fixchr_query is None:
                args.fixchr_query = True
        elif args.fixchr_query is None:
            # autodetect chr naming
            count_with_fix = len([__ for __ in h2["tabix"]["chromosomes"]
                                   if ("chr%s" % str(_)) in args.locations])
            count_no_fix = len([__ for __ in h2["tabix"]["chromosomes"] if str(_) in args.locations])
            logging.info("Query: Number of chromosome names matching with / without renaming : %i / %i " % (
                count_with_fix, count_no_fix))
            if count_with_fix > count_no_fix:
                args.fixchr_query = True
                logging.info("Will fix chromosome names (query).")
            else:
                logging.info("Will not fix chromosome names (query).")
                args.fixchr_query = False

        if args.fixchr_truth or args.preprocessing_norm:
            args.preprocessing_truth = True

        if args.fixchr_query or args.preprocessing_norm:
            args.preprocessing = True

        if args.preprocessing_truth:
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="truth.pp",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            preprocessVCF(args.vcf1, vtf.name, ",".join(args.locations),
                          not args.usefiltered,     # pass_only
                          args.fixchr_truth,        # chrprefix
                          args.preprocessing_norm,  # norm,
                          args.regions_bedfile,
                          args.targets_bedfile,
                          args.ref)
            args.vcf1 = vtf.name
            # get headers again if we preprocessed
            h1 = vcfextract.extractHeadersJSON(args.vcf1)

        if args.preprocessing:
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="query.pp",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            preprocessVCF(args.vcf2, vtf.name, ",".join(args.locations),
                          not args.usefiltered,     # pass_only
                          args.fixchr_query,        # chrprefix
                          args.preprocessing_norm,  # norm,
                          args.regions_bedfile,
                          args.targets_bedfile,
                          args.ref)
            args.vcf2 = vtf.name
            # get headers again if we preprocessed
            h2 = vcfextract.extractHeadersJSON(args.vcf2)

        if not h1["tabix"]:
            raise Exception("Truth file is not Tabix indexed.")

        if not h2["tabix"]:
            raise Exception("Truth file is not Tabix indexed.")

        newlocations = []

        if not h1["tabix"]["chromosomes"]:
            h1["tabix"]["chromosomes"] = []
        if not h2["tabix"]["chromosomes"]:
            h2["tabix"]["chromosomes"] = []

        for xc in args.locations:
            if xc not in h1["tabix"]["chromosomes"]:
                logging.warn("No calls for location %s in truth!" % xc)
            if xc not in h2["tabix"]["chromosomes"]:
                logging.warn("No calls for location %s in query!" % xc)

            if (xc not in h1["tabix"]["chromosomes"]) and (xc not in h2["tabix"]["chromosomes"]):
                logging.warn("Removing location %s because neither input file has calls there." % xc)
            else:
                newlocations.append(xc)

        if not newlocations:
            raise Exception("Location list is empty: the input files do not appear to have variants on any of %s" %
                            str(args.locations))

        args.locations = newlocations

        if args.threads > 1:
            logging.info("Running using %i parallel processes." % args.threads)
            pool = multiprocessing.Pool(int(args.threads))

            # find balanced pieces
            args.pieces = (args.threads + len(args.locations) - 1) / len(args.locations)
            res = runParallel(pool, blocksplitWrapper, args.locations, args)

            if None in res:
                raise Exception("One of the blocksplit processes failed.")

            tempfiles += res

            args.locations = []
            for f in res:
                with open(f) as fp:
                    for l in fp:
                        ll = l.strip().split("\t", 3)
                        if len(ll) < 3:
                            continue
                        xchr = ll[0]
                        start = int(ll[1]) + 1
                        end = int(ll[2])
                        args.locations.append("%s:%i-%i" % (xchr, start, end))
        else:
            pool = None


        # count variants before normalisation
        if not "samples" in h1 or not h1["samples"]:
            raise Exception("Cannot read sample names from truth input file")
        counts_truth = Haplo.quantify.run_quantify(args.vcf1,
                                                   None,
                                                   None,
                                                   {"CONF": args.fp_bedfile} if args.fp_bedfile else None,
                                                   args.ref,
                                                   h1["samples"][0])

        if not "samples" in h2 or not h2["samples"]:
            raise Exception("Cannot read sample names from truth input file")
        counts_query = Haplo.quantify.run_quantify(args.vcf2,
                                                   None,
                                                   None,
                                                   {"CONF": args.fp_bedfile} if args.fp_bedfile else None,
                                                   args.ref,
                                                   h2["samples"][0])
        # do xcmp
        res = runParallel(pool, xcmpWrapper, args.locations, args)
        tempfiles += [x[0] for x in res if x is not None]   # VCFs
        tempfiles += [x[1] for x in res if x is not None and x[1] is not None]   # beds (if any)

        if None in res:
            raise Exception("One of the xcmp jobs failed.")

        tf = tempfile.NamedTemporaryFile(delete=False,
                                         dir=args.scratch_prefix,
                                         prefix="hap.py.result.", suffix=".vcf.gz")
        tf.close()
        tempfiles.append(tf.name)
        output_name = tf.name

        if len(res) == 0:
            raise Exception("Input files/regions do not contain variants (0 haplotype blocks were processed).")

        # concatenate + index
        bedfiles = [x[1] for x in res if x is not None and x[1] is not None]
        if args.write_bed and bedfiles:
            runme = " ".join(["cat"] +
                             bedfiles +
                             [">", args.reports_prefix + ".blocks.bed"])
            logging.info("Concatenating block files: %s..." % runme)
            subprocess.check_call(runme,
                                  shell=True)

        logging.info("Concatenating variants...")
        runme_list = [x[0] for x in res if x is not None]
        if len(runme_list) == 0:
            raise Exception("No outputs to concatenate!")

        fo = Tools.BGZipFile(output_name, True)
        for i, x in enumerate(runme_list):
            f = gzip.GzipFile(x)
            for l in f:
                if i == 0 or not l[0] == "#":
                    fo.write(l)
        fo.close()

        logging.info("Indexing...")
        to_run = "tabix -p vcf %s" % output_name
        logging.info("Running '%s'" % to_run)
        subprocess.check_call(to_run, shell=True)

        if args.write_counts:
            json_name = args.reports_prefix + ".counts.json"
        else:
            tf = tempfile.NamedTemporaryFile(delete=False,
                                             dir=args.scratch_prefix,
                                             prefix="counts.",
                                             suffix=".json")
            tf.close()
            json_name = tf.name

        logging.info("Counting variants...")

        counts = Haplo.quantify.run_quantify(output_name,
                                             json_name,
                                             args.reports_prefix + ".vcf.gz" if args.write_vcf else False,
                                             {"CONF": args.fp_bedfile} if args.fp_bedfile else None,
                                             args.ref)

        df = pandas.DataFrame(counts)
        if args.write_counts:
            df.to_csv(args.reports_prefix + ".counts.csv")

        metrics_output = makeMetricsObject("hap.py.comparison")

        if args.write_counts:
            metrics_output["metrics"].append(dataframeToMetricsTable("raw.counts", df))

        # calculate precision / recall
        simplified_truth_counts = Haplo.quantify.simplify_counts(counts_truth, h1["samples"][0:1])
        simplified_query_counts = Haplo.quantify.simplify_counts(counts_query, h2["samples"][0:1])
        simplified_numbers = Haplo.quantify.simplify_counts(counts)

        for vtype in simplified_numbers.keys():
            for suffix in ["", ".HC"]:
                simplified_numbers[vtype]["METRIC.Recall" + suffix] = 0
                simplified_numbers[vtype]["METRIC.Recall2" + suffix] = 0
                simplified_numbers[vtype]["METRIC.Precision" + suffix] = 0
                simplified_numbers[vtype]["METRIC.Frac_NA" + suffix] = 0

                try:
                    simplified_numbers[vtype]["METRIC.Recall" + suffix] = \
                        float(simplified_numbers[vtype]["TRUTH.TP" + suffix]) / \
                        float(simplified_numbers[vtype]["TRUTH.TP" + suffix] + simplified_numbers[vtype]["TRUTH.FN" + suffix])
                except:
                    pass

                try:
                    simplified_numbers[vtype]["METRIC.Recall2" + suffix] = \
                        float(simplified_numbers[vtype]["TRUTH.TP" + suffix]) / float(simplified_numbers[vtype]["TRUTH.TOTAL"])
                except:
                    pass

                try:
                    simplified_numbers[vtype]["METRIC.Precision" + suffix] = \
                        float(simplified_numbers[vtype]["QUERY.TP" + suffix]) / \
                        float(simplified_numbers[vtype]["QUERY.TP" + suffix] + simplified_numbers[vtype]["QUERY.FP" + suffix])
                except:
                    pass

                try:
                    simplified_numbers[vtype]["METRIC.Frac_NA" + suffix] = \
                        float(simplified_numbers[vtype]["QUERY.UNK" + suffix]) / float(simplified_numbers[vtype]["QUERY.TOTAL"])
                except:
                    pass

                try:
                    simplified_numbers[vtype]["TRUTH.TOTAL.RAW"] = simplified_truth_counts[vtype][h1["samples"][0] + ".TOTAL"]
                except:
                    pass

                try:
                    simplified_numbers[vtype]["QUERY.TOTAL.RAW"] = simplified_query_counts[vtype][h2["samples"][0] + ".TOTAL"]
                except:
                    pass


        pandas.set_option("display.width", 120)
        pandas.set_option("display.max_columns", 1000)
        df = pandas.DataFrame(simplified_numbers).transpose()

        for x in df:
            # everything not a metric is a count
            if not x.startswith("METRIC"):
                df[x] = df[x].astype("int64")

        df[["TRUTH.TOTAL",
            "QUERY.TOTAL",
            "METRIC.Recall.HC",
            "METRIC.Recall2.HC",
            "METRIC.Precision.HC",
            "METRIC.Frac_NA.HC",]].to_csv(args.reports_prefix + ".summary.csv")

        metrics_output["metrics"].append(dataframeToMetricsTable("summary.metrics",
                                         df[["TRUTH.TOTAL",
                                             "QUERY.TOTAL",
                                             "METRIC.Recall.HC",
                                             "METRIC.Recall2.HC",
                                             "METRIC.Precision.HC",
                                             "METRIC.Frac_NA.HC",]]))

        if args.write_counts:
            df.to_csv(args.reports_prefix + ".extended.csv")
            metrics_output["metrics"].append(dataframeToMetricsTable("all.metrics", df))

        essential_numbers = df[["TRUTH.TOTAL",
                                "QUERY.TOTAL",
                                "METRIC.Recall.HC",
                                "METRIC.Precision.HC",
                                "METRIC.Frac_NA.HC",]]

        pandas.set_option('display.max_columns', 500)
        pandas.set_option('display.width', 1000)

        essential_numbers = essential_numbers[essential_numbers.index.isin(
            ["Alleles.SNP",
             "Alleles.INS",
             "Alleles.DEL",
             "Locations.SNP.het",
             "Locations.SNP.homalt",])]

        logging.info("\n" + str(essential_numbers))

        # in default mode, print result summary to stdout
        if not args.quiet and not args.verbose:
            print "Benchmarking Summary:"
            print str(essential_numbers)

        with open(args.reports_prefix + ".metrics.json", "w") as fp:
            json.dump(metrics_output, fp)
    finally:
        if args.delete_scratch:
            for x in tempfiles:
                try:
                    os.remove(x)
                except:
                    pass
        else:
            logging.info("Scratch files kept : %s" % (str(tempfiles)))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
