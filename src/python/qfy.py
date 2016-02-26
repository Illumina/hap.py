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
import multiprocessing
import pandas
import json

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools.metric import makeMetricsObject, dataframeToMetricsTable
from Tools.fastasize import fastaContigLengths
import Haplo.quantify
import Haplo.happyroc


def quantify(args):
    """ Run quantify and write tables """
    vcf_name = args.in_vcf[0]

    if not vcf_name or not os.path.exists(vcf_name):
        raise Exception("Cannot read input VCF.")

    json_name = args.reports_prefix + ".counts.json"

    logging.info("Counting variants...")

    output_vcf = args.reports_prefix + ".vcf.gz"

    roc_table = None

    if args.roc:
        roc_table = args.reports_prefix + ".rocdata.tsv"

    counts = Haplo.quantify.run_quantify(vcf_name,
                                         json_name,
                                         output_vcf if args.write_vcf else False,
                                         {"CONF": args.fp_bedfile} if args.fp_bedfile else None,
                                         args.ref,
                                         threads=args.threads,
                                         output_vtc=args.output_vtc,
                                         qtype=args.type,
                                         roc_val=args.roc,
                                         roc_file=roc_table,
                                         clean_info=not args.preserve_info)

    df = pandas.DataFrame(counts)

    metrics_output = makeMetricsObject("%s.comparison" % args.runner)

    if args.write_counts:
        df.to_csv(args.reports_prefix + ".counts.csv")
        metrics_output["metrics"].append(dataframeToMetricsTable("raw.counts", df))

    # calculate precision / recall
    count_types = []

    simplified_numbers = Haplo.quantify.simplify_counts(counts)

    count_types += simplified_numbers.keys()
    count_types = sorted(list(set(count_types)))

    for vtype in count_types:
        if vtype not in simplified_numbers:
            simplified_numbers[vtype] = {}

        simplified_numbers[vtype]["METRIC.Recall"] = 0
        simplified_numbers[vtype]["METRIC.Recall2"] = 0
        simplified_numbers[vtype]["METRIC.Precision"] = 0
        simplified_numbers[vtype]["METRIC.Frac_NA"] = 0

        try:
            simplified_numbers[vtype]["METRIC.Recall"] = \
                float(simplified_numbers[vtype]["TRUTH.TP"]) / \
                float(simplified_numbers[vtype]["TRUTH.TP"] + simplified_numbers[vtype]["TRUTH.FN"])
        except:
            pass

        try:
            simplified_numbers[vtype]["METRIC.Recall2"] = \
                float(simplified_numbers[vtype]["TRUTH.TP"]) / \
                float(simplified_numbers[vtype]["TRUTH.TOTAL"])
        except:
            pass

        try:
            simplified_numbers[vtype]["METRIC.Precision"] = \
                float(simplified_numbers[vtype]["QUERY.TP"]) / \
                float(simplified_numbers[vtype]["QUERY.TP"] + simplified_numbers[vtype]["QUERY.FP"])
        except:
            pass

        try:
            simplified_numbers[vtype]["METRIC.Frac_NA"] = \
                float(simplified_numbers[vtype]["QUERY.UNK"]) / \
                float(simplified_numbers[vtype]["QUERY.TOTAL"])
        except:
            pass

    pandas.set_option("display.width", 120)
    pandas.set_option("display.max_columns", 1000)
    df = pandas.DataFrame(simplified_numbers).transpose()

    vstring = "%s-%s" % (args.runner, Tools.version)
    vstring += " ".join(sys.argv)

    df.loc[vstring] = 0

    summary_columns = ["TRUTH.TOTAL",
                       "QUERY.TOTAL",
                       "METRIC.Recall",
                       "METRIC.Precision",
                       "METRIC.Frac_NA"]

    for additional_column in ["TRUTH.TOTAL.TiTv_ratio",
                              "QUERY.TOTAL.TiTv_ratio",
                              "TRUTH.TOTAL.het_hom_ratio",
                              "QUERY.TOTAL.het_hom_ratio"]:
        if additional_column in df.columns:
            summary_columns.append(additional_column)

    df[summary_columns].to_csv(args.reports_prefix + ".summary.csv")

    metrics_output["metrics"].append(dataframeToMetricsTable("summary.metrics",
                                                             df[summary_columns]))

    if args.write_counts:
        df.to_csv(args.reports_prefix + ".qfy.extended.csv")
        metrics_output["metrics"].append(dataframeToMetricsTable("all.metrics", df))

    essential_numbers = df[summary_columns]

    pandas.set_option('display.max_columns', 500)
    pandas.set_option('display.width', 1000)

    essential_numbers = essential_numbers[essential_numbers.index.isin(
        ["Locations.SNP", "Locations.INDEL"])]

    logging.info("\n" + str(essential_numbers))

    # in default mode, print result summary to stdout
    if not args.quiet and not args.verbose:
        print "Benchmarking Summary:"
        print str(essential_numbers)

    if args.roc:
        res = Haplo.happyroc.roc(roc_table,
                                 args.roc,
                                 args.roc_filter,
                                 args.reports_prefix + ".roc",
                                 args.roc_reversed)

        for t in res.iterkeys():
            rocdf = pandas.read_table(res[t])
            metrics_output["metrics"].append(dataframeToMetricsTable("roc." + t, rocdf))

    with open(args.reports_prefix + ".metrics.json", "w") as fp:
        json.dump(metrics_output, fp)


def updateArgs(parser):
    """ add common quantification args """
    parser.add_argument("-t", "--type", dest="type", choices=["xcmp", "ga4gh"],
                        help="Annotation format in input VCF file.")

    parser.add_argument("-f", "--false-positives", dest="fp_bedfile",
                        default=None, type=str,
                        help="False positive / confident call regions (.bed or .bed.gz).")

    parser.add_argument("-V", "--write-vcf", dest="write_vcf",
                        default=False, action="store_true",
                        help="Write an annotated VCF.")

    parser.add_argument("-X", "--write-counts", dest="write_counts",
                        default=True, action="store_true",
                        help="Write advanced counts and metrics.")

    parser.add_argument("--no-write-counts", dest="write_counts",
                        default=True, action="store_false",
                        help="Do not write advanced counts and metrics.")

    parser.add_argument("--output-vtc", dest="output_vtc",
                        default=False, action="store_true",
                        help="Write VTC field in the final VCF which gives the counts each position has contributed to.")

    parser.add_argument("--preserve-info", dest="preserve_info", action="store_true", default=False,
                        help="When using XCMP, preserve and merge the INFO fields in truth and query. Useful for ROC computation.")

    parser.add_argument("--roc", dest="roc", default=False,
                        help="Select an INFO feature to produce a ROC on. This works best with "
                             "--no-internal-preprocessing and --no-internal-leftshift since these "
                             "flags preserve the most INFO flags from the input files.")

    parser.add_argument("--roc-filter", dest="roc_filter", default=False,
                        help="Select a filter to ignore when making ROCs.")

    parser.add_argument("--roc-reversed", dest="roc_reversed", default=False,
                        help="Change the meaning of the ROC feature to count the other way around (higher values=bad).")


def main():
    parser = argparse.ArgumentParser("Quantify annotated VCFs")

    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Show version number and exit.")

    parser.add_argument("in_vcf", help="VCF file to quantify", nargs=1)

    updateArgs(parser)

    # generic, keep in sync with hap.py!
    parser.add_argument("-o", "--report-prefix", dest="reports_prefix",
                        default=None,
                        help="Filename prefix for report output.")

    parser.add_argument("-r", "--reference", dest="ref", default=None, help="Specify a reference file.")

    parser.add_argument("--threads", dest="threads",
                        default=multiprocessing.cpu_count(), type=int,
                        help="Number of threads to use.")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    args, unknown_args = parser.parse_known_args()

    args.runner = "qfy.py"

    if not args.ref:
        args.ref = Tools.defaultReference()

    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        level=loglevel)

    # remove some safe unknown args
    unknown_args = [x for x in unknown_args if x not in ["--force-interactive"]]
    if len(sys.argv) < 2 or len(unknown_args) > 0:
        if unknown_args:
            logging.error("Unknown arguments specified : %s " % str(unknown_args))
        parser.print_help()
        exit(0)

    if args.version:
        print "qfy.py %s" % Tools.version
        exit(0)

    quantify(args)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
