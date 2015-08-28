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
# Diploid ROC Computation
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
import logging
import traceback
import argparse

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
import Haplo.happyroc
import Tools.roc
import pandas

def main():
    parser = argparse.ArgumentParser("Hap.py ROC Maker")

    parser.add_argument("input", help="VCF output by hap.py -V, or a csv feature table from som.py.", default=[], nargs=1)
    parser.add_argument("output", help="Output file prefix", default=[], nargs=1)

    parser.add_argument("--roc", dest="roc", required=True,
                        help="Select an INFO feature to produce a ROC on.")

    parser.add_argument("--roc-filter", dest="roc_filter", default=False,
                        help="Select one or more filters (separated by comma or semicolon) to "
                             "ignore when making ROCs.")

    parser.add_argument("--roc-reversed", dest="roc_reversed", action="store_true", default=False,
                        help="Change the meaning of the ROC feature to count the other way around (higher values=bad).")

    parser.add_argument("--roc-label", dest="roc_label", default="tag",
                        help="Select the tag column for CSV files (needs to have TP/FP values).")


    parser.add_argument("--roc-filter-column", dest="roc_filtercolumn", default="FILTER",
                        help="Select the filter column for CSV files (default: FILTER).")

    args = parser.parse_args()
    print args

    if args.input[0].endswith(".vcf") or args.input[0].endswith(".vcf.gz") or args.input[0].endswith(".bcf"):
        Haplo.happyroc.roc(args.input[0], args.roc, args.roc_filter, args.output[0], args.roc_reversed)
    elif args.input[0].endswith(".csv"):
        tbl = pandas.read_csv(args.input[0])
        roctbl = Tools.roc.tableROC(tbl, 
                                    args.roc_label, 
                                    args.roc, 
                                    args.roc_filtercolumn,
                                    args.roc_filter,
                                    args.roc_reversed)
        roctbl.to_csv(args.output[0] + ".roc.tsv", float_format='%.8f', sep="\t", index=False)
    else:
        raise Exception("Unknown extension for the input file (use .csv, .vcf, .vcf.gz or .bcf).")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
