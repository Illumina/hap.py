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

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
import Haplo.happyroc

def main():
    parser = argparse.ArgumentParser("Hap.py ROC Maker")

    parser.add_argument("input", help="VCF output by hap.py -V", default=[], nargs=1)
    parser.add_argument("output", help="VCF output by hap.py -V", default=[], nargs=1)

    parser.add_argument("--roc", dest="roc", required=True,
                        help="Select an INFO feature to produce a ROC on. This works best with "
                             "--no-internal-preprocessing and --no-internal-leftshift since these "
                             "flags preserve the most INFO flags from the input files.")

    parser.add_argument("--roc-filter", dest="roc_filter", default=False,
                        help="Select a filter to ignore when making ROCs.")

    args = parser.parse_args()
    print args

    Haplo.happyroc.roc(args.input[0], args.roc, args.roc_filter, args.output[0])

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
