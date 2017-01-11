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
# 23/11/2015
#
# Calculate header-based statistics for a BAM file
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
import json
import pandas

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
import Tools.bamstats


def main():
    parser = argparse.ArgumentParser("BAM Statistics")

    parser.add_argument("input", nargs="+")

    parser.add_argument("-o", "--output", dest="output", help="Output CSV file name")

    args = parser.parse_args()

    all_results = []
    for inp in args.input:
        bs = Tools.bamstats.bamStats(inp)
        bs["FILE"] = inp
        all_results.append(bs)

    all_results = pandas.concat(all_results)
    all_results.reset_index(inplace=True)

    pandas.set_option("display.width", 10000)
    pandas.set_option("display.max_colwidth", 2048)
    pandas.set_option("display.max_columns", 10000)
    pandas.set_option("display.max_rows", 1000000)
    print str(all_results)
    if args.output:
        all_results.to_csv(args.output)


if __name__ == '__main__':
    main()
