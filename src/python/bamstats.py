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

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
import Tools.bamstats


def main():
    parser = argparse.ArgumentParser("BAM Statistics")

    parser.add_argument("input", nargs=1)

    parser.add_argument("-o", "--output", dest="output", help="Output CSV file name")

    args = parser.parse_args()
    bs = Tools.bamstats.bamStats(args.input[0])
    print str(bs)
    if args.output:
        bs.to_csv(args.output)


if __name__ == '__main__':
    main()
