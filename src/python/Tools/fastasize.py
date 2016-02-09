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

import os


def fastaContigLengths(fastafile):
    """ Return contig lengths in a fasta file
    """
    if not os.path.exists(fastafile + ".fai"):
        raise Exception("Fasta file %s is not indexed" % fastafile)

    fastacontiglengths = {}

    with open(fastafile + ".fai") as fai:
        for l in fai:
            row = l.strip().split("\t")
            fastacontiglengths[row[0]] = int(row[1])

    return fastacontiglengths


def calculateLength(fastacontiglengths, locations):
    """ Calculate total length of contigs overlapping a set of locations """
    if not locations:
        return sum([fastacontiglengths[x] for x in fastacontiglengths.keys()])

    total_length = 0
    for l in locations.split(" "):
        contig, _, pos = l.partition(":")
        if contig not in fastacontiglengths:
            raise Exception("Contig %s is not present in input set %s" %
                            (contig, str(fastacontiglengths)))

        if pos:
            start, _, end = pos.partition("-")
            start = int(start)
            if end:
                length = int(end) - start + 1
            else:
                length = max(0, fastacontiglengths[contig] - start)
        else:
            length = fastacontiglengths[contig]

        total_length += length

    return total_length
