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
import logging
import tempfile
import subprocess
import pipes
import json
import re


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


def fastaNonNContigLengths(fastafile):
    """ Return contig lengths in a fasta file excluding
        Ns in the beginning or end
    """
    if not os.path.exists(fastafile + ".fai"):
        raise Exception("Fasta file %s is not indexed" % fastafile)

    fastacontiglengths = {}

    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()

    try:
        subprocess.check_call("fastainfo %s %s" % (pipes.quote(fastafile), pipes.quote(tf.name)), shell=True)
        with open(tf.name) as f:
            fasta_info = json.load(f)

        for k, v in fasta_info.iteritems():
            fastacontiglengths[k] = int(v["n_trimmed_length"])
    finally:
        os.unlink(tf.name)

    return fastacontiglengths


def calculateLength(fastacontiglengths, locations):
    """ Calculate total length of contigs overlapping a set of locations """
    if not locations:
        return sum([fastacontiglengths[x] for x in fastacontiglengths.keys()])

    total_length = 0
    for l in re.split("[ ,]", locations):
        contig, _, pos = l.partition(":")
        if contig not in fastacontiglengths:
            logging.warn("Contig %s is not present in input set %s; setting length to 0" %
                            (contig, str(fastacontiglengths)))
            length = 0
            continue

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
