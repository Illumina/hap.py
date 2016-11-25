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

import pysam
import pandas
import numpy as np
import logging


def bamStats(bamfile):
    """ Extract average depths + idxstats data from BAM file, return data frame """
    istats = pysam.idxstats(bamfile)
    result = []
    samfile = pysam.Samfile(bamfile, "rb")

    for x in istats.split("\n"):
        xs = x.replace("\n", "").split("\t")
        if len(xs) < 4:
            logging.warn("Ignoring invalid stats line: %s" % x)
            continue
        rec = {
            "CHROM": xs[0],
            "NT": int(xs[1]),
            "MAPPED": int(xs[2]),
            "UNMAPPED": int(xs[3]),
            "READLEN": 0,
            "COVERAGE": 0.0,
        }

        count = 0
        rls = 0.0
        try:
            for read in samfile.fetch(xs[0]):
                rls += float(read.rlen)
                count += 1
                if count > 10000:
                    break

            rls /= count
            rec["READLEN"] = rls
            rec["COVERAGE"] = float(rec["MAPPED"] * rec["READLEN"])/float(rec["NT"])
        except:
            pass
        result.append(rec)

    if result:
        result = pandas.DataFrame(result, columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])
        min_readlen = np.min(result[result["READLEN"] > 0]["READLEN"])
        max_readlen = np.max(result[result["READLEN"] > 0]["READLEN"])
        if min_readlen != max_readlen:
            logging.warn("Read lengths differ within the same BAM file: %s" % str(result["READLEN"].unique()))

        agg_result = [{
            "CHROM": "TOTAL",
            "NT": np.sum(result["NT"]),
            "MAPPED": np.sum(result["MAPPED"]),
            "UNMAPPED": np.sum(result["UNMAPPED"]),
            "READLEN": (max_readlen + min_readlen) / 2.0,
        }]
        agg_result[-1]["COVERAGE"] = np.sum(result["MAPPED"].multiply(result["READLEN"])) / np.sum(result["NT"])

        auto_result = result[result["CHROM"].str.match(r"^(?:chr)?[0-9]+$")]
        min_readlen = np.min(auto_result[auto_result["READLEN"] > 0]["READLEN"])
        max_readlen = np.max(auto_result[auto_result["READLEN"] > 0]["READLEN"])
        agg_result.append({
            "CHROM": "AUTOSOME",
            "NT": np.sum(auto_result["NT"]),
            "MAPPED": np.sum(auto_result["MAPPED"]),
            "UNMAPPED": np.sum(auto_result["UNMAPPED"]),
            "READLEN": (max_readlen + min_readlen) / 2.0,
        })
        agg_result[-1]["COVERAGE"] = np.sum(auto_result["MAPPED"].multiply(auto_result["READLEN"])) / np.sum(auto_result["NT"])

        return pandas.concat([result, pandas.DataFrame(agg_result, columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])]).set_index(["CHROM"])
    else:
        return pandas.DataFrame(columns=["CHROM", "NT", "MAPPED", "UNMAPPED", "READLEN", "COVERAGE"])
