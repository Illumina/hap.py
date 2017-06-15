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

import pandas
import logging
from Tools.vcfextract import vcfExtract, extractHeaders


def extractPiscesSNVFeatures(vcfname, tag, avg_depth=None):
    """ Return a data frame with features collected from the given VCF, tagged by given type
    :param vcfname: name of the VCF file
    :param tag: type of variants
    :param avg_depth: average chromosome depths from BAM file
    """
    features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                "I.DP",
                "I.EVS",
                "S.1.GT",
                "S.1.GQ",
                "S.1.AD",
                "S.1.DP",
                "S.1.VF",
                "S.1.NL",
                "S.1.SB",
                "S.1.NC",
                "S.1.AQ",
                "S.1.GQX"]

    cols = ["CHROM", "POS", "REF", "ALT",
            "FILTER", "GQX", "EVS",
            "T_DP", "T_DP_RATE",
            "T_AF",
            "tag"]

    vcfheaders = list(extractHeaders(vcfname))

    evs_featurenames = {}
    for l in vcfheaders:
        if '##snv_scoring_features' in l:
            try:
                xl = str(l).split('=', 1)
                xl = xl[1].split(",")
                for i, n in enumerate(xl):
                    evs_featurenames[i] = n
                    cols.append("E." + n)
                    logging.info("Scoring feature %i : %s" % (i, n))
            except:
                logging.warn("Could not parse scoring feature names from Pisces output")

    records = []

    if not avg_depth:
        avg_depth = {}

        for l in vcfheaders:
            x = str(l).lower()
            x = x.replace("##meandepth_", "##maxdepth_")
            x = x.replace("##depth_", "##maxdepth_")
            if '##maxdepth_' in x:
                p, _, l = l.partition("_")
                xl = str(l).split('=')
                xchr = xl[0]
                avg_depth[xchr] = float(xl[1])
                logging.info("%s depth from VCF header is %f" % (xchr, avg_depth[xchr]))

    has_warned = {}

    for vr in vcfExtract(vcfname, features):
        rec = {}
        for i, ff in enumerate(features):
            rec[ff] = vr[i]

        # read VQSR value, if it's not present, set to -1 (old versions of Pisces)
        try:
            rec["I.VQSR"] = float(rec["I.VQSR"])
        except:
            rec["I.VQSR"] = -1.0

        # read EVS value, if it's not present, set to -1 (old versions of Pisces)
        if "I.SomaticEVS" in rec:
            try:
                rec["I.EVS"] = float(rec["I.SomaticEVS"])
            except:
                rec["I.EVS"] = -1.0
        else:
            try:
                rec["I.EVS"] = float(rec["I.EVS"])
            except:
                rec["I.EVS"] = -1.0

        # fix missing features
        for q in ["S.1.NC", "S.1.AQ"]:
            if q not in rec or rec[q] is None:
                rec[q] = 0
                if not ("feat:" + q) in has_warned:
                    logging.warn("Missing feature %s" % q)
                    has_warned["feat:" + q] = True

        rec["tag"] = tag

        t_DP = float(rec["S.1.DP"])
        t_VF = float(rec["S.1.VF"])
        GQX = float(rec["S.1.GQX"])

        t_DP_ratio = 0

        if avg_depth:
            try:
                t_DP_ratio = t_DP / float(avg_depth[rec["CHROM"]])
            except:
                if not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
        elif "DPnorm" not in has_warned:
            logging.warn("Cannot normalize depths.")
            has_warned["DPnorm"] = True

        # Gather the computed data into a dict
        qrec = {
            "CHROM": rec["CHROM"],
            "POS": int(rec["POS"]),
            "REF": rec["REF"],
            "ALT": ",".join(rec["ALT"]),
            "FILTER": ",".join(rec["FILTER"]),
            "GQX": GQX,
            "EVS": rec["I.EVS"],
            "T_DP": t_DP,
            "T_DP_RATE": t_DP_ratio,
            "T_AF": t_VF,
            "tag": tag
        }

        records.append(qrec)

    if records:
        df = pandas.DataFrame(records, columns=cols)
    else:
        df = pandas.DataFrame(columns=cols)

    return df


def extractPiscesIndelFeatures(vcfname, tag, avg_depth=None):
    """ Return a data frame with features collected from the given VCF, tagged by given type
    :param vcfname: name of the VCF file
    :param tag: type of variants
    :param avg_depth: average chromosome depths from BAM file
    """

    return extractPiscesSNVFeatures(vcfname, tag, avg_depth)
