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


def extractStrelkaSNVFeatures(vcfname, tag, avg_depth=None):
    """ Return a data frame with features collected from the given VCF, tagged by given type
    :param vcfname: name of the VCF file
    :param tag: type of variants
    :param avg_depth: average chromosome depths from BAM file
    """
    features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                "I.NT", "I.SOMATIC", "I.QSS_NT",
                "I.VQSR", "I.EVS", "I.EVSF", "I.SomaticEVS",
                "I.SGT", "I.MQ", "I.MQ0",
                "I.SNVSB", "I.ReadPosRankSum",
                "S.1.SDP", "S.2.SDP",
                "S.1.FDP", "S.2.FDP",
                "S.1.DP", "S.2.DP",
                "S.1.AU", "S.2.AU",
                "S.1.CU", "S.2.CU",
                "S.1.GU", "S.2.GU",
                "S.1.TU", "S.2.TU"]

    cols = ["CHROM", "POS", "REF", "ALT",
            "NT", "NT_REF", "QSS_NT", "FILTER", "SomaticEVS", "EVS", "VQSR",
            "N_FDP_RATE", "T_FDP_RATE", "N_SDP_RATE", "T_SDP_RATE",
            "N_DP", "T_DP", "N_DP_RATE", "T_DP_RATE",
            "N_AF", "T_AF",
            "MQ", "MQ0",
            "SNVSB",
            "ReadPosRankSum", "tag"]

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
                logging.warn("Could not parse scoring feature names from Strelka output")

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

        # read VQSR value, if it's not present, set to -1 (old versions of Strelka)
        try:
            rec["I.VQSR"] = float(rec["I.VQSR"])
        except:
            rec["I.VQSR"] = -1.0

        # read EVS value, if it's not present, set to -1 (old versions of Strelka)
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
        for q in ["I.QSS_NT", "I.MQ", "I.MQ0",
                  "I.SNVSB", "I.ReadPosRankSum", "S.1.SDP", "S.2.SDP",
                  "S.1.FDP", "S.2.FDP",
                  "S.1.DP", "S.2.DP",
                  "S.1.AU", "S.2.AU",
                  "S.1.CU", "S.2.CU",
                  "S.1.GU", "S.2.GU",
                  "S.1.TU", "S.2.TU"]:
            if q not in rec or rec[q] is None:
                rec[q] = 0
                if not ("feat:" + q) in has_warned:
                    logging.warn("Missing feature %s" % q)
                    has_warned["feat:" + q] = True

        rec["tag"] = tag

        NT = rec["I.NT"]
        NT_is_ref = int(NT == "ref")
        QSS_NT = int(rec["I.QSS_NT"])

        try:
            MQ = float(rec["I.MQ"])
        except:
            MQ = None

        try:
            MQ_ZERO = float(rec["I.MQ0"])
        except:
            MQ_ZERO = None

        n_FDP = float(rec["S.1.FDP"])
        t_FDP = float(rec["S.2.FDP"])
        n_SDP = float(rec["S.1.SDP"])
        t_SDP = float(rec["S.2.SDP"])
        n_DP = float(rec["S.1.DP"])
        t_DP = float(rec["S.2.DP"])

        n_FDP_ratio = n_FDP / n_DP if n_DP != 0 else 0
        t_FDP_ratio = t_FDP / t_DP if t_DP != 0 else 0

        n_SDP_ratio = n_SDP / (n_DP + n_SDP) if (n_DP + n_SDP) != 0 else 0
        t_SDP_ratio = t_SDP / (t_DP + t_SDP) if (t_DP + t_SDP) != 0 else 0

        n_DP_ratio = 0
        t_DP_ratio = 0

        if avg_depth:
            try:
                n_DP_ratio = n_DP / float(avg_depth[rec["CHROM"]])
                t_DP_ratio = t_DP / float(avg_depth[rec["CHROM"]])
            except:
                if not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
        elif "DPnorm" not in has_warned:
            logging.warn("Cannot normalize depths.")
            has_warned["DPnorm"] = True

        # Ref and alt allele counts for tier1 and tier2
        allele_ref = rec["REF"]
        try:
            t_allele_ref_counts = map(float, rec['S.2.' + allele_ref + 'U'])
        except:
            t_allele_ref_counts = [0, 0]

        alleles_alt = rec["ALT"]

        try:
            t_allele_alt_counts = [0, 0]
            for a in alleles_alt:
                for i in range(2):
                    t_allele_alt_counts[i] += float(rec['S.2.' + a + 'U'][i])
        except:
            t_allele_alt_counts = [0, 0]

        # Compute the tier1 and tier2 alt allele rates.
        if t_allele_alt_counts[0] + t_allele_ref_counts[0] == 0:
            t_tier1_allele_rate = 0
        else:
            t_tier1_allele_rate = t_allele_alt_counts[0] / float(t_allele_alt_counts[0] + t_allele_ref_counts[0])

        try:
            n_allele_ref_counts = map(float, rec['S.1.' + allele_ref + 'U'])
        except:
            n_allele_ref_counts = [0, 0]

        alleles_alt = rec["ALT"]

        try:
            n_allele_alt_counts = [0, 0]
            for a in alleles_alt:
                for i in range(2):
                    n_allele_alt_counts[i] += float(rec['S.1.' + a + 'U'][i])
        except:
            n_allele_alt_counts = [0, 0]

        # Compute the tier1 and tier2 alt allele rates.
        if n_allele_alt_counts[0] + n_allele_ref_counts[0] == 0:
            n_tier1_allele_rate = 0
        else:
            n_tier1_allele_rate = n_allele_alt_counts[0] / float(n_allele_alt_counts[0] + n_allele_ref_counts[0])

        try:
            snvsb = rec["I.SNVSB"]
        except:
            snvsb = 0

        try:
            rprs = rec["I.ReadPosRankSum"]
        except:
            rprs = 0

        # Gather the computed data into a dict
        qrec = {
            "CHROM": rec["CHROM"],
            "POS": int(rec["POS"]),
            "REF": rec["REF"],
            "ALT": ",".join(rec["ALT"]),
            "FILTER": ",".join(rec["FILTER"]),
            "NT": NT,
            "NT_REF": NT_is_ref,
            "QSS_NT": QSS_NT,
            "VQSR": rec["I.VQSR"],
            "EVS": rec["I.EVS"],
            "N_FDP_RATE": n_FDP_ratio,
            "T_FDP_RATE": t_FDP_ratio,
            "N_SDP_RATE": n_SDP_ratio,
            "T_SDP_RATE": t_SDP_ratio,
            "N_DP": n_DP,
            "T_DP": t_DP,
            "N_DP_RATE": n_DP_ratio,
            "T_DP_RATE": t_DP_ratio,
            "N_AF": n_tier1_allele_rate,
            "T_AF": t_tier1_allele_rate,
            "MQ": MQ,
            "MQ0": MQ_ZERO,
            "SNVSB": snvsb,
            "ReadPosRankSum": rprs,
            "tag": tag
        }
        # ESF features
        try:
            for i, v in enumerate(rec["I.EVSF"]):
                if i in evs_featurenames:
                    try:
                        qrec["E." + evs_featurenames[i]] = float(v)
                    except:
                        # failure to parse
                        pass
        except:
            pass
        for k, v in evs_featurenames.iteritems():
            if not "E." + v in qrec:
                qrec["E." + v] = 0

        records.append(qrec)

    if records:
        df = pandas.DataFrame(records, columns=cols)
    else:
        df = pandas.DataFrame(columns=cols)

    return df


def extractStrelkaIndelFeatures(vcfname, tag, avg_depth=None):
    """ Return a data frame with features collected from the given VCF, tagged by given type
    :param vcfname: name of the VCF file
    :param tag: type of variants
    :param avg_depth: average chromosome depths from BAM file
    """
    features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                "I.NT", "I.SOMATIC", "I.QSI_NT", "I.EVS", "I.EVSF", "I.SomaticEVS",
                "I.SGT", "I.RC", "I.RU",
                "I.IC", "I.IHP",
                "I.MQ", "I.MQ0",
                "S.1.DP", "S.2.DP",
                "S.1.TAR", "S.2.TAR",
                "S.1.TIR", "S.2.TIR",
                "S.1.TOR", "S.2.TOR",
                "S.1.BCN50", "S.2.BCN50",
                "S.1.FDP50", "S.2.FDP50",
                ]

    cols = ["CHROM",
            "POS",
            "REF",
            "ALT",
            "LENGTH",
            "INDELTYPE",
            "FILTER",
            "NT",
            "NT_REF",
            "EVS",
            "QSI_NT",
            "N_DP",
            "T_DP",
            "N_DP_RATE",
            "T_DP_RATE",
            "N_BCN",
            "T_BCN",
            "N_FDP",
            "T_FDP",
            "N_AF",
            "T_AF",
            "SGT",
            "RC",
            "RU",
            "RU_LEN",
            "IC",
            "IHP",
            "MQ",
            "MQ0",
            "tag"]

    records = []

    vcfheaders = list(extractHeaders(vcfname))

    evs_featurenames = {}

    for l in vcfheaders:
        if '##indel_scoring_features' in l:
            try:
                xl = str(l).split('=', 1)
                xl = xl[1].split(",")
                for i, n in enumerate(xl):
                    evs_featurenames[i] = n
                    cols.append("E." + n)
                    logging.info("Scoring feature %i : %s" % (i, n))
            except:
                logging.warn("Could not parse scoring feature names from Strelka output")

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
        rec["tag"] = tag

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
        for q in ["I.QSI_NT", "I.RC", "I.IC", "I.IHP",
                  "S.1.DP", "S.2.DP",
                  "S.1.BCN50", "S.2.BCN50",
                  "S.1.FDP50", "S.2.FDP50"]:
            if q not in rec or rec[q] is None:
                rec[q] = 0
                if not ("feat:" + q) in has_warned:
                    logging.warn("Missing feature %s" % q)
                    has_warned["feat:" + q] = True

        for q in ["S.1.TAR", "S.2.TAR",
                  "S.1.TIR", "S.2.TIR",
                  "S.1.TOR", "S.2.TOR"]:
            if q not in rec or rec[q] is None:
                rec[q] = [0, 0]
                if not ("feat:" + q) in has_warned:
                    logging.warn("Missing feature %s" % q)
                    has_warned["feat:" + q] = True

        NT = rec["I.NT"]
        NT_is_ref = int(NT == "ref")
        QSI_NT = int(rec["I.QSI_NT"])

        n_DP = float(rec["S.1.DP"])
        t_DP = float(rec["S.2.DP"])

        in_del = 0

        max_len = len(rec["REF"])
        min_len = len(rec["REF"])

        for a in rec["ALT"]:
            if len(a) > len(rec["REF"]):
                in_del |= 1
            else:
                in_del |= 2
            min_len = min(len(a), min_len)
            max_len = max(len(a), max_len)

        ilen = max_len - min_len

        n_DP_ratio = 0
        t_DP_ratio = 0

        if avg_depth:
            try:
                n_DP_ratio = n_DP / float(avg_depth[rec["CHROM"]])
                t_DP_ratio = t_DP / float(avg_depth[rec["CHROM"]])
            except:
                if not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
        elif "DPnorm" not in has_warned:
            logging.warn("Cannot normalize depths.")
            has_warned["DPnorm"] = True

        # extract observed AF from strelka counts. TIR = ALT; TAR = REF
        try:
            n_af = float(rec["S.1.TIR"][0]) / (float(rec["S.1.TIR"][0]) + float(rec["S.1.TAR"][0]))
        except:
            n_af = 0

        try:
            t_af = float(rec["S.2.TIR"][0]) / (float(rec["S.2.TIR"][0]) + float(rec["S.2.TAR"][0]))
        except:
            t_af = 0

        # Gather the computed data into a dict
        qrec = {
            "CHROM": rec["CHROM"],
            "POS": int(rec["POS"]),
            "REF": rec["REF"],
            "ALT": ",".join(rec["ALT"]),
            "LENGTH": ilen,
            "INDELTYPE": in_del,
            "FILTER": ",".join(rec["FILTER"]),
            "NT": NT,
            "NT_REF": NT_is_ref,
            "QSI_NT": QSI_NT,
            "N_DP": n_DP,
            "T_DP": t_DP,
            "N_DP_RATE": n_DP_ratio,
            "T_DP_RATE": t_DP_ratio,
            "N_AF": n_af,
            "T_AF": t_af,
            "SGT": rec["I.SGT"],
            "tag": tag
        }

        # fields with defaults
        fields = [
            {"n": "EVS", "s": "I.EVS", "def": 0, "t": float},
            {"n": "VQSR", "s": "I.VQSR", "def": 0, "t": float},
            {"n": "RC", "s": "I.RC", "def": 0, "t": int},
            {"n": "RU", "s": "I.RU", "def": ""},
            {"n": "RU_LEN", "s": "I.RU", "def": 0, "t": len},
            {"n": "IC", "s": "I.IC", "def": 0, "t": int},
            {"n": "IHP", "s": "I.IHP", "def": 0, "t": int},
            {"n": "MQ", "s": "I.MQ", "def": 0.0, "t": float},
            {"n": "MQ0", "s": "I.MQ0", "def": 0.0, "t": float},
            {"n": "N_BCN", "s": "S.1.BCN50", "def": 0.0, "t": float},
            {"n": "T_BCN", "s": "S.2.BCN50", "def": 0.0, "t": float},
            {"n": "N_FDP", "s": "S.1.FDP50", "def": 0.0, "t": float},
            {"n": "T_FDP", "s": "S.2.FDP50", "def": 0.0, "t": float},
        ]

        for fd in fields:
            try:
                res = rec[fd["s"]]
                if "t" in fd:
                    res = fd["t"](res)
            except:
                res = fd["def"]

            qrec[fd["n"]] = res

        # ESF features
        try:
            for i, v in enumerate(rec["I.EVSF"]):
                if i in evs_featurenames:
                    try:
                        qrec["E." + evs_featurenames[i]] = float(v)
                    except:
                        # failure to parse
                        pass
        except:
            pass

        for k, v in evs_featurenames.iteritems():
            if not "E." + v in qrec:
                qrec["E." + v] = 0

        records.append(qrec)

    if records:
        df = pandas.DataFrame(records, columns=cols)
    else:
        df = pandas.DataFrame(columns=cols)

    return df
