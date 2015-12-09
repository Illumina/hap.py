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

import pandas
import logging

from Tools.vcfextract import vcfExtract, extractHeaders


def extractStrelkaSNVFeatures(vcfname, tag, avg_depth=None):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                    "I.NT", "I.SOMATIC", "I.QSS_NT", "I.VQSR",
                    "I.SGT", "I.MQ", "I.MQ0", "I.PNOISE", "I.PNOISE2",
                    "I.SNVSB", "I.ReadPosRankSum",
                    "S.1.SDP", "S.2.SDP",
                    "S.1.FDP", "S.2.FDP",
                    "S.1.DP", "S.2.DP",
                    "S.1.AU", "S.2.AU",
                    "S.1.CU", "S.2.CU",
                    "S.1.GU", "S.2.GU",
                    "S.1.TU", "S.2.TU"]

        records = []

        if not avg_depth:
            avg_depth = {}

            for l in list(extractHeaders(vcfname)):
                x = str(l).lower()
                x = x.replace("meandepth_", "maxdepth_")
                if '##maxdepth_' in x:
                    xl = str(l).split('=')
                    xchr = xl[0][12:]
                    avg_depth[xchr] = float(xl[1])
                    logging.info("%s depth from VCF header is %f" % (xchr, avg_depth[xchr]))

        has_warned = {}

        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, ff in enumerate(features):
                rec[ff] = vr[i]

            # fix missing features
            for q in ["I.QSS_NT", "I.MQ", "I.MQ0", "I.PNOISE", "I.PNOISE2", "I.VQSR",
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

            n_FDP_ratio = n_FDP/n_DP if n_DP != 0 else 0
            t_FDP_ratio = t_FDP/t_DP if t_DP != 0 else 0

            n_SDP_ratio = n_SDP/(n_DP + n_SDP) if (n_DP + n_SDP) != 0 else 0
            t_SDP_ratio = t_SDP/(t_DP + t_SDP) if (t_DP + t_SDP) != 0 else 0

            n_DP_ratio = 0
            t_DP_ratio = 0

            if avg_depth:
                if rec["CHROM"] in avg_depth:
                    n_DP_ratio = n_DP/float(avg_depth[rec["CHROM"]])
                    t_DP_ratio = t_DP/float(avg_depth[rec["CHROM"]])
                elif not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
            elif "DPnorm" not in has_warned:
                logging.warn("Cannot normalize depths.")
                has_warned["DPnorm"] = True

            # Ref and alt allele counts for tier1 and tier2
            allele_ref = rec["REF"]
            t_allele_ref_counts = map(float, rec['S.2.' + allele_ref + 'U'])

            alleles_alt = rec["ALT"]

            if alleles_alt == ['.']:
                t_allele_alt_counts = [0, 0]
            else:
                t_allele_alt_counts = [0, 0]
                for a in alleles_alt:
                    for i in range(2):
                        t_allele_alt_counts[i] += float(rec['S.2.' + a + 'U'][i])

            # Compute the tier1 and tier2 alt allele rates.
            if t_allele_alt_counts[0] + t_allele_ref_counts[0] == 0:
                t_tier1_allele_rate = 0
            else:
                t_tier1_allele_rate = t_allele_alt_counts[0] / float(t_allele_alt_counts[0] + t_allele_ref_counts[0])

            if t_allele_alt_counts[1] + t_allele_ref_counts[1] == 0:
                t_tier2_allele_rate = 0
            else:
                t_tier2_allele_rate = t_allele_alt_counts[1] / float(t_allele_alt_counts[1] + t_allele_ref_counts[1])

            n_allele_ref_counts = map(float, rec['S.1.' + allele_ref + 'U'])

            alleles_alt = rec["ALT"]

            if alleles_alt == ['.']:
                n_allele_alt_counts = [0, 0]
            else:
                n_allele_alt_counts = [0, 0]
                for a in alleles_alt:
                    for i in range(2):
                        n_allele_alt_counts[i] += float(rec['S.1.' + a + 'U'][i])

            # Compute the tier1 and tier2 alt allele rates.
            if n_allele_alt_counts[0] + n_allele_ref_counts[0] == 0:
                n_tier1_allele_rate = 0
            else:
                n_tier1_allele_rate = n_allele_alt_counts[0] / float(n_allele_alt_counts[0] + n_allele_ref_counts[0])

            if n_allele_alt_counts[1] + n_allele_ref_counts[1] == 0:
                n_tier2_allele_rate = 0
            else:
                n_tier2_allele_rate = n_allele_alt_counts[1] / float(n_allele_alt_counts[1] + n_allele_ref_counts[1])

            try:
                pnoise = rec["I.PNOISE"]
            except:
                pnoise = 0

            try:
                pnoise2 = rec["I.PNOISE2"]
            except:
                pnoise2 = 0

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
                "N_FDP_RATE": n_FDP_ratio,
                "T_FDP_RATE": t_FDP_ratio,
                "N_SDP_RATE": n_SDP_ratio,
                "T_SDP_RATE": t_SDP_ratio,
                "N_DP": n_DP,
                "T_DP": t_DP,
                "N_DP_RATE": n_DP_ratio,
                "T_DP_RATE": t_DP_ratio,
                "T_TIER1_ALT_RATE": t_tier1_allele_rate,
                "T_TIER2_ALT_RATE": t_tier2_allele_rate,
                "N_TIER1_ALT_RATE": n_tier1_allele_rate,
                "N_TIER2_ALT_RATE": n_tier2_allele_rate,
                "MQ_SCORE": MQ,
                "MQ_ZERO_RATE": MQ_ZERO,
                "PNOISE": pnoise,
                "PNOISE2": pnoise2,
                "SNVSB": snvsb,
                "ReadPosRankSum": rprs,
                "tag": tag
            }
            records.append(qrec)

        cols = ["CHROM", "POS", "REF", "ALT",
                "NT", "NT_REF", "QSS_NT", "FILTER", "VQSR",
                "N_FDP_RATE", "T_FDP_RATE", "N_SDP_RATE", "T_SDP_RATE",
                "N_DP", "T_DP", "N_DP_RATE", "T_DP_RATE",
                "T_TIER1_ALT_RATE", "T_TIER2_ALT_RATE", "N_TIER1_ALT_RATE", "N_TIER2_ALT_RATE",
                "MQ_SCORE", "MQ_ZERO_RATE", "PNOISE", "PNOISE2", "SNVSB",
                "ReadPosRankSum", "tag"]

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df


def extractStrelkaIndelFeatures(vcfname, tag, avg_depth=None):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                    "I.NT", "I.SOMATIC", "I.QSI_NT", "I.EQSI", "I.ESF",
                    "I.SGT", "I.RC", "I.RU",
                    "I.IC", "I.IHP",
                    "I.MQ", "I.MQ0",
                    "S.1.DP", "S.2.DP",
                    "S.1.TAR", "S.2.TAR",
                    "S.1.TIR", "S.2.TIR",
                    "S.1.TOR", "S.2.TOR",
                    "S.1.AF", "S.2.AF",
                    "S.1.OF", "S.2.OF",
                    "S.1.SOR", "S.2.SOR",
                    "S.1.FS", "S.2.FS",
                    "S.1.BSA", "S.2.BSA",
                    "S.1.RR", "S.2.RR",
                    "S.1.BCN50", "S.2.BCN50",
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
                "VQSR",
                "EQSI",
                "QSI_NT",
                "N_DP",
                "T_DP",
                "N_DP_RATE",
                "T_DP_RATE",
                "N_AF",
                "T_AF",
                "N_OF",
                "T_OF",
                "N_SOR",
                "T_SOR",
                "N_FS",
                "T_FS",
                "N_BSA",
                "T_BSA",
                "N_RR",
                "T_RR",
                "N_BCN",
                "T_BCN",
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

        vqsr_featurenames = {}

        for l in vcfheaders:
            if '##vqsr_features' in l:
                try:
                    xl = str(l).split('=', 1)
                    xl = xl[1].split(",")
                    for x in xl:
                        i, n = x.split(":", 1)
                        i = int(i)
                        vqsr_featurenames[i] = n
                        cols.append("VQSR." + n)
                        logging.info("VQSR feature %i : %s" % (i, n))
                except:
                    logging.warn("Could not parse VQSR feature names from Strelka output")

        if not avg_depth:
            avg_depth = {}

            for l in vcfheaders:
                x = str(l).lower()
                x = x.replace("meandepth_", "maxdepth_")
                if '##maxdepth_' in x:
                    xl = str(l).split('=')
                    xchr = xl[0][12:]
                    avg_depth[xchr] = float(xl[1])
                    logging.info("%s depth from VCF header is %f" % (xchr, avg_depth[xchr]))

        has_warned = {}
        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, ff in enumerate(features):
                rec[ff] = vr[i]
            rec["tag"] = tag

            # fix missing features
            for q in ["I.QSI_NT", "I.RC", "I.IC", "I.IHP", "I.EQSI",
                      "S.1.DP", "S.2.DP",
                      "S.1.OF", "S.2.OF",
                      "S.1.RR", "S.2.RR",
                      "S.1.FS", "S.2.FS",
                      "S.1.BSA", "S.2.BSA",
                      "S.1.BCN50", "S.2.BCN50",
                      "S.1.AF", "S.2.AF" ]:
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
                if rec["CHROM"] in avg_depth:
                    n_DP_ratio = n_DP/float(avg_depth[rec["CHROM"]])
                    t_DP_ratio = t_DP/float(avg_depth[rec["CHROM"]])
                elif not rec["CHROM"] in has_warned:
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
                "SGT": rec["I.SGT"],
                "tag": tag
            }

            # fields with defaults
            fields = [
                { "n": "EQSI", "s": "I.EQSI", "def": 0, "t": float },
                { "n": "VQSR", "s": "I.EQSI", "def": 0, "t": float },
                { "n": "RC", "s": "I.RC", "def": 0, "t": int },
                { "n": "RU", "s": "I.RU", "def": "" },
                { "n": "RU_LEN", "s": "I.RU", "def": 0, "t": len },
                { "n": "IC", "s": "I.IC", "def": 0, "t": int },
                { "n": "IHP", "s": "I.IHP", "def": 0, "t": int },
                { "n": "MQ", "s": "I.MQ", "def": 0.0, "t": float },
                { "n": "MQ0", "s": "I.MQ0", "def": 0.0, "t": float },
                { "n": "N_AF", "s": "S.1.AF", "def": 0.0, "t": float },
                { "n": "T_AF", "s": "S.2.AF", "def": 0.0, "t": float },
                { "n": "N_OF", "s": "S.1.OF", "def": 0.0, "t": float },
                { "n": "T_OF", "s": "S.2.OF", "def": 0.0, "t": float },
                { "n": "N_SOR", "s": "S.1.SOR", "def": 0.0, "t": float },
                { "n": "T_SOR", "s": "S.2.SOR", "def": 0.0, "t": float },
                { "n": "N_FS", "s": "S.1.FS", "def": 0.0, "t": float },
                { "n": "T_FS", "s": "S.2.FS", "def": 0.0, "t": float },
                { "n": "N_BSA", "s": "S.1.BSA", "def": 0.0, "t": float },
                { "n": "T_BSA", "s": "S.2.BSA", "def": 0.0, "t": float },
                { "n": "N_RR", "s": "S.1.RR", "def": 0.0, "t": float },
                { "n": "T_RR", "s": "S.2.RR", "def": 0.0, "t": float },
                { "n": "N_BCN", "s": "S.1.BCN50", "def": 0.0, "t": float },
                { "n": "T_BCN", "s": "S.2.BCN50", "def": 0.0, "t": float },
            ]

            for fd in fields:
                try:
                    res = rec[fd["s"]]
                    if "t" in fd:
                        res = fd["t"](res)
                except:
                    res = fd["def"]

                qrec[fd["n"]] = res

            # VQSR features
            try:
                for i, v in enumerate(rec["I.ESF"]):
                    if i in vqsr_featurenames:
                        try:
                            qrec["VQSR." + vqsr_featurenames[i]] = float(v)
                        except:
                            # failure to parse
                            pass
            except:
                pass
            for k, v in vqsr_featurenames.iteritems():
                if not "VQSR." + v in qrec:
                    qrec["VQSR." + v] = 0

            records.append(qrec)

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df
