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

"""
Date:   2/10/2015
Author: Peter Krusche <pkrusche@illumina.com>
"""

import pandas
import logging
import re

from Tools.vcfextract import vcfExtract, extractHeadersJSON

def extractMutectSNVFeatures(vcfname, tag, avg_depth=None):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        records = []

        if not avg_depth:
            logging.warn("No average depths available, normalized depth features cannot be calculated")

        hdrs = extractHeadersJSON(vcfname)

        tsn = ""
        nsn = ""

        t_sample = "S.1."
        n_sample = "S.2."

        try:
            samples = hdrs["samples"]
            for f in hdrs["fields"]:
                if f["key"] == "GATKCommandLine" and f["values"]["ID"].lower() == "mutect":
                    clopts = f["values"]["CommandLineOptions"]
                    # ... tumor_sample_name=HCC2218_tumour ... normal_sample_name=HCC2218_normal
                    m = re.search("tumor_sample_name=([^\s]+)", clopts)
                    if m:
                        tsn = m.group(1)
                        for i, x in enumerate(samples):
                            if x == tsn:
                                t_sample = "S.%i." % (i+1)
                                break
                    m = re.search("normal_sample_name=([^\s]+)", clopts)
                    if m:
                        nsn = m.group(1)
                        for i, x in enumerate(samples):
                            if x == nsn:
                                n_sample = "S.%i." % (i+1)
                                break

        except:
            logging.warn("Unable to detect tumour / normal sample order from VCF header")

        logging.info("Normal sample name : %s (prefix %s) / tumour sample name : %s (prefix %s)" % (nsn, n_sample,
                                                                                                    tsn, t_sample))

        features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                    "I.DB", "I.TLOD",
                    n_sample + "GT", t_sample + "GT",
                    n_sample + "DP", t_sample + "DP",
                    n_sample + "AD", t_sample + "AD",
                    n_sample + "BQ", t_sample + "BQ",
                    n_sample + "FA", t_sample + "FA",
                    n_sample + "SS", t_sample + "SS"]

        has_warned = {"feat:I.DB": 1}

        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, ff in enumerate(features):
                rec[ff] = vr[i]

            for q in [n_sample + "GT", t_sample + "GT"]:
                if not q in rec or rec[q] is None:
                    rec[q] = "."
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True

            # fix missing features
            for q in ["I.DB",
                      n_sample + "DP", t_sample + "DP",
                      n_sample + "AD", t_sample + "AD",
                      n_sample + "BQ", t_sample + "BQ",
                      n_sample + "FA", t_sample + "FA",
                      n_sample + "SS", t_sample + "SS"]:
                if not q in rec or rec[q] is None:
                    rec[q] = 0
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True
                else:
                    if q.endswith("FA"):
                        try:
                            rec[q] = float(rec[q])
                        except ValueError:
                            rec[q] = float("NaN")
                    elif q.endswith("AD"):
                        if type(rec[q]) is not list:
                            if not has_warned["AD_PARSE_FAIL"]:
                                logging.warn("Cannot parse AD: %s" % str(rec[q]))
                                has_warned["AD_PARSE_FAIL"] = True
                                rec[q] = [0] * (1 + len(rec["ALT"]))

                            for xx in range(0, 1 + len(rec["ALT"])):
                                if len(rec[q]) <= xx:
                                    rec[q].append(0)
                                else:
                                    try:
                                        rec[q][xx] = float(rec[q][xx])
                                    except ValueError:
                                        rec[q][xx] = 0
                    else:
                        try:
                            rec[q] = int(rec[q])
                        except ValueError:
                            rec[q] = -1

            rec["tag"] = tag
            rec["TLOD"] = float(rec["I.TLOD"])

            n_DP        = float(rec[n_sample + "DP"])
            t_DP        = float(rec[t_sample + "DP"])

            n_DP_ratio = 0
            t_DP_ratio = 0

            if avg_depth:
                if rec["CHROM"] in avg_depth:
                    n_DP_ratio      = n_DP/float(avg_depth[rec["CHROM"]])
                    t_DP_ratio      = t_DP/float(avg_depth[rec["CHROM"]])
                elif not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
            elif not "DPnorm" in has_warned:
                logging.warn("Cannot normalize depths.")
                has_warned["DPnorm"] = True

            n_allele_ref_count = rec[n_sample + "AD"][0]
            alleles_alt = rec["ALT"]

            if alleles_alt == ['.']:
                n_allele_alt_count = 0
            else:
                n_allele_alt_count = 0
                for a in xrange(0, len(alleles_alt)):
                    n_allele_alt_count += float(rec[n_sample + "AD"][a + 1])

            if n_allele_alt_count + n_allele_ref_count == 0:
                n_allele_rate = 0
            else:
                n_allele_rate = n_allele_alt_count / float(n_allele_alt_count + n_allele_ref_count)

            t_allele_ref_count = rec[t_sample + "AD"][0]
            alleles_alt = rec["ALT"]

            if alleles_alt == ['.']:
                t_allele_alt_count = 0
            else:
                t_allele_alt_count = 0
                for a in xrange(0, len(alleles_alt)):
                    t_allele_alt_count += float(rec[t_sample + "AD"][a + 1])

            if t_allele_alt_count + t_allele_ref_count == 0:
                t_allele_rate = 0
            else:
                t_allele_rate = t_allele_alt_count / float(t_allele_alt_count + t_allele_ref_count)

            # Gather the computed data into a dict
            qrec = {
                "CHROM": rec["CHROM"],
                "POS": int(rec["POS"]),
                "REF": rec["REF"],
                "ALT": ",".join(rec["ALT"]),
                "FILTER": ",".join(rec["FILTER"]),
                "DBSNP": rec["I.DB"],
                "N_DP": n_DP,
                "T_DP": t_DP,
                "TLOD": rec["I.TLOD"],
                "N_DP_RATE" : n_DP_ratio,
                "T_DP_RATE" : t_DP_ratio,
                "N_GT": rec[n_sample + "GT"],
                "T_GT": rec[t_sample + "GT"],
                "N_AD": rec[n_sample + "AD"],
                "T_AD": rec[t_sample + "AD"],
                "N_BQ": rec[n_sample + "BQ"],
                "T_BQ": rec[t_sample + "BQ"],
                "N_FA": rec[n_sample + "FA"],
                "T_FA": rec[t_sample + "FA"],
                "N_SS": rec[n_sample + "SS"],
                "T_SS": rec[t_sample + "SS"],
                "N_ALT_RATE": n_allele_rate,
                "T_ALT_RATE": t_allele_rate,
                "tag" : tag
            }
            records.append(qrec)

        cols = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "FILTER",
            "DBSNP",
            "N_DP",
            "T_DP",
            "N_DP_RATE",
            "T_DP_RATE",
            "N_GT",
            "T_GT",
            "N_AD",
            "T_AD",
            "N_BQ",
            "T_BQ",
            "N_FA",
            "T_FA",
            "N_SS",
            "T_SS",
            "TLOD",
            "N_ALT_RATE",
            "T_ALT_RATE",
            "tag"]

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df


def extractMutectIndelFeatures(vcfname, tag, avg_depth=None):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        records = []

        if not avg_depth:
            logging.warn("No average depths available, normalized depth features cannot be calculated")

        hdrs = extractHeadersJSON(vcfname)

        tsn = ""
        nsn = ""

        t_sample = "S.1."
        n_sample = "S.2."

        try:
            samples = hdrs["samples"]
            for f in hdrs["fields"]:
                if f["key"] == "GATKCommandLine" and f["values"]["ID"].lower() == "mutect":
                    clopts = f["values"]["CommandLineOptions"]
                    # ... tumor_sample_name=HCC2218_tumour ... normal_sample_name=HCC2218_normal
                    m = re.search("tumor_sample_name=([^\s]+)", clopts)
                    if m:
                        tsn = m.group(1)
                        for i, x in enumerate(samples):
                            if x == tsn:
                                t_sample = "S.%i." % (i+1)
                                break
                    m = re.search("normal_sample_name=([^\s]+)", clopts)
                    if m:
                        nsn = m.group(1)
                        for i, x in enumerate(samples):
                            if x == nsn:
                                n_sample = "S.%i." % (i+1)
                                break

        except:
            logging.warn("Unable to detect tumour / normal sample order from VCF header")

        logging.info("Normal sample name : %s (prefix %s) / tumour sample name : %s (prefix %s)" % (nsn, n_sample,
                                                                                                    tsn, t_sample))
        has_warned = {}

        ##FORMAT=<ID=MM,Number=2,Type=Float,Description="Average # of mismatches per ref-/consensus indel-supporting read">
        ##FORMAT=<ID=MQS,Number=2,Type=Float,Description="Average mapping qualities of ref-/consensus indel-supporting reads">
        ##FORMAT=<ID=NQSBQ,Number=2,Type=Float,Description="Within NQS window: average quality of bases in ref-/consensus indel-supporting reads">
        ##FORMAT=<ID=NQSMM,Number=2,Type=Float,Description="Within NQS window: fraction of mismatching bases in ref/consensus indel-supporting reads">
        ##FORMAT=<ID=REnd,Number=2,Type=Integer,Description="Median/mad of indel offsets from the ends of the reads">
        ##FORMAT=<ID=RStart,Number=2,Type=Integer,Description="Median/mad of indel offsets from the starts of the reads">
        ##FORMAT=<ID=SC,Number=4,Type=Integer,Description="Strandness: counts of forward-/reverse-aligned reference and indel-supporting reads (FwdRef,RevRef,FwdIndel,RevIndel)">

        features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                    "I.TLOD",
                    n_sample + "GT", t_sample + "GT",
                    n_sample + "DP", t_sample + "DP",
                    n_sample + "AD", t_sample + "AD",
                    n_sample + "MM", t_sample + "MM",
                    n_sample + "MQS", t_sample + "MQS",
                    n_sample + "NQSBQ", t_sample + "NQSBQ",
                    n_sample + "NQSMM", t_sample + "NQSMM",
                    n_sample + "RStart", t_sample + "RStart",
                    n_sample + "REnd", t_sample + "REnd",
                    n_sample + "SC", t_sample + "SC"]

        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, ff in enumerate(features):
                rec[ff] = vr[i]

            for q in [n_sample + "GT", t_sample + "GT"]:
                if not q in rec or rec[q] is None:
                    rec[q] = "."
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True

            # fix missing features
            for q in [n_sample + "GT", t_sample + "GT",
                      n_sample + "DP", t_sample + "DP",
                      n_sample + "AD", t_sample + "AD",
                      n_sample + "MM", t_sample + "MM",
                      n_sample + "MQS", t_sample + "MQS",
                      n_sample + "NQSBQ", t_sample + "NQSBQ",
                      n_sample + "NQSMM", t_sample + "NQSMM",
                      n_sample + "RStart", t_sample + "RStart",
                      n_sample + "REnd", t_sample + "REnd",
                      n_sample + "SC", t_sample + "SC"]:
                if not q in rec or rec[q] is None:
                    rec[q] = 0
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True
                else:
                    if q.endswith("AD") or q.endswith("MM") or q.endswith("MQS") or \
                       q.endswith("NQSBQ") or q.endswith("NQSMM") or \
                       q.endswith("REnd") or q.endswith("RStart"):
                        if type(rec[q]) is not list:
                            if not has_warned[q + "_PARSE_FAIL"]:
                                logging.warn("Cannot parse %s: %s" % (q, str(rec[q])))
                                has_warned[q + "_PARSE_FAIL"] = True
                                rec[q] = [-1, -1]
                            for xx in range(2):
                                if len(rec[q]) <= xx:
                                    rec[q].append(-1)
                                else:
                                    try:
                                        rec[q][xx] = float(rec[q][xx])
                                    except ValueError:
                                        rec[q][xx] = -1
                    elif q.endswith("SC"):
                        if type(rec[q]) is not list:
                            if not has_warned[q + "_PARSE_FAIL"]:
                                logging.warn("Cannot parse %s: %s" % (q, str(rec[q])))
                                has_warned[q + "_PARSE_FAIL"] = True
                                rec[q] = [-1, -1, -1, -1]
                        else:
                            for xx in range(4):
                                if len(rec[q]) <= xx:
                                    rec[q].append(-1)
                                else:
                                    try:
                                        rec[q][xx] = float(rec[q][xx])
                                    except ValueError:
                                        rec[q][xx] = -1
                    else:
                        try:
                            rec[q] = int(rec[q])
                        except ValueError:
                            rec[q] = -1

            rec["tag"] = tag
            rec["TLOD"] = float(rec["I.TLOD"])

            n_DP        = float(rec[n_sample + "DP"])
            t_DP        = float(rec[t_sample + "DP"])

            n_DP_ratio = 0
            t_DP_ratio = 0

            if avg_depth:
                if rec["CHROM"] in avg_depth:
                    n_DP_ratio      = n_DP/float(avg_depth[rec["CHROM"]])
                    t_DP_ratio      = t_DP/float(avg_depth[rec["CHROM"]])
                elif not rec["CHROM"] in has_warned:
                    logging.warn("Cannot normalize depths on %s" % rec["CHROM"])
                    has_warned[rec["CHROM"]] = True
            elif not "DPnorm" in has_warned:
                logging.warn("Cannot normalize depths.")
                has_warned["DPnorm"] = True

            n_allele_ref_count = rec[n_sample + "AD"][0]
            alleles_alt = rec["ALT"]

            if alleles_alt == ['.']:
                n_allele_alt_count = 0
            else:
                n_allele_alt_count = 0
                for a in xrange(1, len(rec[n_sample + "AD"])):
                    n_allele_alt_count += float(rec[n_sample + "AD"][a])

            if n_allele_alt_count + n_allele_ref_count == 0:
                n_allele_rate = 0
            else:
                n_allele_rate = n_allele_alt_count / float(n_allele_alt_count + n_allele_ref_count)

            t_allele_ref_count = rec[t_sample + "AD"][0]
            alleles_alt = rec["ALT"]

            if alleles_alt == ['.']:
                t_allele_alt_count = 0
            else:
                t_allele_alt_count = 0
                for a in xrange(1, len(rec[t_sample + "AD"])):
                    t_allele_alt_count += float(rec[t_sample + "AD"][a])

            if t_allele_alt_count + t_allele_ref_count == 0:
                t_allele_rate = 0
            else:
                t_allele_rate = t_allele_alt_count / float(t_allele_alt_count + t_allele_ref_count)

            # Gather the computed data into a dict
            qrec = {
                "CHROM": rec["CHROM"],
                "POS": int(rec["POS"]),
                "REF": rec["REF"],
                "ALT": ",".join(rec["ALT"]),
                "FILTER": ",".join(rec["FILTER"]),
                "N_DP": n_DP,
                "T_DP": t_DP,
                "N_DP_RATE" : n_DP_ratio,
                "T_DP_RATE" : t_DP_ratio,
                "N_GT": rec[n_sample + "GT"],
                "T_GT": rec[t_sample + "GT"],
                "N_AD": rec[n_sample + "AD"],
                "T_AD": rec[t_sample + "AD"],
                "N_ALT_RATE": n_allele_rate,
                "T_ALT_RATE": t_allele_rate,
                "N_MM": n_sample + "MM",
                "T_MM": t_sample + "MM",
                "N_MQS": n_sample + "MQS",
                "T_MQS": t_sample + "MQS",
                "N_NQSBQ": n_sample + "NQSBQ",
                "T_NQSBQ": t_sample + "NQSBQ",
                "N_NQSMM": n_sample + "NQSMM",
                "T_NQSMM": t_sample + "NQSMM",
                "N_RStart": n_sample + "RStart",
                "T_RStart": t_sample + "RStart",
                "N_REnd": n_sample + "REnd",
                "T_REnd": t_sample + "REnd",
                "N_SC": n_sample + "SC",
                "T_SC": t_sample + "SC",
                "tag" : tag,
                "TLOD": rec["TLOD"]
            }
            records.append(qrec)

        cols = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "FILTER",
            "DBSNP",
            "N_DP",
            "T_DP",
            "N_DP_RATE",
            "T_DP_RATE",
            "N_GT",
            "T_GT",
            "N_AD",
            "T_AD",
            "N_ALT_RATE",
            "T_ALT_RATE",
            "N_MM",
            "T_MM",
            "N_MQS",
            "T_MQS",
            "N_NQSBQ",
            "T_NQSBQ",
            "N_NQSMM",
            "T_NQSMM",
            "N_RStart",
            "T_RStart",
            "N_REnd",
            "T_REnd",
            "N_SC",
            "T_SC",
            "tag",
            "TLOD"]

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df
