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
                    "I.DB", "I.TLOD", "I.NLOD", "I.ECNT",
                    "I.HCNT", "I.MAX_ED", "I.MIN_ED",
                    n_sample + "GT", t_sample + "GT",
                    n_sample + "DP", t_sample + "DP",
                    n_sample + "QSS", t_sample + "QSS",
                    n_sample + "AD", t_sample + "AD"]

        has_warned = {}

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
            for q in ["I.DB", "I.TLOD", "I.NLOD", "I.ECNT",
                      "I.HCNT", "I.MAX_ED", "I.MIN_ED",
                      n_sample + "GT", t_sample + "GT",
                      n_sample + "DP", t_sample + "DP",
                      n_sample + "QSS", t_sample + "QSS",
                      n_sample + "AD", t_sample + "AD"]:
                if not q in rec or rec[q] is None:
                    rec[q] = 0
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True
                else:
                    # list features
                    if q.endswith("AD") or q.endswith("QSS"):
                        if type(rec[q]) is not list:
                            if not q + "_PARSE_FAIL" in has_warned:
                                logging.warn("Cannot parse %s: %s" % (q, str(rec[q])))
                                has_warned[q + "_PARSE_FAIL"] = True
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
            TLOD = float(rec["I.TLOD"])
            NLOD = float(rec["I.NLOD"])

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
                "TLOD": TLOD,
                "NLOD": NLOD,
                "N_DP": n_DP,
                "T_DP": t_DP,
                "N_DP_RATE" : n_DP_ratio,
                "T_DP_RATE" : t_DP_ratio,
                "N_GT": rec[n_sample + "GT"],
                "T_GT": rec[t_sample + "GT"],
                "N_AD": rec[n_sample + "AD"],
                "T_AD": rec[t_sample + "AD"],
                "N_QSS": rec[n_sample + "QSS"],
                "T_QSS": rec[t_sample + "QSS"],
                "N_AF": n_allele_rate,
                "T_AF": t_allele_rate,
                "ECNT": rec["I.ECNT"],
                "HCNT": rec["I.HCNT"],
                "MAX_ED": rec["I.MAX_ED"],
                "MIN_ED": rec["I.MIN_ED"],
                "tag" : tag
            }
            records.append(qrec)

        cols = ["CHROM", "POS", "REF", "ALT",
                "FILTER", "TLOD", "NLOD", "DBSNP",
                "N_DP", "T_DP", "N_DP_RATE", "T_DP_RATE", "N_GT", "T_GT",
                "N_AD", "T_AD", "N_QSS", "T_QSS",
                "N_AF", "T_AF",
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

        features = ["CHROM", "POS", "REF", "ALT", "FILTER",
                    "I.DB", "I.TLOD", "I.NLOD", "I.ECNT",
                    "I.HCNT", "I.MAX_ED", "I.MIN_ED",
                    "I.RPA", "I.RU", # indel only
                    n_sample + "GT", t_sample + "GT",
                    n_sample + "DP", t_sample + "DP",
                    n_sample + "QSS", t_sample + "QSS",
                    n_sample + "AD", t_sample + "AD"]

        has_warned = {}

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
            for q in ["I.DB", "I.TLOD", "I.NLOD", "I.ECNT",
                      "I.HCNT", "I.MAX_ED", "I.MIN_ED",
                      "I.RPA", "I.RU",
                      n_sample + "GT", t_sample + "GT",
                      n_sample + "DP", t_sample + "DP",
                      n_sample + "QSS", t_sample + "QSS",
                      n_sample + "AD", t_sample + "AD"]:
                if not q in rec or rec[q] is None:
                    rec[q] = 0
                    if not ("feat:" + q) in has_warned:
                        logging.warn("Missing feature %s" % q)
                        has_warned["feat:" + q] = True
                else:
                    # list features
                    if q.endswith("AD") or q.endswith("QSS") or q.endswith("RPA"):
                        if type(rec[q]) is not list:
                            if not q + "_PARSE_FAIL" in has_warned:
                                logging.warn("Cannot parse %s: %s" % (q, str(rec[q])))
                                has_warned[q + "_PARSE_FAIL"] = True
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
            TLOD = float(rec["I.TLOD"])
            NLOD = float(rec["I.NLOD"])

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
                "TLOD": TLOD,
                "NLOD": NLOD,
                "N_DP": n_DP,
                "T_DP": t_DP,
                "N_DP_RATE" : n_DP_ratio,
                "T_DP_RATE" : t_DP_ratio,
                "N_GT": rec[n_sample + "GT"],
                "T_GT": rec[t_sample + "GT"],
                "N_AD": rec[n_sample + "AD"],
                "T_AD": rec[t_sample + "AD"],
                "N_QSS": rec[n_sample + "QSS"],
                "T_QSS": rec[t_sample + "QSS"],
                "N_AF": n_allele_rate,
                "T_AF": t_allele_rate,
                "ECNT": rec["I.ECNT"],
                "HCNT": rec["I.HCNT"],
                "MAX_ED": rec["I.MAX_ED"],
                "MIN_ED": rec["I.MIN_ED"],
                "I.RPA": rec["I.RPA"],
                "I.RU": rec["I.RU"],
                "tag" : tag
            }
            records.append(qrec)

        cols = ["CHROM", "POS", "REF", "ALT",
                "FILTER", "TLOD", "NLOD", "DBSNP",
                "N_DP", "T_DP", "N_DP_RATE", "T_DP_RATE", "N_GT", "T_GT",
                "N_AD", "T_AD", "N_QSS", "T_QSS",
                "N_AF", "T_AF",
                "tag"]

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df
