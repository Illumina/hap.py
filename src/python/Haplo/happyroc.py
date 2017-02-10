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
# 01/04/2015
#
# Calculate diploid SNP and Indel ROCs
#
# This is a little different from the somatic case that is handled in Tools/roc.py
# -- we don't want to read in any of the data since we have a lot more of it.
#

import os
import re
import logging
import pandas
import numpy as np
import itertools

from Tools import ci

RESULT_ALLCOLUMNS = ["Type",
                     "Subtype",
                     "Subset",
                     "Filter",
                     "Genotype",
                     "QQ.Field",
                     "QQ",
                     "METRIC.Recall",
                     "METRIC.Precision",
                     "METRIC.Frac_NA",
                     "METRIC.F1_Score",
                     "FP.gt",
                     "FP.al",
                     "Subset.Size",
                     "Subset.IS_CONF.Size",
                     "Subset.Level"]

RESULT_ALLDTYPES = [str,
                    str,
                    str,
                    str,
                    str,
                    str,
                    str,
                    float,
                    float,
                    float,
                    float,
                    float,
                    float,
                    float,
                    float,
                    int]

for count_type in ["TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
                   "QUERY.TOTAL", "QUERY.TP", "QUERY.FP",
                   "QUERY.UNK"]:
    RESULT_ALLCOLUMNS.append(count_type)
    RESULT_ALLCOLUMNS.append(count_type + ".ti")
    RESULT_ALLCOLUMNS.append(count_type + ".tv")
    RESULT_ALLCOLUMNS.append(count_type + ".het")
    RESULT_ALLCOLUMNS.append(count_type + ".homalt")
    RESULT_ALLCOLUMNS.append(count_type + ".TiTv_ratio")
    RESULT_ALLCOLUMNS.append(count_type + ".het_hom_ratio")
    RESULT_ALLDTYPES += [float] * 7


def roc(roc_table, output_path,
        filter_handling=None,
        ci_alpha=0.05,
        total_region_size=None):
    """ Calculate SNP and indel ROC.

    Return a dictionary of variant types and corresponding files.

    :param filter_handling: can be None, "PASS" or "ALL" to filter rows based on the "Filter" column.
                            this is necessary because vcfeval doesn't preserve filter information
                            in GA4GH output mode, so we need to remove the corresponding rows here
    :param ci_alpha: Jeffrey's CI confidence level for recall, precision, na
    :param total_region_size: correct Subset.Size for "*" region if a subset was selected in hap.py

    """
    result = {}
    header = None
    with open(roc_table) as rt:
        for l in rt:
            l = l.strip()
            if not header:
                header = l.split("\t")
            else:
                rec = {}
                for k, v in itertools.izip(header, l.split("\t")):
                    rec[k] = v

                if filter_handling:
                    try:
                        if rec["Filter"] != filter_handling:
                            continue
                    except:
                        pass

                try:
                    if rec["Type"] in ["SNP", "INDEL"] \
                       and rec["Filter"] == "ALL" \
                       and rec["Subset"] == "*" \
                       and rec["Genotype"] == "*" \
                       and rec["Subtype"] == "*" \
                       and rec["QQ"] != "*":  # this is the ROC score field
                        roc = "Locations." + rec["Type"]
                        if roc not in result:
                            result[roc] = [rec]
                        else:
                            result[roc].append(rec)
                except:
                    pass

                try:
                    if rec["Type"] in ["SNP", "INDEL"] \
                       and rec["Filter"] == "PASS" \
                       and rec["Subset"] == "*" \
                       and rec["Genotype"] == "*" \
                       and rec["Subtype"] == "*" \
                       and rec["QQ"] != "*":  # this is the ROC score field
                        roc = "Locations." + rec["Type"] + ".PASS"
                        if roc not in result:
                            result[roc] = [rec]
                        else:
                            result[roc].append(rec)
                except:
                    pass

                try:
                    if rec["Type"] in ["SNP", "INDEL"] \
                       and rec["Filter"] == "SEL" \
                       and rec["Subset"] == "*" \
                       and rec["Genotype"] == "*" \
                       and rec["Subtype"] == "*" \
                       and rec["QQ"] != "*":  # this is the ROC score field
                        roc = "Locations." + rec["Type"] + ".SEL"
                        if roc not in result:
                            result[roc] = [rec]
                        else:
                            result[roc].append(rec)
                except:
                    pass

                roc = "all"
                if roc not in result:
                    result[roc] = [rec]
                else:
                    result[roc].append(rec)

    if "all" not in result:
        # minimal empty DF
        minidata = [{"Type": "SNP", "Subtype": "*", "Filter": "ALL", "Genotype": "*", "Subset": "*", "QQ": "*"} for _ in xrange(2)]
        minidata[1]["Type"] = "INDEL"
        result["all"] = pandas.DataFrame(minidata, columns=RESULT_ALLCOLUMNS)
        for i, c in enumerate(RESULT_ALLCOLUMNS):
            result["all"][c] = result["all"][c].astype(RESULT_ALLDTYPES[i], raise_on_error=False)

    for k, v in result.items():
        result[k] = _postprocessRocData(pandas.DataFrame(v, columns=RESULT_ALLCOLUMNS))

        # compute ratios
        for count_type in ["TRUTH.TOTAL", "TRUTH.FN", "TRUTH.TP", "QUERY.FP",
                           "QUERY.TP", "QUERY.TOTAL", "QUERY.UNK"]:
            result[k][count_type + ".TiTv_ratio"] = pandas.to_numeric(result[k][count_type + ".ti"], errors="coerce") / pandas.to_numeric(result[k][count_type + ".tv"], errors="coerce")
            result[k][count_type + ".het_hom_ratio"] = pandas.to_numeric(result[k][count_type + ".het"], errors="coerce") / pandas.to_numeric(result[k][count_type + ".homalt"], errors="coerce")
            result[k][count_type + ".TiTv_ratio"].replace([np.inf, -np.inf], np.nan, inplace=True)
            result[k][count_type + ".het_hom_ratio"].replace([np.inf, -np.inf], np.nan, inplace=True)

        if 0 < ci_alpha < 1:
            logging.info("Computing recall CIs for %s" % k)
            rc, rc_min, rc_max = ci.binomialCI( result[k]["TRUTH.TP"].values,
                                               (result[k]["TRUTH.TP"] + result[k]["TRUTH.FN"]).values,
                                               ci_alpha)
            result[k]["METRIC.Recall.Lower"] = rc_min
            result[k]["METRIC.Recall.Upper"] = rc_max

            logging.info("Computing precision CIs for %s" % k)
            pc, pc_min, pc_max = ci.binomialCI( result[k]["QUERY.TP"].values,
                                               (result[k]["QUERY.TP"] + result[k]["QUERY.FP"]).values,
                                               ci_alpha)
            result[k]["METRIC.Precision.Lower"] = pc_min
            result[k]["METRIC.Precision.Upper"] = pc_max

            logging.info("Computing Frac_NA CIs for %s" % k)
            fna, fna_min, fna_max = ci.binomialCI(result[k]["QUERY.UNK"].values,
                                                  result[k]["QUERY.TOTAL"].values,
                                                  ci_alpha)
            result[k]["METRIC.Frac_NA.Lower"] = fna_min
            result[k]["METRIC.Frac_NA.Upper"] = fna_max

        # write correct subset.size
        if total_region_size is not None:
            result[k].loc[result[k]["Subset"] == "*", "Subset.Size"] = total_region_size

        vt = re.sub("[^A-Za-z0-9\\.\\-_]", "_", k, flags=re.IGNORECASE)
        if output_path:
            result[k].to_csv(output_path + "." + vt + ".csv.gz", index=False,
                             compression="gzip")

    return result


def _postprocessRocData(roctable):
    """ post-process ROC data by correcting the types,
        sorting and changing the ordering of the columns
    """
    if roctable.empty:
        return roctable

    def safe_int(f):
        try:
            return int(f)
        except:
            return 0

    typeconv = [("METRIC.Recall", None),
                ("METRIC.Precision", None),
                ("METRIC.Frac_NA", None),
                ("TRUTH.TP", safe_int),
                ("TRUTH.FN", safe_int),
                ("QUERY.TP", safe_int),
                ("QUERY.FP", safe_int),
                ("QUERY.UNK", safe_int),
                ("QUERY.TOTAL", safe_int),
                ("TRUTH.TOTAL", safe_int),
                ("FP.al", safe_int),
                ("FP.gt", safe_int),
                ("TRUTH.TOTAL.TiTv_ratio", None),
                ("TRUTH.TOTAL.het_hom_ratio", None),
                ("TRUTH.FN.TiTv_ratio", None),
                ("TRUTH.FN.het_hom_ratio", None),
                ("TRUTH.TP.TiTv_ratio", None),
                ("TRUTH.TP.het_hom_ratio", None),
                ("METRIC.F1_Score", None),
                ("QUERY.FP.TiTv_ratio", None),
                ("QUERY.FP.het_hom_ratio", None),
                ("QUERY.TP.TiTv_ratio", None),
                ("QUERY.TOTAL.TiTv_ratio", None),
                ("QUERY.TOTAL.het_hom_ratio", None),
                ("QUERY.TP.het_hom_ratio", None),
                ("QUERY.UNK.TiTv_ratio", None),
                ("QUERY.UNK.het_hom_ratio", None)]

    for col, c in typeconv:
        roctable[col] = pandas.to_numeric(roctable[col], errors="coerce")
        if c:
            roctable[col] = roctable[col].apply(c)

    roctable.sort_values(["Type", "Subtype", "Subset", "Filter", "Genotype", "QQ.Field", "QQ"], inplace=True)

    return roctable[RESULT_ALLCOLUMNS]
