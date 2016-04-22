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
# 01/04/2015
#
# Calculate diploid SNP and Indel ROCs
#
# This is a little different from the somatic case that is handled in Tools/roc.py
# -- we don't want to read in any of the data since we have a lot more of it.
#

import os
import re
import subprocess
import logging
import pandas
import itertools

def _rm(f):
    """ Quietly delete """
    try:
        os.unlink(f)
    except:
        pass


def _run(cmd):
    """ Run something, log output
    """
    logging.info(cmd)
    po = subprocess.Popen(cmd,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

    o, e = po.communicate()

    po.wait()

    rc = po.returncode

    if rc != 0:
        raise Exception("Error running ROC. Return code was %i\n" % rc)

    logging.info(o)
    logging.info(e)

    return o


def roc(roc_table, output_path, filter_handling=None):
    """ Calculate SNP and indel ROC.

    Return a dictionary of variant types and corresponding files.

    :param filter_handling: can be None, "PASS" or "ALL" to filter rows based on the "Filter" column.
                            this is necessary because vcfeval doesn't preserve filter information
                            in GA4GH output mode, so we need to remve the corresponding rows here

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
                        roc = "Locations." + rec["type"]
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
                        roc = "Locations." + rec["type"] + ".PASS"
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

    for k, v in result.items():
        result[k] = _postprocessRocData(pandas.DataFrame(v, columns=header))
        vt = re.sub("[^A-Za-z0-9\\.\\-_]", "_", k, flags=re.IGNORECASE)
        if output_path:
            result[k].to_csv(output_path + "." + vt + ".csv", index=False)

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
        roctable[col] = roctable[col].convert_objects(convert_numeric=True)
        if c:
            roctable[col] = roctable[col].apply(c)

    roctable.sort(["Type", "Subtype", "Subset", "Filter", "Genotype", "QQ.Field", "QQ"], inplace=True)

    permuted_columns = ["Type",
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
                        "TRUTH.TOTAL",
                        "TRUTH.TP",
                        "TRUTH.FN",
                        "QUERY.TOTAL",
                        "QUERY.TP",
                        "QUERY.FP",
                        "QUERY.UNK",
                        "FP.gt",
                        "FP.al",
                        "TRUTH.TOTAL.TiTv_ratio",
                        "TRUTH.TOTAL.het_hom_ratio",
                        "TRUTH.FN.TiTv_ratio",
                        "TRUTH.FN.het_hom_ratio",
                        "TRUTH.TP.TiTv_ratio",
                        "TRUTH.TP.het_hom_ratio",
                        "QUERY.FP.TiTv_ratio",
                        "QUERY.FP.het_hom_ratio",
                        "QUERY.TP.TiTv_ratio",
                        "QUERY.TOTAL.TiTv_ratio",
                        "QUERY.TOTAL.het_hom_ratio",
                        "QUERY.TP.het_hom_ratio",
                        "QUERY.UNK.TiTv_ratio",
                        "QUERY.UNK.het_hom_ratio"]

    return roctable[permuted_columns]
