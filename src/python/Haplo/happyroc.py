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


def roc(roc_table, output_path):
    """ Calculate SNP and indel ROC.

    Return a dictionary of variant types and corresponding files.

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
                try:
                    if rec["type"] in ["SNP", "INDEL"] \
                       and rec["filter"] == "ALL" \
                       and rec["subset"] == "*" \
                       and rec["genotype"] == "*" \
                       and rec["subtype"] == "*" \
                       and rec[header[5]] != "*":  # this is the ROC score field
                        roc = "Locations." + rec["type"]
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
        result[k] = pandas.DataFrame(v, columns=header)
        vt = k.replace("f:", "")
        vt = re.sub("[^A-Za-z0-9\\.\\-_]", "_", vt, flags=re.IGNORECASE)
        if output_path:
            result[k].to_csv(output_path + "." + vt + ".csv")

    return result


def postprocessRocData(roctable, qq):
    """ post-process ROC data by adding on Ti/Tv ratios, het-hom ratios,
    """

    for vt in ["TRUTH.TP", "TRUTH.FN", "QUERY.FP", "QUERY.UNK", "TRUTH.TOTAL", "QUERY.TOTAL"]:
        for var in ["TiTv_ratio", "het_hom_ratio"]:
            roctable[vt + "." + var] = "."

    unique_types = roctable["Type"].unique()
    # unique_subsets = roctable["Subset"].unique()
    unique_subsets = "*"
    unique_filters = roctable["Filter"].unique()

    # try:
    #     unique_quals = roctable[qq].unique()
    # except:
    unique_quals = "*"

    for t in unique_types:
        for s in unique_subsets:
            for f in unique_filters:
                for q in unique_quals:
                    logging.info(str((t, s, f, q)))
                    for vt in ["TRUTH.TP", "TRUTH.FN", "QUERY.FP", "QUERY.UNK", "TRUTH.TOTAL", "QUERY.TOTAL"]:
                        selection = (roctable["Type"] == t) & (roctable["Subset"] == s) & (roctable["Filter"] == f) & (roctable[qq] == q)
                        unique_gts = roctable[selection]["Genotype"].unique()
                        unique_sts = roctable[selection]["Subtype"].unique()

                        for st in unique_sts:
                            try:
                                hets = roctable.ix[selection & (roctable["Genotype"] == "het") & (roctable["Subtype"] == st), vt]
                                homs = roctable.ix[selection & (roctable["Genotype"] == "homalt") & (roctable["Subtype"] == st), vt]
                                assert hets.size == 1
                                assert homs.size == 1
                                hets = float(hets.values[0])
                                homs = float(homs.values[0])
                                roctable.ix[selection & (roctable["Genotype"] == "*") & (roctable["Subtype"] == st), vt + ".het_hom_ratio"] = hets / homs
                            except:
                                pass

                        if t == "SNP":
                            for gt in unique_gts:
                                try:
                                    tis = roctable.ix[selection & (roctable["Genotype"] == gt) & (roctable["Subtype"] == "ti"), vt]
                                    tvs = roctable.ix[selection & (roctable["Genotype"] == gt) & (roctable["Subtype"] == "tv"), vt]
                                    assert tis.size == 1
                                    assert tvs.size == 1
                                    tis = float(tis.values[0])
                                    tvs = float(tvs.values[0])
                                    roctable.ix[selection & (roctable["Genotype"] == gt) & (roctable["Subtype"] == "*"), vt + ".TiTv_ratio"] = tis / tvs
                                except:
                                    pass

    return roctable
