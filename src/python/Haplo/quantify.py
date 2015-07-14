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
# Process raw counts coming out of quantify

import os
import tempfile
import subprocess
import copy
import json
import logging

import Tools


def run_quantify(filename, json_name=None, write_vcf=False, regions=None,
                 reference=Tools.defaultReference(), sample_column="*"):
    """Run quantify and return parsed JSON

    :filename: the VCF file name
    :json_name: output file name (if None, will use a temp file)
    :write_vcf: write annotated VCF (give filename)
    :regions: dictionary of region names and file names
    :reference: reference fasta path
    :sample_column: column to use (or * for all)
    :returns: parsed counts JSON
    """

    if not json_name:
        json_name = tempfile.NamedTemporaryFile().name

    run_str = "quantify '%s:%s' -o '%s'" % (filename, sample_column, json_name)
    run_str += " -r '%s'" % reference

    if write_vcf:
        if not write_vcf.endswith(".vcf.gz"):
            write_vcf += ".vcf.gz"
        run_str += " -v '%s'" % write_vcf

    if regions:
        for k, v in regions.iteritems():
            run_str += " -R '%s:%s'" % (k, v)

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stdout",
                                      suffix=".log")

    logging.info("Running '%s'" % run_str)
    subprocess.check_call(run_str, shell=True, stdout=tfo, stderr=tfe)
    tfo.close()
    tfe.close()
    with open(tfo.name) as f:
        for l in f:
            logging.info("[stdout] " + l.replace("\n", ""))
    os.unlink(tfo.name)
    with open(tfe.name) as f:
        for l in f:
            logging.info("[stderr] " + l.replace("\n", ""))
    os.unlink(tfe.name)

    if write_vcf:
        to_run = "tabix -p vcf '%s'" % write_vcf
        logging.info("Running '%s'" % to_run)
        subprocess.check_call(to_run, shell=True)
    return json.load(open(json_name))


def simplify_counts(counts, snames=None):
    """ Return simplified counts from quantify  """

    if not snames:
        snames = ["TRUTH", "QUERY"]

    counter_dict = {}
    if type(snames) is list:
        for sn in snames:
            counter_dict[sn + ".TOTAL"] = 0
            counter_dict[sn + ".TP"] = 0
            counter_dict[sn + ".FP"] = 0
            counter_dict[sn + ".FN"] = 0
            counter_dict[sn + ".UNK"] = 0
            counter_dict[sn + ".TP.HC"] = 0
            counter_dict[sn + ".FP.HC"] = 0
            counter_dict[sn + ".FN.HC"] = 0
            counter_dict[sn + ".UNK.HC"] = 0

    simplified_numbers = {
        "Alleles" : copy.copy(counter_dict),
        "Alleles.SNP" : copy.copy(counter_dict),
        "Alleles.INS" : copy.copy(counter_dict),
        "Alleles.DEL" : copy.copy(counter_dict),
        "Alleles.other" : copy.copy(counter_dict),
        "Locations" : copy.copy(counter_dict),
        "Locations.het" : copy.copy(counter_dict),
        "Locations.homalt" : copy.copy(counter_dict),
        "Locations.hetalt" : copy.copy(counter_dict),
        "Locations.unknown" : copy.copy(counter_dict),
        "Locations.SNP" : copy.copy(counter_dict),
        "Locations.SNP.het" : copy.copy(counter_dict),
        "Locations.SNP.homalt" : copy.copy(counter_dict),
        "Locations.SNP.hetalt" : copy.copy(counter_dict),
        "Locations.SNP.unknown" : copy.copy(counter_dict),
        "Locations.INDEL" : copy.copy(counter_dict),
        "Locations.INDEL.het" : copy.copy(counter_dict),
        "Locations.INDEL.homalt" : copy.copy(counter_dict),
        "Locations.INDEL.hetalt" : copy.copy(counter_dict),
        "Locations.INDEL.unknown" : copy.copy(counter_dict),
    }

    for k1, v in counts.iteritems():
        if k1.startswith("all") or v is None:
            continue
        # k1 is made in quantify.cpp
        # type + ":" + hapmatch + ":" + kind + ":" + tag_string
        vtype, hapmatch, kind, tags, tq = k1.split(":", 5)

        if type(snames) is list and not tq in snames:
            logging.warn("Ignoring invalid key %s" % k1)
            continue

        if not vtype in ["TP", "FP", "FN"]:
            vtype = "UNK"

        # these don't make sense to count
        if vtype == "FN" and tq == "QUERY":
            continue

        # turn FP outside of confident regions into UNK
        if vtype == "FP" and not "CONF" in tags and kind == "missing":
            vtype = "UNK"

        for k2, v2 in v.iteritems():
            if "__" in k2:
                atype, _, gt = k2.partition("__")
            else:
                atype = k2
                gt = None

            if atype == "alleles":
                atype = "Alleles"
            elif atype == "locations":
                atype = "Locations"
            elif atype == "snp":
                if gt:
                    atype = "Locations.SNP"
                else:
                    atype = "Alleles.SNP"
            elif atype in ["ins", "del"]: # for alleles, split out ins / del
                if gt:  # for locations, everything not a SNP is an Indel
                    atype = "Locations.INDEL"
                else:
                    atype = "Alleles." + atype.upper()
            elif atype in ["het", "hetalt", "homalt", "unknown_gt"]: # for alleles, split out ins / del
                gt = atype
                atype = "Locations"
            else:
                if gt:
                    atype = "Locations.INDEL"
                else:
                    atype = "Alleles.other"

            if gt and gt not in ["het", "homalt", "hetalt"]:
                gt = "unknown"

            key1 = atype
            if gt:
                key1 += "." + gt

            keys1 = [key1]
            # add to SNP / indel locations
            for qprefix in ["Locations.SNP", "Locations.INDEL"]:
                if key1.startswith(qprefix):
                    keys1.append(qprefix)

            keys2 = [tq + ".TOTAL", tq + "." + vtype]

            # haplotype-matched FP/FN become TP
            if hapmatch:
                keys2.append(tq + ".TP.HC")
            else:
                keys2.append(tq + "." + vtype + ".HC")

            for key1 in keys1:
                for key2 in keys2:
                    try:
                        simplified_numbers[key1][key2] += v2
                    except:
                        simplified_numbers[key1][key2] = v2

    return simplified_numbers
