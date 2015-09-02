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


def _locations_tmp_bed_file(locations):
    """ turn a list of locations into a bed file """
    if type(locations) is str:
        locations = locations.split(",")
    if type(locations) is not list:
        raise Exception("Invalid list of locations (must be str or list): %s" % str(locations))

    llocations = []

    for l in locations:
        xchr, _, _pos = l.partition(":")
        start, _, end = _pos.partition("-")
        if not xchr:
            raise Exception("Invalid chromosome name in %s" % str(l))
        try:
            start = int(start)
        except:
            start = 0
        
        try:
            end = int(end)
        except:
            end = 2**31-1

        llocations.append((xchr, start, end))

    locations = sorted(llocations)

    tf = tempfile.NamedTemporaryFile(delete=False)
    for xchr, start, end in locations:
        print >>tf, "%s\t%i\t%i" % (xchr, start-1, end)
    tf.close()

    return tf.name



def run_quantify(filename, json_name=None, write_vcf=False, regions=None,
                  reference=Tools.defaultReference(), sample_column="*",
                  locations=None):
    """Run quantify and return parsed JSON

    :filename: the VCF file name
    :json_name: output file name (if None, will use a temp file)
    :write_vcf: write annotated VCF (give filename)
    :regions: dictionary of stratification region names and file names
    :reference: reference fasta path
    :sample_column: column to use (or * for all)
    :location: a location to use
    :returns: parsed counts JSON
    """

    if not json_name:
        json_name = tempfile.NamedTemporaryFile().name

    run_str = "quantify '%s:%s' -o '%s'" % (filename.replace(" ", "\\ "), sample_column, json_name)
    run_str += " -r '%s'" % reference.replace(" ", "\\ ")

    if write_vcf:
        if not write_vcf.endswith(".vcf.gz"):
            write_vcf += ".vcf.gz"
        run_str += " -v '%s'" % write_vcf

    if regions:
        for k, v in regions.iteritems():
            run_str += " -R '%s:%s'" % (k, v)

    location_file = None
    if locations:
        location_file = _locations_tmp_bed_file(locations)
        run_str += " --only '%s'" % location_file

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stdout",
                                      suffix=".log")

    logging.info("Running '%s'" % run_str)

    try:
        subprocess.check_call(run_str, shell=True, stdout=tfo, stderr=tfe)
    except:
        tfo.close()
        tfe.close()
        with open(tfo.name) as f:
            for l in f:
                logging.error("[stdout] " + l.replace("\n", ""))
        os.unlink(tfo.name)
        with open(tfe.name) as f:
            for l in f:
                logging.error("[stderr] " + l.replace("\n", ""))
        os.unlink(tfe.name)
        if location_file:
            os.unlink(location_file)
        raise

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
    if location_file:
        os.unlink(location_file)

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

    simplified_numbers = {}
    for k1, v in counts.iteritems():
        if k1.startswith("all") or v is None:
            continue
        # k1 is made in quantify.cpp
        # type + ":" + kind + ":" + tag_string + ":" + sample
        try:
            vtype, kind, tags, tq = k1.split(":", 4)
        except:
            logging.warn("Invalid k1 : %s" % k1)
            continue

        if type(snames) is list and tq not in snames:
            logging.warn("Ignoring invalid key %s" % k1)
            continue

        for k2, v2 in v.iteritems():
            # k2 is created in quantify.cpp, observation type, then allele types separated by __
            ct, _, vt = k2.partition("__")

            if ct == "nuc":
                if vt == "s":
                    altype = "SNP"
                elif vt == "i":
                    altype = "INS"
                elif vt == "d":
                    altype = "DEL"
                else:
                    raise Exception("Invalid nucleotide type: %s" % vt)
                keys1 = ["Nucleotides", "Nucleotides." + altype]
            elif ct == "al":
                if vt == "s":
                    altype = "SNP"
                elif vt == "i":
                    altype = "INS"
                elif vt == "d":
                    altype = "DEL"
                else:
                    altype = "COMPLEX"
                keys1 = ["Alleles", "Alleles." + altype]
            elif not (vt == "nc" or ct == "homref" or vt == "r"):  # ignore non-called locations in a sample
                if vt == "s" or vt == "rs":
                    altype = "SNP"
                else:
                    altype = "INDEL"
                xkeys1 = ["Locations", "Locations." + altype]
                if vt != "s":
                    xkeys1.append("Locations.detailed." + vt)

                keys1 = copy.copy(xkeys1)
                for k in xkeys1:
                    if k != "Locations":
                        keys1.append(k + "." + ct)
            elif vt == "nc" or ct == "homref" or vt == "r":
                if vt == "nc":
                    keys1 = ["Records.nocall"]
                else:
                    keys1 = ["Records.homref"]
            else:
                keys1 = ["Records.unknown"]

            keys2 = [tq + ".TOTAL", tq + "." + vtype]

            for key1 in keys1:
                if key1 not in simplified_numbers:
                    simplified_numbers[key1] = copy.copy(counter_dict)
                for key2 in keys2:
                    try:
                        simplified_numbers[key1][key2] += v2
                    except:
                        simplified_numbers[key1][key2] = v2

    return simplified_numbers
