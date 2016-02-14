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
            end = 2 ** 31 - 1

        llocations.append((xchr, start, end))

    locations = sorted(llocations)

    tf = tempfile.NamedTemporaryFile(delete=False)
    for xchr, start, end in locations:
        print >> tf, "%s\t%i\t%i" % (xchr, start - 1, end)
    tf.close()

    return tf.name


def run_quantify(filename, json_name=None, write_vcf=False, regions=None,
                 reference=Tools.defaultReference(),
                 locations=None, threads=1):
    """Run quantify and return parsed JSON

    :param filename: the VCF file name
    :param json_name: output file name (if None, will use a temp file)
    :param write_vcf: write annotated VCF (give filename)
    :type write_vcf: str
    :param regions: dictionary of stratification region names and file names
    :param reference: reference fasta path
    :param locations: a location to use
    :returns: parsed counts JSON
    """

    if not json_name:
        json_name = tempfile.NamedTemporaryFile().name

    run_str = "quantify '%s' -o '%s'" % (filename.replace(" ", "\\ "), json_name)
    run_str += " -r '%s'" % reference.replace(" ", "\\ ")
    run_str += " --threads %i" % threads

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
    """ Return simplified counts from quantify
    :param counts: counts from running quantify
    :param snames: sample names
    """

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
            elif not (vt == "nc" or ct == "homref" or ct == "fail" or vt == "r"):  # ignore non-called locations in a sample
                # process locations
                if ct in ["Transitions", "Transversions"]:
                    # these are additional counts. Every SNP we see in
                    # here is also seen separately below.
                    keys1 = ["Locations.SNP." + ct]
                else:
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
            elif ct == "fail":
                keys1 = ["Records.fail"]
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

    for vt in ["Locations.SNP", "Locations.INDEL"]:
        if vt not in simplified_numbers:
            continue
        for sample in snames:
            for ct in ["TP", "FP", "FN", "UNK", "TOTAL"]:
                homalt_count = 0
                try:
                    homalt_count = simplified_numbers[vt + ".homalt"][sample + "." + ct]
                except:
                    pass
                het_count = 0
                try:
                    het_count = simplified_numbers[vt + ".het"][sample + "." + ct]
                except:
                    pass
                hetalt_count = 0
                try:
                    hetalt_count = simplified_numbers[vt + ".hetalt"][sample + "." + ct]
                except:
                    pass

                if homalt_count != 0:
                    simplified_numbers[vt][sample + "." + ct + ".het_hom_ratio"] = float(het_count) / float(homalt_count)
                    simplified_numbers[vt][sample + "." + ct + ".hethetalt_hom_ratio"] = float(het_count + hetalt_count) / float(homalt_count)

                if vt == "Locations.SNP":
                    ti_count = 0
                    try:
                        ti_count = simplified_numbers[vt + ".Transitions"][sample + "." + ct]
                    except:
                        pass

                    tv_count = 0
                    try:
                        tv_count = simplified_numbers[vt + ".Transversions"][sample + "." + ct]
                    except:
                        pass

                    if tv_count > 0:
                        simplified_numbers[vt][sample + "." + ct + ".TiTv_ratio"] = float(ti_count) / float(tv_count)

    return simplified_numbers
