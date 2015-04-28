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
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt

import os
import logging
import subprocess
import tempfile
import gzip
import re
import time
import json


def field(val):
    """ extract field into result, guess type """
    if "," in val:
        val = map(field, val.split(","))
    else:
        done = False
        try:
            val = int(val)
            done = True
        except:
            pass

        if done:
            return val
        try:
            val = float(val)
        except:
            pass
    return val


def getInfo(istr):
    """Split VCF INFO String"""
    spi = istr.split(";")
    res = {}
    for x in spi:
        ax = [q.strip() for q in x.split("=", 1)]

        if len(ax) == 1:
            res[ax[0]] = True
        else:
            res[ax[0]] = field(ax[1])
    return res


def getFormats(fstr, fsample):
    """Split VCF INFO String"""
    spf = fstr.split(":")
    sps = fsample.split(":")

    res = {}
    for i, x in enumerate(spf):
        res[x] = field(sps[i])
    return res


def splitIndex(ffield):
    """If a field is given as an index into a list, split these"""
    rm = re.match("(.*)\[([0-9]+)\]", ffield)

    if rm:
        return rm.group(1), int(rm.group(2))
    else:
        return ffield, None


def vcfExtract(vcfname, features, filterfun=None):
    """ Given a list of VCF features, get tab-separated list from VCF file

    :param vcfname: the vcf file name
    :param features: list of features to extract
    :param filterfun: optional filter function to ignore certain lines in the VCF

    """

    if vcfname.endswith(".gz"):
        ff = gzip.GzipFile(vcfname)
    else:
        ff = open(vcfname)

    feature_index = [splitIndex(f) for f in features]

    start = time.time()
    last_time = start
    nrecords = 0
    previous_end = 0
    for line in ff:
        if line.startswith("#"):
            continue
        line = line.replace("\n", "")
        if filterfun and filterfun(line):
            continue

        nrecords += 1
        lstart = time.time()

        if lstart - last_time > 10:
            last_time = lstart
            total = lstart-start
            # noinspection PyBroadException
            try:
                tpr = 1000000.0*total/float(nrecords)
            except:
                tpr = -1
            logging.info("Since start: %i records %.2f seconds, %.2f us/record." % (nrecords, total, tpr))

        spl = line.split("\t")
        current = []
        curinfo = None
        curformats = {}
        for i, f in enumerate(features):
            if f.lower().startswith("chr"):
                current.append(spl[0])
            elif f.lower().startswith("pos"):
                current.append(int(spl[1]))
            elif f.lower().startswith("id"):
                current.append(spl[2])
            elif f.lower().startswith("ref"):
                current.append(spl[3])
            elif f.lower().startswith("alt"):
                val = spl[4].split(",")
                if feature_index[i][1] is not None:
                    if feature_index[i][1] < len(val):
                        val = val[feature_index[i][1]]
                    else:
                        val = None
                current.append(val)
            elif f.lower().startswith("qual"):
                try:
                    current.append(float(spl[5]))
                except:
                    current.append(None)
            elif f.lower().startswith("fil"):
                if spl[6] == "PASS" or spl[6] == ".":
                    val = []
                else:
                    val = spl[6].split(",")

                if feature_index[i][1] is not None:
                    if feature_index[i][1] < len(val):
                        val = val[feature_index[i][1]]
                    else:
                        val = None

                current.append(val)
            elif f.startswith("I."):
                if curinfo is None:
                    curinfo = getInfo(spl[7])
                val = None
                try:
                    ff, ii = feature_index[i]
                    val = curinfo[ff[2:]]

                    if ii is not None:
                        if ii < len(val):
                            val = val[ii]
                        else:
                            val = None
                except:
                    pass
                current.append(val)
            elif f.startswith("S."):
                ff, ii = feature_index[i]
                dx = ff.split(".", 3)
                sample = int(dx[1])
                field = dx[2]

                val = None
                try:
                    if not sample in curformats:
                        curformats[sample] = getFormats(spl[8], spl[8+sample])
                    val = curformats[sample][field]

                    if ii is not None:
                        if ii < len(val):
                            val = val[ii]
                        else:
                            val = None
                except:
                    pass
                current.append(val)
            else:
                current.append(f)
        yield current


def extractHeaders(vcfname):
    """ Read the header lines from a VCF file """
    if vcfname.endswith(".gz"):
        ff = gzip.GzipFile(vcfname)
    else:
        ff = open(vcfname)

    for l in ff:
        if l.startswith("#"):
            yield l.replace("\n", "")
        else:
            break


def extractHeadersJSON(vcfname):
    """ Extract the VCF header and turn into JSON
    :param vcfname: VCF file name
    :return: VCF header in JSON format
    """
    tf = tempfile.NamedTemporaryFile(delete=False)
    tf.close()
    vfh = {}

    try:
        sp = subprocess.Popen("vcfhdr2json %s %s" % (vcfname, tf.name),
                              shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = sp.communicate()

        if sp.returncode != 0:
            raise Exception("Samtools call failed: %s / %s" % (o, e))

        vfh = json.load(open(tf.name))
    finally:
        try:
            os.unlink(tf.name)
        except:
            pass

    return vfh

