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
                if "type" in rec:
                    roc = "Locations." + rec["type"]
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
