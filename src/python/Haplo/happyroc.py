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
import tempfile
import subprocess
import logging


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


def roc(roc_table, feature, filter_name, output_path, rreversed):
    """ Calculate SNP and indel ROC.

    Return a dictionary of variant types and corresponding files.

    """
    # tf = tempfile.NamedTemporaryFile(delete=False)
    # tf.close()

    # files = {}
    result = {}
    # try:
    #     def getfile(vtype, ltype):
    #         fname = output_path.replace(" ", "\\ ") + \
    #             "." + vtype.lower() + "." + ltype.lower() + ".data"
    #         if fname not in files:
    #             files[fname] = {
    #                 "vtype": vtype,
    #                 "ltype": ltype,
    #                 "file": open(fname, "w")
    #             }
    #             print >>files[fname]["file"], "CHROM\tPOS\ttype\tltype\tlabel\t%s\tfilter" % feature
    #         return files[fname]["file"]

    #     def output_line(l):
    #         ll = l.split("\t")
    #         f1 = getfile(ll[2], "all")
    #         f2 = getfile(ll[2], ll[3])
    #         f1.write(l)
    #         f2.write(l)

    #     # split / distribute ROC table recorts
    #     with open(roc_table) as rt:
    #         # skip header
    #         next(rt)
    #         for l in rt:
    #             output_line(l)

    #     cmdline = "roc -t label -v %s -f filter --verbose -R %i" % (feature, 1 if rreversed else 0)
    #     if filter_name:
    #         cmdline += " -n '%s'" % filter_name

    #     for n, ff in files.iteritems():
    #         ff["file"].close()
    #         fname = output_path.replace(" ", "\\ ") + \
    #             "." + ff["vtype"].lower() + "." + ff["ltype"].lower() + \
    #             ".tsv"
    #         result["Locations." + ff["vtype"].upper() + ("." + ff["ltype"].lower()
    #                                                      if ff["ltype"] not in ["", "all"] else "")] = fname
    #         cmdlines = cmdline + " -o %s %s" % (fname, n)
    #         _run(cmdlines)

    # finally:
    #     tf.close()
    #     _rm(tf.name)
    #     _rm(tf.name + ".indel")
    #     _rm(tf.name + ".fp")
    #     _rm(tf.name + ".tp")
    #     _rm(tf.name + ".fn")
    #     # remove intermediate data -- might add an option to keep
    #     for n in files.iterkeys():
    #         _rm(n)
    return result
