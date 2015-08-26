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


def roc(vcf, feature, filter_name, output_path):
    """ Calculate SNP and indel ROC. """
    tf = tempfile.NamedTemporaryFile(delete=False)

    try:
        cmdline = "bcftools query -i 'INFO/type==\"FP\"' -f '%%INFO/Q_VT\t%%INFO/%s\t%%INFO/type\t%%FILTER\\n' " \
                  "%s -o %s" % (feature, vcf.replace(" ", "\\ "), tf.name + ".fp")
        logging.info("Running: %s" % cmdline)
        subprocess.check_call(cmdline, shell=True)
        cmdline = "bcftools query -i 'INFO/type==\"TP\"' -f '%%INFO/T_VT\t%%INFO/%s\t%%INFO/type\t%%FILTER\\n' " \
                  "%s -o %s" % (feature, vcf.replace(" ", "\\ "), tf.name + ".tp")
        logging.info("Running: %s" % cmdline)
        subprocess.check_call(cmdline, shell=True)

        tf_i = open(tf.name + ".indel", "w")
        print >>tf, "type\t%s\tlabel\tfilter"
        print >>tf_i, "type\t%s\tlabel\tfilter"
        with open(tf.name + ".tp") as tp:
            for l in tp:
                if l.startswith("INDEL"):
                    tf_i.write(l)
                else:
                    tf.write(l)
        with open(tf.name + ".fp") as fp:
            for l in fp:
                if l.startswith("INDEL"):
                    tf_i.write(l)
                else:
                    tf.write(l)
        tf.close()
        tf_i.close()

        cmdline = "roc -t label -v %s -f filter --verbose " % (feature)
        if filter_name:
            cmdline += " -n %s" % filter_name

        cmdlines = cmdline + " -o %s %s" % (output_path.replace(" ", "\\ ") + ".snp.tsv", tf.name)
        logging.info("Running %s" % cmdlines)
        subprocess.check_call(cmdlines, shell=True)

        cmdlinei = cmdline + " -o %s %s" % (output_path.replace(" ", "\\ ") + ".indel.tsv", tf.name + ".indel")
        logging.info("Running %s" % cmdlinei)
        subprocess.check_call(cmdlinei, shell=True)
    finally:
        tf.close()
        _rm(tf.name)
        _rm(tf.name + ".indel")
        _rm(tf.name + ".fp")
        _rm(tf.name + ".tp")
