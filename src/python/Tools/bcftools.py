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

import os
import subprocess
import logging
import pandas
import tempfile
import gzip

import Tools


def runBcftools(*args):
    """ Run bcftools, return output
    """
    runme = "bcftools %s" % " ".join([a.replace(" ", "\\ ") for a in args])
    logging.info(runme)
    po = subprocess.Popen(runme,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

    o, e = po.communicate()

    po.wait()

    rc = po.returncode

    if rc != 0:
        raise Exception("Error running BCFTOOLS; please check if your file has issues using vcfcheck"
                        ". Return code was %i, output: %s / %s \n" % (rc, o, e))

    return o


def parseStats(output, colname="count"):
    """ Parse BCFTOOLS Stats Output """

    result = {}
    for x in output.split("\n"):
        if x.startswith("SN"):
            vx = x.split("\t")
            name = vx[2].replace("number of ", "").replace(":", "")
            count = int(vx[3])
            result[name] = count

    result = pandas.DataFrame(list(result.iteritems()), columns=["type", colname])
    return result


def countVCFRows(filename):
    """ Count the number of rows in a VCF
    :param filename: VCF file name
    :return: number of rows
    """
    if filename.endswith(".gz"):
        f = gzip.open(filename, "r")
    else:
        f = open(filename, "r")

    count = 0
    for s in f:
        if not s.startswith("#"):
            count += 1
    return count


def makeIndex(inputfile, output_name=None):
    """ Create bgzipped copy and make sure we have an up-to-date index for the given VCF file.

    Will re-write the given VCF and return a new filename (or output_name if given)

    :inputfile: input file name
    :output_name: output file name
    :returns: name of file that was written,
              inputfile if no file was written, and output_name if not None
    """
    if not output_name:
        of = tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False)
        of.close()
        output_name = of.name

    runBcftools("view",
                "-o", output_name,
                "-O", "z",
                inputfile)
    runBcftools("index", "-t", output_name)

    return output_name


def concatenateParts(output, *args):
    """ Concatenate BCF files


    Trickier than it sounds because when there are many files we might run into
    various limits like the number of open files, or the length of a command line.

    This function will bcftools concat in a tree-like fashion to avoid this.
    """
    to_delete = []
    try:
        if output.endswith(".bcf"):
            outputformat = "b"
            outputext = ".bcf"
        else:
            outputformat = "z"
            outputext = ".vcf.gz"
        if len(args) < 10:
            for x in args:
                if not os.path.exists(x + ".tbi") and not os.path.exists(x + ".csi"):
                    to_delete.append(x + ".csi")
                    runBcftools("index", "-f", x)
            cmdlist = ["concat", "-a", "-O", outputformat, "-o", output] + list(args)
            runBcftools(*cmdlist)
        else:
            # block in chunks (TODO: make parallel)
            tf1 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
            tf2 = tempfile.NamedTemporaryFile(suffix=outputext, delete=False)
            to_delete.append(tf1.name)
            to_delete.append(tf2.name)
            to_delete.append(tf1.name + ".csi")
            to_delete.append(tf2.name + ".csi")
            half1 = [tf1.name] + list(args[:len(args)/2])
            half2 = [tf2.name] + list(args[len(args)/2:])
            concatenateParts(*half1)
            runBcftools("index", tf1.name)
            concatenateParts(*half2)
            runBcftools("index", tf2.name)
            concatenateParts(output, tf1.name, tf2.name)
    finally:
        for f in to_delete:
            try:
                os.unlink(f)
            except:
                pass


# noinspection PyShadowingBuiltins
def preprocessVCF(input, output, location="",
                  pass_only=True,
                  chrprefix=True, norm=False,
                  regions=None, targets=None,
                  reference=Tools.defaultReference(),
                  filters_only=None):
    """ Preprocess a VCF + create index

    :param input: the input VCF / BCF / ...
    :param output: the output VCF
    :param location: optional location string -- comma separated
    :param pass_only: only return passing variants
    :param chrprefix: fix chromosome prefix
    :param norm: run through bcftools norm to leftshift indels
    :param regions: specify a subset of regions (traversed using tabix index, which must exist)
    :param targets: specify a subset of target regions (streaming traversal)
    :param reference: reference fasta file to use
    :param filters_only: require a set of filters (overridden by pass_only)
    """
    vargs = ["view", input]

    if pass_only:
        vargs += ["-f", "PASS,."]
    elif filters_only:
        vargs += ["-f", filters_only]

    if chrprefix:
        vargs += ["|", "perl", "-pe", "'s/^([0-9XYM])/chr$1/'", "|", "perl", "-pe", "'s/chrMT/chrM/'", "|", "bcftools", "view"]

    if targets:
        vargs += ["-T", targets, "|", "bcftools", "view"]

    if location:
        vargs += ["-t", location, "|", "bcftools", "view"]

    if output.endswith(".vcf.gz"):
        int_suffix = ".vcf.gz"
    else:
        int_suffix = ".bcf"
    tff = tempfile.NamedTemporaryFile(delete=False, suffix=int_suffix)
    try:
        # anything needs tabix? if so do an intermediate stage where we
        # index first
        if regions:
            if int_suffix == ".vcf.gz":
                vargs += ["-o", tff.name, "-O", "z"]
                runBcftools(*vargs)
                runBcftools("index", "-t", tff.name)
            else:
                vargs += ["-o", tff.name, "-O", "b"]
                runBcftools(*vargs)
                runBcftools("index", tff.name)
            vargs = ["view", tff.name, "-R", regions]

        if norm:
            vargs += ["|", "bcftools", "norm", "-f",  reference, "-c", "x", "-D"]

        vargs += ["-o", output]
        if int_suffix == ".vcf.gz":
            vargs += ["-O", "z"]
            istabix = True
        else:
            vargs += ["-O", "b"]
            istabix = False
        runBcftools(*vargs)
        if istabix:
            runBcftools("index", "-t", output)
        else:
            runBcftools("index", output)
    finally:
        try:
            os.unlink(tff.name)
        except:
            pass
        try:
            os.unlink(tff.name + ".tbi")
        except:
            pass
        try:
            os.unlink(tff.name + ".csi")
        except:
            pass


def bedOverlapCheck(filename):
    """ Check for overlaps / out of order in a bed file """
    if filename.endswith(".gz"):
        f = gzip.open(filename, "r")
    else:
        f = open(filename, "r")
    last = -1
    lines = 1
    thischr = None
    for line in f:
        l = line.split("\t")
        if len(l) < 3:
            continue
        if thischr is not None and thischr != l[0]:
            last = -1
        thischr = l[0]
        if (last-1) > int(l[1]):
            logging.warn("%s has overlapping regions at %s:%i (line %i)" % (filename, l[0], int(l[1]), lines))
            return 1
        last = int(l[2])
        lines += 1
    return 0
