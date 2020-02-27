#!/usr/bin/env python
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
# 9/9/2014
#
# Preprocessing for a VCF file
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import sys
import os
import argparse
import logging
import subprocess
import multiprocessing
import tempfile
import time
import pipes

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools import vcfextract
from Tools.bcftools import preprocessVCF
from Tools.fastasize import fastaContigLengths
from Tools.bcftools import runBcftools
import Haplo.partialcredit


def hasChrPrefix(chrlist):
    """ returns if list of chr names has a chr prefix or not """

    noprefix = map(str, range(23)) + ["X", "Y", "MT"]
    withprefix = ["chr" + x for x in map(str, range(23)) + ["X", "Y", "M"]]

    count_noprefix = len(list(set(noprefix) & set(chrlist)))
    count_prefix = len(list(set(withprefix) & set(chrlist)))

    # None == undecided
    if count_prefix == count_noprefix:
        return None

    return count_noprefix < count_prefix


def preprocess(vcf_input,
               vcf_output,
               reference,
               locations=None,
               filters=None,
               fixchr=None,
               regions=None,
               targets=None,
               leftshift=True,
               decompose=True,
               bcftools_norm=False,
               windowsize=10000,
               threads=1,
               gender=None,
               somatic_allele_conversion=False,
               sample="SAMPLE",
               filter_nonref=True,
               convert_gvcf_to_vcf=False):
    """ Preprocess a single VCF file

    :param vcf_input: input file name
    :param vcf_output: output file name
    :param reference: reference fasta name
    :param locations: list of locations or None
    :param filters: list of filters to apply ("*" to only allow PASS)
    :param fixchr: None for auto, or True/False -- fix chr prefix to match reference
    :param regions: regions bed file
    :param targets: targets bed file
    :param leftshift: left-shift variants
    :param decompose: decompose variants
    :param bcftools_norm: use bcftools_norm
    :param windowsize: normalisation window size
    :param threads: number of threads to for preprcessing
    :param gender: the sex of the sample ("male" / "female" / "auto" / None)
    :param somatic_allele_conversion: convert somatic alleles -- False / half / het / hemi / hom
    :param sample: when using somatic_allele_conversion, name of the output sample
    :param filter_nonref: remove any variants genotyped as <NON_REF>

    :return: the sex if auto-determined (otherwise the same value as sex parameter)
    """

    tempfiles = []
    try:
        # If the input is in BCF format, we can continue to
        # process it in bcf
        # if it is in .vcf.gz, don't try to convert it to
        # bcf because there are a range of things that can
        # go wrong there (e.g. undefined contigs and bcftools
        # segfaults)
        if vcf_input.endswith(".bcf") or vcf_output.endswith(".bcf"):
            int_suffix = ".bcf"
            int_format = "b"
            if not vcf_input.endswith(".bcf") and vcf_output.endswith(".bcf"):
                logging.warn("Turning vcf into bcf can cause problems when headers are not consistent with all "
                             "records in the file. I will run vcfcheck to see if we will run into trouble. "
                             "To save time in the future, consider converting your files into bcf using bcftools before"
                             " running pre.py.")
        else:
            int_suffix = ".vcf.gz"
            int_format = "z"

        # HAP-317 always check for BCF errors since preprocessing tools now require valid headers
        mf = subprocess.check_output("vcfcheck %s --check-bcf-errors 1" % pipes.quote(vcf_input), shell=True)

        if gender == "auto":
            logging.info(mf)
            if "female" in mf:
                gender = "female"
            else:
                gender = "male"

        h = vcfextract.extractHeadersJSON(vcf_input)
        reference_contigs = set(fastaContigLengths(reference).keys())
        reference_has_chr_prefix = hasChrPrefix(reference_contigs)

        allfilters = []
        for f in h["fields"]:
            try:
                if f["key"] == "FILTER":
                    allfilters.append(f["values"]["ID"])
            except:
                logging.warn("ignoring header: %s" % str(f))

        required_filters = None
        if filters:
            fts = filters.split(",")
            required_filters = ",".join(list(set(["PASS", "."] + [x for x in allfilters if x not in fts])))

        if fixchr is None:
            try:
                if not h["tabix"]:
                    logging.warn("input file is not tabix indexed, consider doing this in advance for performance reasons")
                    vtf = tempfile.NamedTemporaryFile(delete=False,
                                                      suffix=int_suffix)
                    vtf.close()
                    tempfiles.append(vtf.name)
                    runBcftools("view", "-o", vtf.name, "-O", int_format, vcf_input)
                    runBcftools("index", vtf.name)
                    h2 = vcfextract.extractHeadersJSON(vcf_input)
                    chrlist = h2["tabix"]["chromosomes"]
                else:
                    chrlist = h["tabix"]["chromosomes"]
                vcf_has_chr_prefix = hasChrPrefix(chrlist)

                if reference_has_chr_prefix and not vcf_has_chr_prefix:
                    fixchr = True
            except:
                logging.warn("Guessing the chr prefix in %s has failed." % vcf_input)

        if leftshift or decompose: # all these require preprocessing
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              suffix=int_suffix)
            vtf.close()
            tempfiles.append(vtf.name)
            vtf = vtf.name
        else:
            vtf = vcf_output

        preprocessVCF(vcf_input,
                      vtf,
                      locations,
                      filters == "*",
                      fixchr,
                      bcftools_norm,
                      regions,
                      targets,
                      reference,
                      required_filters,
                      somatic_allele_conversion=somatic_allele_conversion,
                      sample=sample,
                      filter_nonref=filter_nonref,
                      convert_gvcf=convert_gvcf_to_vcf,
                      num_threads=threads)

        if leftshift or decompose or gender == "male":
            Haplo.partialcredit.partialCredit(vtf,
                                              vcf_output,
                                              reference,
                                              locations,
                                              threads=threads,
                                              window=windowsize,
                                              leftshift=leftshift,
                                              decompose=decompose,
                                              haploid_x=gender == "male")
    finally:
        for t in tempfiles:
            try:
                os.unlink(t)
            except:
                pass

    return gender

def preprocessWrapper(args):
    """ wrapper for running in parallel """

    starttime = time.time()
    logging.info("Preprocessing %s" % args.input)

    if args.pass_only:
        filtering = "*"
    else:
        filtering = args.filters_only

    if args.bcf and not args.output.endswith(".bcf"):
        args.output += ".bcf"

    preprocess(args.input,
               args.output,
               args.ref,
               args.locations,
               filtering,
               args.fixchr,
               args.regions_bedfile,
               args.targets_bedfile,
               args.preprocessing_leftshift,
               args.preprocessing_decompose,
               args.preprocessing_norm,
               args.window,
               args.threads,
               args.gender,
               args.somatic_allele_conversion,
               convert_gvcf_to_vcf=args.convert_gvcf_to_vcf)

    elapsed = time.time() - starttime
    logging.info("preprocess for %s -- time taken %.2f" % (args.input, elapsed))


def updateArgs(parser):
    """ update command line parser with preprocessing args """

    parser.add_argument('--location', '-l', dest='locations', required=False, default=None,
                        help="Comma-separated list of locations [use naming after preprocessing], "
                             "when not specified will use whole VCF.")

    parser.add_argument("--pass-only", dest="pass_only", action="store_true", default=False,
                        help="Keep only PASS variants.")

    parser.add_argument("--filters-only", dest="filters_only", default="",
                        help="Specify a comma-separated list of filters to apply (by default all filters"
                             " are ignored / passed on.")

    parser.add_argument("-R", "--restrict-regions", dest="regions_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (sparse) regions (using -R in bcftools).")

    parser.add_argument("-T", "--target-regions", dest="targets_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (dense) regions (using -T in bcftools).")

    # preprocessing steps
    parser.add_argument("-L", "--leftshift", dest="preprocessing_leftshift", action="store_true",
                        default=True,
                        help="Left-shift variants safely.")
    parser.add_argument("--no-leftshift", dest="preprocessing_leftshift", action="store_false",
                        help="Do not left-shift variants safely.")

    parser.add_argument("--decompose", dest="preprocessing_decompose", action="store_true",
                        default=True,
                        help="Decompose variants into primitives. This results in more granular counts.")
    parser.add_argument("-D", "--no-decompose", dest="preprocessing_decompose", action="store_false",
                        help="Do not decompose variants into primitives.")

    parser.add_argument("--bcftools-norm", dest="preprocessing_norm", action="store_true", default=False,
                        help="Enable preprocessing through bcftools norm -c x -D (requires external "
                             " preprocessing to be switched on).")

    parser.add_argument("--fixchr", dest="fixchr", action="store_true", default=None,
                        help="Add chr prefix to VCF records where necessary (default: auto, attempt to match reference).")

    parser.add_argument("--no-fixchr", dest="fixchr", action="store_false",
                        help="Do not add chr prefix to VCF records (default: auto, attempt to match reference).")

    parser.add_argument("--bcf", dest="bcf", action="store_true", default=False,
                        help="Use BCF internally. This is the default when the input file"
                             " is in BCF format already. Using BCF can speed up temp file access, "
                             " but may fail for VCF files that have broken headers or records that "
                             " don't comply with the header.")

    parser.add_argument("--somatic", dest="somatic_allele_conversion", action="store_true", default=False,
                        help="Assume the input file is a somatic call file and squash "
                             "all columns into one, putting all FORMATs into INFO + use "
                             " half genotypes (see also --set-gt). This will replace all sample "
                             "columns and replace them with a single one.")

    parser.add_argument("--set-gt", dest="somatic_allele_conversion", choices=["half", "hemi", "het", "hom", "first"],
                        help="This is used to treat Strelka somatic files "
                             "Possible values for this parameter: "
                             "half / hemi / het / hom / half "
                             "to assign one of the following genotypes to the "
                             "resulting sample:  1 | 0/1 | 1/1 | ./1. This will replace all sample "
                             "columns and replace them with a single one.")

    parser.add_argument('--filter-nonref', dest='filter_nonref', action="store_true", default=False,
                        help='Remove any variants genotyped as <NON_REF>.')  
                        
    parser.add_argument('--convert-gvcf-to-vcf', dest='convert_gvcf_to_vcf', action="store_true", default=False,
                        help='Convert the input set from genome VCF format to a VCF before processing.')
                        
    # genotype handling on chrX.
    parser.add_argument("--gender", dest="gender", choices=["male", "female", "auto", "none"], default="auto",
                        help="Specify sex. This determines how haploid calls on chrX get treated: for male samples,"
                             " all non-ref calls (in the truthset only when running through hap.py) are given a 1/1 genotype.")


def main():
    parser = argparse.ArgumentParser("VCF preprocessor")

    # input
    parser.add_argument("input", help="VCF file to process.", default=[], nargs=1)
    parser.add_argument("output", help="Output filename.", default=[], nargs=1)

    updateArgs(parser)

    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Show version number and exit.")

    parser.add_argument("-r", "--reference", dest="ref", help="Specify a reference file.",
                        default=Tools.defaultReference())

    parser.add_argument("-w", "--window-size", dest="window",
                        default=10000, type=int,
                        help="Preprocessing window size (variants further apart than that size are not expected to interfere).")

    parser.add_argument("--threads", dest="threads",
                        default=multiprocessing.cpu_count(), type=int,
                        help="Number of threads to use.")

    if Tools.has_sge:
        parser.add_argument("--force-interactive", dest="force_interactive",
                            default=False, action="store_true",
                            help="Force running interactively (i.e. when JOB_ID is not in the environment)")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    args, unknown_args = parser.parse_known_args()

    if not Tools.has_sge:
        args.force_interactive = True

    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        level=loglevel)

    # remove some safe unknown args
    unknown_args = [x for x in unknown_args if x not in ["--force-interactive"]]
    if len(sys.argv) < 2 or len(unknown_args) > 0:
        if unknown_args:
            logging.error("Unknown arguments specified : %s " % str(unknown_args))
        parser.print_help()
        exit(0)

    if args.version:
        print "pre.py %s" % Tools.version  # noqa:E999
        exit(0)

    args.input = args.input[0]
    args.output = args.output[0]

    preprocessWrapper(args)


if __name__ == "__main__":
    main()
