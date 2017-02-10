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
# Diploid VCF File Comparison
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
import traceback
import subprocess
import multiprocessing
import gzip
import tempfile
import time

scriptDir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools import vcfextract
from Tools import bcftools
from Tools.parallel import runParallel, getPool
from Tools.bcftools import preprocessVCF, bedOverlapCheck
from Tools.fastasize import fastaContigLengths
import Haplo.blocksplit
import Haplo.xcmp
import Haplo.vcfeval
import Haplo.scmp
import Haplo.quantify
import Haplo.partialcredit
import Haplo.gvcf2bed

import qfy
import pre


def main():
    parser = argparse.ArgumentParser("Haplotype Comparison")

    # input
    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Show version number and exit.")

    parser.add_argument("-r", "--reference", dest="ref", default=None, help="Specify a reference file.")

    # output
    parser.add_argument("-o", "--report-prefix", dest="reports_prefix",
                        default=None,
                        help="Filename prefix for report output.")
    parser.add_argument("--scratch-prefix", dest="scratch_prefix",
                        default=None,
                        help="Directory for scratch files.")
    parser.add_argument("--keep-scratch", dest="delete_scratch",
                        default=True, action="store_false",
                        help="Filename prefix for scratch report output.")


    # add quantification args
    qfy.updateArgs(parser)

    # control preprocessing
    pre.updateArgs(parser)
    parser.add_argument("--preprocess-truth", dest="preprocessing_truth", action="store_true", default=False,
                        help="Preprocess truth file with same settings as query (default is to accept truth in original format).")
    parser.add_argument("--usefiltered-truth", dest="usefiltered_truth", action="store_true", default=False,
                        help="Use filtered variant calls in truth file (by default, only PASS calls in the truth file are used)")
    parser.add_argument("--preprocessing-window-size", dest="preprocess_window",
                        default=10000, type=int,
                        help="Preprocessing window size (variants further apart than that size are not expected to interfere).")
    parser.add_argument("--adjust-conf-regions", dest="preprocessing_truth_confregions", action="store_true", default=True,
                        help="Adjust confident regions to include variant locations.")
    parser.add_argument("--no-adjust-conf-regions", dest="preprocessing_truth_confregions", action="store_false",
                        help="Adjust confident regions to include variant locations.")

    parser.add_argument("--quantify-homref-calls", dest="quantify_homref_calls", action="store_true", default=False,
                        help="When a confident region file is given we can also quantify the homref calls "
                             "(i.e. \"true negatives\") on a gVCF query file.")

    # detailed control of comparison
    parser.add_argument("--unhappy", "--no-haplotype-comparison", dest="no_hc", action="store_true", default=False,
                        help="Disable haplotype comparison (only count direct GT matches as TP).")

    parser.add_argument("-w", "--window-size", dest="window",
                        default=50, type=int,
                        help="Minimum distance between variants such that they fall into the same superlocus.")

    # xcmp-specific stuff
    parser.add_argument("--xcmp-enumeration-threshold", dest="max_enum",
                        default=16768, type=int,
                        help="Enumeration threshold / maximum number of sequences to enumerate per block.")

    parser.add_argument("--xcmp-expand-hapblocks", dest="hb_expand",
                        default=30, type=int,
                        help="Expand haplotype blocks by this many basepairs left and right.")
    parser.add_argument("--threads", dest="threads",
                        default=multiprocessing.cpu_count(), type=int,
                        help="Number of threads to use.")

    parser.add_argument("--engine", dest="engine",
                        default="xcmp", choices=["xcmp", "vcfeval", "scmp-somatic"],
                        help="Comparison engine to use.")

    parser.add_argument("--engine-vcfeval-path", dest="engine_vcfeval", required=False,
                        default=Haplo.vcfeval.findVCFEval(),
                        help="This parameter should give the path to the \"rtg\" executable. "
                             "The default is %s" % Haplo.vcfeval.findVCFEval())
    parser.add_argument("--engine-vcfeval-template", dest="engine_vcfeval_template", required=False,
                        help="Vcfeval needs the reference sequence formatted in its own file format "
                             "(SDF -- run rtg format -o ref.SDF ref.fa). You can specify this here "
                             "to save time when running hap.py with vcfeval. If no SDF folder is "
                             "specified, hap.py will create a temporary one.")

    if Tools.has_sge:
        parser.add_argument("--force-interactive", dest="force_interactive",
                            default=False, action="store_true",
                            help="Force running interactively (i.e. when JOB_ID is not in the environment)")

    parser.add_argument("_vcfs", help="Two VCF files.", default=[], nargs="*")

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
        exit(1)

    if args.version:
        print "Hap.py %s" % Tools.version
        exit(0)

    if args.roc:
        args.write_vcf = True

    # sanity-check regions bed file (HAP-57)
    if args.regions_bedfile:
        logging.info("Checking input regions.")
        if bedOverlapCheck(args.regions_bedfile):
            raise Exception("The regions bed file (specified using -R) has overlaps, this will not work with xcmp."
                            " You can either use -T, or run the file through bedtools merge")

    if args.fp_bedfile and not os.path.exists(args.fp_bedfile):
        raise Exception("FP/confident call region bed file does not exist.")

    if not args.force_interactive and "JOB_ID" not in os.environ:
        parser.print_help()
        raise Exception("Please qsub me so I get approximately 1 GB of RAM per thread.")

    if not args.ref:
        args.ref = Tools.defaultReference()

    if not os.path.exists(args.ref):
        raise Exception("Please specify a valid reference path using -r.")

    if not args.reports_prefix:
        raise Exception("Please specify an output prefix using -o ")

    if not os.path.exists(os.path.dirname(os.path.abspath(args.reports_prefix))):
        raise Exception("The output path does not exist. Please specify a valid output path and prefix using -o")

    if os.path.basename(args.reports_prefix) == "" or os.path.isdir(args.reports_prefix):
        raise Exception("The output path should specify a file name prefix. Please specify a valid output path "
                        "and prefix using -o. For example, -o /tmp/test will create files named /tmp/test* .")

    # noinspection PyProtectedMember
    if not args._vcfs or len(args._vcfs) != 2:
        raise Exception("Please specify exactly two input VCFs.")

    # noinspection PyProtectedMember
    args.vcf1 = args._vcfs[0]
    # noinspection PyProtectedMember
    args.vcf2 = args._vcfs[1]

    if not os.path.exists(args.vcf1):
        raise Exception("Input file %s does not exist." % args.vcf1)
    if not os.path.exists(args.vcf2):
        raise Exception("Input file %s does not exist." % args.vcf2)

    tempfiles = []

    # turn on allele conversion
    if args.engine == "scmp-somatic" and args.somatic_allele_conversion == False:
        args.somatic_allele_conversion = True

    # somatic allele conversion should also switch off decomposition
    if args.somatic_allele_conversion == True and "--decompose" not in sys.argv:
        args.preprocessing_decompose = False

    # xcmp/scmp support bcf; others don't
    if args.engine in ["xcmp", "scmp-somatic"] and (args.bcf or (args.vcf1.endswith(".bcf") and args.vcf2.endswith(".bcf"))):
        internal_format_suffix = ".bcf"
    else:
        internal_format_suffix = ".vcf.gz"

    try:
        logging.info("Comparing %s and %s" % (args.vcf1, args.vcf2))

        logging.info("Preprocessing truth: %s" % args.vcf1)
        starttime = time.time()

        ttf = tempfile.NamedTemporaryFile(delete=False,
                                          dir=args.scratch_prefix,
                                          prefix="truth.pp",
                                          suffix=internal_format_suffix)
        ttf.close()

        if args.engine.endswith("somatic") and \
           args.preprocessing_truth and \
           (args.preprocessing_leftshift or args.preprocessing_norm or args.preprocessing_decompose):
            args.preprocessing_truth = False
            logging.info("Turning off pre.py preprocessing for somatic comparisons")

        tempfiles.append(ttf.name)
        tempfiles.append(ttf.name + ".csi")
        tempfiles.append(ttf.name + ".tbi")
        args.gender = pre.preprocess(args.vcf1,
                                     ttf.name,
                                     args.ref,
                                     args.locations,
                                     None if args.usefiltered_truth else "*",  # filters
                                     args.fixchr,
                                     args.regions_bedfile,
                                     args.targets_bedfile,
                                     args.preprocessing_leftshift if args.preprocessing_truth else False,
                                     args.preprocessing_decompose if args.preprocessing_truth else False,
                                     args.preprocessing_norm if args.preprocessing_truth else False,
                                     args.preprocess_window,
                                     args.threads,
                                     args.gender,
                                     args.somatic_allele_conversion)

        args.vcf1 = ttf.name

        if args.fp_bedfile and args.preprocessing_truth_confregions:
            conf_temp = Haplo.gvcf2bed.gvcf2bed(args.vcf1, args.ref, args.fp_bedfile, args.scratch_prefix)
            tempfiles.append(conf_temp)
            args.strat_regions.append("CONF_VARS:" + conf_temp)

        h1 = vcfextract.extractHeadersJSON(args.vcf1)

        elapsed = time.time() - starttime
        logging.info("preprocess for %s -- time taken %.2f" % (args.vcf1, elapsed))

        # once we have preprocessed the truth file we can resolve the locations
        # doing this here improves the time for query preprocessing below
        reference_contigs = set(fastaContigLengths(args.ref).keys())

        if not args.locations:
            # default set of locations is the overlap between truth and reference
            args.locations = list(reference_contigs & set(h1["tabix"]["chromosomes"]))
            if not args.locations:
                raise Exception("Truth and reference have no chromosomes in common!")
        elif type(args.locations) is not list:
            args.locations = args.locations.split(",")

        args.locations = sorted(args.locations)

        logging.info("Preprocessing query: %s" % args.vcf2)
        starttime = time.time()

        if args.pass_only:
            filtering = "*"
        else:
            filtering = args.filters_only

        qtf = tempfile.NamedTemporaryFile(delete=False,
                                          dir=args.scratch_prefix,
                                          prefix="query.pp",
                                          suffix=internal_format_suffix)
        qtf.close()
        tempfiles.append(qtf.name)
        tempfiles.append(qtf.name + ".csi")
        tempfiles.append(qtf.name + ".tbi")

        if args.engine.endswith("somatic") and \
           (args.preprocessing_leftshift or args.preprocessing_norm or args.preprocessing_decompose):
            args.preprocessing_leftshift = False
            args.preprocessing_norm = False
            args.preprocessing_decompose = False
            logging.info("Turning off pre.py preprocessing (query) for somatic comparisons")

        pre.preprocess(args.vcf2,
                       qtf.name,
                       args.ref,
                       str(",".join(args.locations)),
                       filtering,
                       args.fixchr,
                       args.regions_bedfile,
                       args.targets_bedfile,
                       args.preprocessing_leftshift,
                       args.preprocessing_decompose,
                       args.preprocessing_norm,
                       args.preprocess_window,
                       args.threads,
                       args.gender,
                       args.somatic_allele_conversion)  # same gender as truth above

        args.vcf2 = qtf.name
        h2 = vcfextract.extractHeadersJSON(args.vcf2)

        if args.quantify_homref_calls:
            conf_temp = Haplo.gvcf2bed.gvcf2bed(args.vcf2, args.ref, None, args.scratch_prefix)
            tempfiles.append(conf_temp)
            args.strat_regions.append("QUERY_HOMREF:" + conf_temp)

        elapsed = time.time() - starttime
        logging.info("preprocess for %s -- time taken %.2f" % (args.vcf2, elapsed))

        if not h1["tabix"]:
            raise Exception("Truth file is not indexed after preprocesing.")

        if not h2["tabix"]:
            raise Exception("Query file is not indexed after preprocessing.")

        for _xc in args.locations:
            if _xc not in h2["tabix"]["chromosomes"]:
                logging.warn("No calls for location %s in query!" % _xc)

        pool = getPool(args.threads)
        if args.threads > 1 and args.engine == "xcmp":
            logging.info("Running using %i parallel processes." % args.threads)

            # find balanced pieces
            # cap parallelism at 64 since otherwise bcftools concat below might run out
            # of file handles
            args.pieces = min(args.threads, 64)
            res = runParallel(pool, Haplo.blocksplit.blocksplitWrapper, args.locations, args)

            if None in res:
                raise Exception("One of the blocksplit processes failed.")

            tempfiles += res

            args.locations = []
            for f in res:
                with open(f) as fp:
                    for l in fp:
                        ll = l.strip().split("\t", 3)
                        if len(ll) < 3:
                            continue
                        xchr = ll[0]
                        start = int(ll[1]) + 1
                        end = int(ll[2])
                        args.locations.append("%s:%i-%i" % (xchr, start, end))

        # count variants before normalisation
        if "samples" not in h1 or not h1["samples"]:
            raise Exception("Cannot read sample names from truth VCF file")

        if "samples" not in h2 or not h2["samples"]:
            raise Exception("Cannot read sample names from query VCF file")

        tf = tempfile.NamedTemporaryFile(delete=False,
                                         dir=args.scratch_prefix,
                                         prefix="hap.py.result.",
                                         suffix=internal_format_suffix)
        tf.close()
        tempfiles.append(tf.name)
        tempfiles.append(tf.name + ".tbi")
        tempfiles.append(tf.name + ".csi")
        output_name = tf.name

        if args.engine == "xcmp":
            # do xcmp
            logging.info("Using xcmp for comparison")
            res = runParallel(pool, Haplo.xcmp.xcmpWrapper, args.locations, args)
            tempfiles += [x for x in res if x is not None]  # VCFs

            if None in res:
                raise Exception("One of the xcmp jobs failed.")

            if len(res) == 0:
                raise Exception("Input files/regions do not contain variants (0 haplotype blocks were processed).")

            # concatenate + index
            logging.info("Concatenating variants...")
            runme_list = [x for x in res if x is not None]
            if len(runme_list) == 0:
                raise Exception("No outputs to concatenate!")

            logging.info("Concatenating...")
            bcftools.concatenateParts(output_name, *runme_list)
            logging.info("Indexing...")
            bcftools.runBcftools("index", output_name)
            # passed to quantify
            args.type = "xcmp"
            # xcmp extracts whichever field we're using into the QQ info field
            args.roc_header = args.roc
            args.roc = "IQQ"
        elif args.engine == "vcfeval":
            tempfiles += Haplo.vcfeval.runVCFEval(args.vcf1, args.vcf2, output_name, args)
            # passed to quantify
            args.type = "ga4gh"
        elif args.engine == "scmp-somatic":
            tempfiles += Haplo.scmp.runSCmp(args.vcf1, args.vcf2, output_name, args)
            # passed to quantify
            args.type = "ga4gh"
        else:
            raise Exception("Unknown comparison engine: %s" % args.engine)

        args.in_vcf = [output_name]
        args.runner = "hap.py"
        qfy.quantify(args)

    finally:
        if args.delete_scratch:
            for x in tempfiles:
                try:
                    os.remove(x)
                except:
                    pass
        else:
            logging.info("Scratch files kept : %s" % (str(tempfiles)))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        traceback.print_exc(file=Tools.LoggingWriter(logging.ERROR))
        exit(1)
