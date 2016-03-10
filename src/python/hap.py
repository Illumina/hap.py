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

scriptDir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(scriptDir, '..', 'lib', 'python27')))

import Tools
from Tools import vcfextract
from Tools.bcftools import preprocessVCF, bedOverlapCheck
from Tools.parallel import runParallel
from Tools.fastasize import fastaContigLengths
import Haplo.blocksplit
import Haplo.xcmp
import Haplo.vcfeval
import Haplo.quantify

import qfy


def main():
    parser = argparse.ArgumentParser("Haplotype Comparison")

    # input
    parser.add_argument('--location', '-l', dest='locations', required=False, default=None,
                        help='Add a location to the compare list (when not given, will use chr1-22, chrX, chrY).')

    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help="Show version number and exit.")

    parser.add_argument("-P", "--include-nonpass", dest="usefiltered", action="store_true", default=False,
                        help="Use to include failing query variants in comparison.")

    parser.add_argument("--include-nonpass-truth", dest="usefiltered_truth", action="store_true", default=False,
                        help="Include failing variants from the truth dataset.")

    parser.add_argument("-R", "--restrict-regions", dest="regions_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (sparse) regions (using -R in bcftools).")

    parser.add_argument("-T", "--target-regions", dest="targets_bedfile",
                        default=None, type=str,
                        help="Restrict analysis to given (dense) regions (using -T in bcftools).")

    parser.add_argument("-r", "--reference", dest="ref", default=None, help="Specify a reference file.")

    # output
    parser.add_argument("-o", "--report-prefix", dest="reports_prefix",
                        default=None,
                        help="Filename prefix for report output.")

    # DEPRECATED: we don't write bed files after 0.2.9
    parser.add_argument("-B", "--write-bed", dest="write_bed",
                        default=False, action="store_true",
                        help="This option is deprecated. BED files will not be written anymore.")

    # add quantification args
    qfy.updateArgs(parser)

    parser.add_argument("--scratch-prefix", dest="scratch_prefix",
                        default=None,
                        help="Directory for scratch files.")

    parser.add_argument("--keep-scratch", dest="delete_scratch",
                        default=True, action="store_false",
                        help="Filename prefix for scratch report output.")

    # detailed control of comparison
    parser.add_argument("--preprocess-truth", dest="preprocessing_truth", action="store_true", default=False,
                        help="Preprocess truth file using bcftools.")

    parser.add_argument("--external-preprocessing", dest="preprocessing", action="store_true", default=False,
                        help="Perform VCF preprocessing using bcftools.")

    parser.add_argument("--bcftools-norm", dest="preprocessing_norm", action="store_true", default=False,
                        help="Enable preprocessing through bcftools norm -c x -D (requires external "
                             " preprocessing to be switched on).")

    parser.add_argument("-N", "--numeric-chromosomes", dest="numeric_chrs", action="store_true", default=None,
                        help="Use numeric chromosome names for truth and query. This is a shortcut for "
                             "-l 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y "
                             "--no-fixchr-truth --no-fixchr-query")

    parser.add_argument("-C", "--no-numeric-chromosomes", dest="numeric_chrs", action="store_false",
                        help="Use chr-prefixed chromosome names for truth and query. This is a shortcut for "
                             "-l chr1,...,chrY"
                             "--fixchr-truth --fixchr-query")

    parser.add_argument("--fixchr-truth", dest="fixchr_truth", action="store_true", default=None,
                        help="Add chr prefix to truth file (default: auto).")

    parser.add_argument("--fixchr-query", dest="fixchr_query", action="store_true", default=None,
                        help="Add chr prefix to query file (default: auto).")

    parser.add_argument("--no-fixchr-truth", dest="fixchr_truth", action="store_false",
                        help="Disable chr replacement for truth (default: auto).")

    parser.add_argument("--no-fixchr-query", dest="fixchr_query", action="store_false",
                        help="Add chr prefix to query file (default: auto).")

    parser.add_argument("--partial-credit", dest="partial_credit", action="store_true", default=None,
                        help="give credit for partially matched variants. "
                             "this is equivalent to --internal-leftshift and --internal-preprocessing.")

    parser.add_argument("--no-partial-credit", dest="partial_credit", action="store_false", default=None,
                        help="Give credit for partially matched variants. "
                             "This is equivalent to --internal-leftshift and --no-internal-preprocessing.")

    parser.add_argument("--internal-leftshift", dest="int_preprocessing_ls", action="store_true", default=None,
                        help="Switch off xcmp's internal VCF leftshift preprocessing.")

    parser.add_argument("--internal-preprocessing", dest="int_preprocessing", action="store_true", default=None,
                        help="Switch off xcmp's internal VCF leftshift preprocessing.")

    parser.add_argument("--no-internal-leftshift", dest="int_preprocessing_ls", action="store_false", default=None,
                        help="Switch off xcmp's internal VCF leftshift preprocessing.")

    parser.add_argument("--no-internal-preprocessing", dest="int_preprocessing", action="store_false", default=None,
                        help="Switch off xcmp's internal VCF leftshift preprocessing.")

    parser.add_argument("--no-haplotype-comparison", dest="no_hc", action="store_true", default=False,
                        help="Disable haplotype comparison (only count direct GT matches as TP).")

    parser.add_argument("--unhappy", dest="unhappy", action="store_true", default=False,
                        help="Combination of --no-haplotype-comparison --no-internal-preprocessing "
                             "--no-internal-leftshift.")

    parser.add_argument("--no-auto-index", dest="auto_index", action="store_false", default=True,
                        help="Disable automatic index creation for input files. "
                             "The index is only necessary at this stage if we want to auto-detect locations. "
                             "When used with -l, and when it is known that there are variants at all given locations "
                             "this is not needed and can be switched off to save time.")

    parser.add_argument("-w", "--window-size", dest="window",
                        default=50, type=int,
                        help="Minimum distance between two variants such that they fall into different haplotype "
                             "blocks")

    parser.add_argument("--enumeration-threshold", dest="max_enum",
                        default=16768, type=int,
                        help="Enumeration threshold / maximum number of sequences to enumerate per block.")

    parser.add_argument("-e", "--expand-hapblocks", dest="hb_expand",
                        default=30, type=int,
                        help="Expand haplotype blocks by this many basepairs left and right.")
    parser.add_argument("--threads", dest="threads",
                        default=multiprocessing.cpu_count(), type=int,
                        help="Number of threads to use.")

    parser.add_argument("--engine", dest="engine",
                        default="xcmp", choices=["xcmp", "vcfeval"],
                        help="Comparison engine to use.")

    parser.add_argument("--engine-vcfeval-path", dest="engine_vcfeval", required=False,
                        help="This parameter should give the path to the \"rtg\" executable.")
    parser.add_argument("--engine-vcfeval-template", dest="engine_vcfeval_template", required=False,
                        help="Vcfeval needs the reference sequence formatted in its own file format "
                             "(SDF -- run rtg format -o ref.SDF ref.fa).")

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
        exit(0)

    if args.version:
        print "Hap.py %s" % Tools.version
        exit(0)

    if args.write_bed:
        logging.warn("The -B / --write-bed switches are deprecated in versions 0.2.9+, ")

    if args.roc:
        args.write_vcf = True

    # disable all clever matching
    if args.unhappy:
        args.int_preprocessing = False
        args.int_preprocessing_ls = False
        args.no_hc = True

    # Counting with partial credit
    elif args.partial_credit:
        # partial_credit switch is overridden by --no-* switches
        args.int_preprocessing = True
        args.int_preprocessing_ls = True
    elif args.partial_credit is None:
        # in the default setting, we enable partial credit but only override the
        # preprocessing settings if they haven't been specified
        if args.int_preprocessing is None:
            args.int_preprocessing = True
        if args.int_preprocessing_ls is None:
            args.int_preprocessing_ls = True
    elif args.partial_credit is not None:  # explicitly set to false
        args.int_preprocessing = False
        args.int_preprocessing_ls = True

    if args.int_preprocessing is None:
        args.int_preprocessing = False
    if args.int_preprocessing_ls is None:
        args.int_preprocessing_ls = False

    logging.info("Preprocessing settings: %s / %s / %s" % ("leftshift" if args.int_preprocessing_ls else "no-leftshift",
                                                           "splitting" if args.int_preprocessing else "raw calls",
                                                           "haplocompare" if not args.no_hc else "no-haplocompare"))

    # sanity-check regions bed file (HAP-57)
    if args.regions_bedfile:
        logging.info("Checking input regions.")
        if bedOverlapCheck(args.regions_bedfile):
            raise Exception("The regions bed file (specified using -R) has overlaps, this will not work with xcmp."
                            " You can either use -T, or run the file through bedtools merge")
        args.preprocessing_truth = True
        args.preprocessing = True

    if args.targets_bedfile or args.engine != "xcmp":
        args.preprocessing_truth = True
        args.preprocessing = True

    if args.fp_bedfile and not os.path.exists(args.fp_bedfile):
        raise Exception("FP/confident call region bed file does not exist.")

    tempfiles = []

    try:
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

        logging.info("Comparing %s and %s" % (args.vcf1, args.vcf2))

        # detect numeric chromosome names
        if args.numeric_chrs is None:
            cts = fastaContigLengths(args.ref)
            cts = set(cts.keys())
            numeric_names = set(map(str, range(1, 23)) + ["X", "Y", "M"])
            non_numeric_names = set(["chr" + x for x in numeric_names])
            numeric_names &= cts
            non_numeric_names &= cts
            numeric_names = len(list(numeric_names))
            non_numeric_names = len(list(non_numeric_names))
            if numeric_names != 0 and non_numeric_names == 0:
                args.numeric_chrs = True
                logging.info("Auto-detected numeric chromosome names")
            elif numeric_names == 0 and non_numeric_names != 0:
                args.numeric_chrs = False
                logging.info("Auto-detected chr-prefixed chromosome names")

        if args.numeric_chrs:
            args.fixchr_truth = False
            args.fixchr_query = False
        elif args.numeric_chrs is not None:
            args.fixchr_truth = True
            args.fixchr_query = True

        h1 = vcfextract.extractHeadersJSON(args.vcf1)
        if args.auto_index and not h1["tabix"]:
            logging.info("Creating indexed version of %s -- consider creating an index beforehand to save time here." %
                         args.vcf1)
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="truth.ix",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            tempfiles.append(vtf.name + ".tbi")
            args.vcf1 = Tools.bcftools.makeIndex(args.vcf1, vtf.name)
            h1 = vcfextract.extractHeadersJSON(args.vcf1)

        h2 = vcfextract.extractHeadersJSON(args.vcf2)
        if args.auto_index and not h2["tabix"]:
            logging.info("Creating indexed version of %s -- consider creating an index beforehand to save time here." %
                         args.vcf2)
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="query.ix",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            tempfiles.append(vtf.name + ".tbi")
            args.vcf2 = Tools.bcftools.makeIndex(args.vcf2, vtf.name)
            h2 = vcfextract.extractHeadersJSON(args.vcf2)

        ref_check = True
        try:
            happy_ref = args.ref
            v1r = [_h for _h in h1["fields"] if _h["key"] == "reference"]
            v2r = [_h for _h in h2["fields"] if _h["key"] == "reference"]
            if args.verbose:
                logging.info("References used: hap.py: %s / truth: %s / "
                             "query: %s" % (str(happy_ref), str(v1r), str(v2r)))

            v1_ref = ";".join([str(xxy["value"]) for xxy in v1r]).replace("file://", "")
            v2_ref = ";".join([str(xxy["value"]) for xxy in v2r]).replace("file://", "")

            if happy_ref == v1_ref and v1_ref == v2_ref:
                ref_check = True

            rids_vh = set()
            rids_v1 = set()
            rids_v2 = set()
            for refid in ["hg19", "hg38", "grc37", "grc38"]:
                if refid in happy_ref.lower():
                    rids_vh.add(refid)
                if refid in v1_ref.lower():
                    rids_v1.add(refid)
                if refid in v2_ref.lower():
                    rids_v2.add(refid)

            rids_v1 = sorted(list(rids_v1))
            rids_v2 = sorted(list(rids_v2))
            rids_vh = sorted(list(rids_vh))

            to_cmp = None
            if rids_v1:
                to_cmp = rids_v1
            if rids_v2:
                to_cmp = rids_v2
            if rids_vh:
                to_cmp = rids_vh
            if to_cmp and rids_v1 and rids_v1 != to_cmp:
                ref_check = False
            if to_cmp and rids_v2 and rids_v2 != to_cmp:
                ref_check = False
            if to_cmp and rids_vh and rids_vh != to_cmp:
                ref_check = False

        except:
            pass

        if not ref_check:
            logging.warn("Reference sequence check failed! "
                         "Please ensure that truth and query VCF use the same reference sequence as "
                         "hap.py. XCMP may fail if this is not the case, and the results will not be "
                         " accurate.")

        if args.locations is None or len(args.locations) == 0:
            # all chromosomes
            if args.numeric_chrs:
                args.locations = [x for x in map(str, range(1, 23))]
            else:
                args.locations = ["chr" + x for x in map(str, range(1, 23))]

        if type(args.locations) is not list and args.locations is not None:
            # noinspection PyUnresolvedReferences
            args.locations = args.locations.split(",")

        # HAP-143 fix the case where no chromosomes are in truth or query
        try:
            if not h1["tabix"]["chromosomes"]:
                h1["tabix"]["chromosomes"] = []
        except:
            pass
        try:
            if not h2["tabix"]["chromosomes"]:
                h2["tabix"]["chromosomes"] = []
        except:
            pass

        if not h1["tabix"]:
            args.preprocessing_truth = True
            logging.warn("Truth file is not Tabix indexed. Switching on pre-processing + chr name conversion.")
            if args.fixchr_truth is None:
                args.fixchr_truth = True
        elif args.fixchr_truth is None:
            logging.info(str(h1["tabix"]))
            # autodetect chr naming
            count_with_fix = len([__ for __ in h1["tabix"]["chromosomes"]
                                  if ("chr%s" % str(__)) in args.locations])
            count_no_fix = len([__ for __ in h1["tabix"]["chromosomes"] if str(__) in args.locations])
            logging.info("Truth: Number of chromosome names matching with / without renaming : %i / %i " % (
                count_with_fix, count_no_fix))
            if count_with_fix > count_no_fix:
                args.fixchr_truth = True
                logging.info("Will fix chromosome names (truth).")
            else:
                logging.info("Will not fix chromosome names (truth).")
                args.fixchr_truth = False

        if not h2["tabix"]:
            args.preprocessing = True
            logging.warn("Query file is not Tabix indexed. Switching on pre-processing + chr name conversion.")
            # don't overwrite setting, but if it's None, replace with True to be sure
            if args.fixchr_query is None:
                args.fixchr_query = True
        elif args.fixchr_query is None:
            # autodetect chr naming
            count_with_fix = len([__ for __ in h2["tabix"]["chromosomes"]
                                  if ("chr%s" % str(__)) in args.locations])
            count_no_fix = len([__ for __ in h2["tabix"]["chromosomes"] if str(__) in args.locations])
            logging.info("Query: Number of chromosome names matching with / without renaming : %i / %i " % (
                count_with_fix, count_no_fix))
            if count_with_fix > count_no_fix:
                args.fixchr_query = True
                logging.info("Will fix chromosome names (query).")
            else:
                logging.info("Will not fix chromosome names (query).")
                args.fixchr_query = False

        if args.fixchr_truth or args.preprocessing_norm:
            args.preprocessing_truth = True

        if args.fixchr_query or args.preprocessing_norm:
            args.preprocessing = True

        if args.preprocessing_truth:
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="truth.pp",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            preprocessVCF(args.vcf1, vtf.name, ",".join(args.locations),
                          not args.usefiltered_truth,  # pass_only
                          args.fixchr_truth,  # chrprefix
                          args.preprocessing_norm,  # norm,
                          args.regions_bedfile,
                          args.targets_bedfile,
                          args.ref)
            args.vcf1 = vtf.name
            # get headers again if we preprocessed
            h1 = vcfextract.extractHeadersJSON(args.vcf1)

        if args.preprocessing:
            vtf = tempfile.NamedTemporaryFile(delete=False,
                                              dir=args.scratch_prefix,
                                              prefix="query.pp",
                                              suffix=".vcf.gz")
            vtf.close()
            tempfiles.append(vtf.name)
            preprocessVCF(args.vcf2, vtf.name, ",".join(args.locations),
                          False,  # query filters are handled further down in matching
                          args.fixchr_query,  # chrprefix
                          args.preprocessing_norm,  # norm,
                          args.regions_bedfile,
                          args.targets_bedfile,
                          args.ref)
            args.vcf2 = vtf.name
            # get headers again if we preprocessed
            h2 = vcfextract.extractHeadersJSON(args.vcf2)

        if not h1["tabix"]:
            raise Exception("Truth file is not Tabix indexed.")

        if not h2["tabix"]:
            raise Exception("Query file is not Tabix indexed.")

        newlocations = []

        if not h1["tabix"]["chromosomes"]:
            h1["tabix"]["chromosomes"] = []
        if not h2["tabix"]["chromosomes"]:
            h2["tabix"]["chromosomes"] = []

        for _xc in args.locations:
            xc = _xc.split(":")[0]
            if xc not in h1["tabix"]["chromosomes"]:
                logging.warn("No calls for location %s in truth!" % xc)
            if xc not in h2["tabix"]["chromosomes"]:
                logging.warn("No calls for location %s in query!" % xc)

            if (xc not in h1["tabix"]["chromosomes"]) and (xc not in h2["tabix"]["chromosomes"]):
                logging.warn("Removing location %s because neither input file has calls there." % xc)
            else:
                newlocations.append(_xc)

        if not newlocations:
            raise Exception("Location list is empty: the input files do not appear to have variants on any of %s" %
                            str(args.locations))

        args.locations = newlocations

        if args.threads > 1 and args.engine == "xcmp":
            logging.info("Running using %i parallel processes." % args.threads)
            pool = multiprocessing.Pool(int(args.threads))

            # find balanced pieces
            args.pieces = (args.threads + len(args.locations) - 1) / len(args.locations)
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
        else:
            pool = None

        # count variants before normalisation
        if "samples" not in h1 or not h1["samples"]:
            raise Exception("Cannot read sample names from truth VCF file")

        if "samples" not in h2 or not h2["samples"]:
            raise Exception("Cannot read sample names from query VCF file")

        tf = tempfile.NamedTemporaryFile(delete=False,
                                         dir=args.scratch_prefix,
                                         prefix="hap.py.result.", suffix=".vcf.gz")
        tf.close()
        tempfiles.append(tf.name)
        output_name = tf.name

        if args.engine == "xcmp":
            # do xcmp
            logging.info("Using xcmp for comparison")
            res = runParallel(pool, Haplo.xcmp.xcmpWrapper, args.locations, args)
            tempfiles += [x[0] for x in res if x is not None]  # VCFs
            tempfiles += [x[1] for x in res if x is not None and x[1] is not None]  # beds (if any)

            if None in res:
                raise Exception("One of the xcmp jobs failed.")

            if len(res) == 0:
                raise Exception("Input files/regions do not contain variants (0 haplotype blocks were processed).")

            # concatenate + index
            logging.info("Concatenating variants...")
            runme_list = [x[0] for x in res if x is not None]
            if len(runme_list) == 0:
                raise Exception("No outputs to concatenate!")

            fo = Tools.BGZipFile(output_name, True)
            for i, x in enumerate(runme_list):
                f = gzip.GzipFile(x)
                for l in f:
                    if i == 0 or not l[0] == "#":
                        fo.write(l)
            fo.close()

            logging.info("Indexing...")
            to_run = "tabix -p vcf %s" % output_name.replace(" ", "\\ ")
            logging.info("Running '%s'" % to_run)
            subprocess.check_call(to_run, shell=True)
            # passed to quantify
            args.type = "xcmp"
        elif args.engine == "vcfeval":
            tempfiles += Haplo.vcfeval.runVCFEval(args.vcf1, args.vcf2, output_name, args)
            # passed to quantify
            args.type = "ga4gh"
        else:
            raise Exception("Unknown comparison engine: %s" % args.engine)

        logging.info("Counting variants...")

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
