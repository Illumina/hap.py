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
import errno
import logging
import subprocess
from datetime import date

# noinspection PyUnresolvedReferences
try:
    import Haplo.version as vs  # pylint: disable=E0401,E0611
    version = vs.__version__
    has_sge = vs.has_sge
except ImportError:
    logging.warn("No version found. Please follow the installation instructions.")
    version = "unknown"
    has_sge = False


def defaultReference():
    to_try = ['/opt/hap.py-data/hg19.fa']
    try:
        to_try.insert(0, os.environ["HGREF"])
    except:
        pass

    try:
        to_try.insert(0, os.environ["HG19"])
    except:
        pass

    for x in to_try:
        if os.path.exists(x):
            return x
    logging.warn("No reference file found at default locations. You can set the environment"
                 " variable 'HGREF' or 'HG19' to point to a suitable Fasta file.")
    return None


def which(program):
    def is_exe(xfpath):
        return os.path.isfile(xfpath) and os.access(xfpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def init():
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)

    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    paths = ["bin"]

    for p in paths:
        pp = os.path.join(base, p)
        if not os.path.exists(pp):
            raise Exception("Dependency path %s not found" % pp)
        os.environ["PATH"] = pp + os.pathsep + os.environ["PATH"]

    executables = [
        "blocksplit",
        "hapenum",
        "dipenum",
        "hapcmp",
        "xcmp",
        "bcftools",
        "samtools",
    ]

    for x in executables:
        if not which(x):
            raise Exception("Dependency %s not found" % x)

    os.environ['DYLD_LIBRARY_PATH'] = os.path.join(base, "lib")
    os.environ['LD_LIBRARY_PATH'] = os.path.join(base, "lib")


init()

# safely import here
# noinspection PyUnresolvedReferences
import pysam  # noqa: F401
# noinspection PyUnresolvedReferences
import pandas  # noqa: F401


class LoggingWriter(object):
    """ Helper class to write tracebacks to log file
    """
    def __init__(self, level):
        self.level = level

    def write(self, message):
        message = message.replace("\n", "")
        if message:
            logging.log(self.level, message)


def writeVCFHeader(filelike, extrainfo="", chrprefix="chr"):
    """ write a VCF header
    """
    header = ["CHROM",
              "POS",
              "ID",
              "REF",
              "ALT",
              "QUAL",
              "FILTER",
              "INFO",
              "FORMAT",
              "SIMPLE"]

    infos = ['##fileformat=VCFv4.1',
             '##reference=hg19',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']

    if extrainfo:
        if type(extrainfo) is list:
            infos += extrainfo
        else:
            infos += extrainfo.split("\n")

    contigs = [["1", "249250621"], ["2", "243199373"], ["3", "198022430"], ["4", "191154276"],
               ["5", "180915260"], ["6", "171115067"], ["7", "159138663"], ["8", "146364022"],
               ["9", "141213431"], ["10", "135534747"], ["11", "135006516"],
               ["12", "133851895"], ["13", "115169878"], ["14", "107349540"], ["15", "102531392"],
               ["16", "90354753"], ["17", "81195210"], ["18", "78077248"], ["19", "59128983"],
               ["20", "63025520"], ["21", "48129895"], ["22", "51304566"], ["X", "155270560"]]

    meta = ["##fileDate=" + date.today().isoformat(),
            "##source=HaploCompare",
            "##source_version=%s" % version]

    for i in infos:
        filelike.write(i + "\n")
    for c in contigs:
        filelike.write('##contig=<ID=' + chrprefix + c[0] + ',length=' + c[1] + '>\n')
    for i in meta:
        filelike.write(i + "\n")

    filelike.write("#" + "\t".join(header) + "\n")


def mkdir_p(path):
    """ mkdir -p path """
    try:
        os.makedirs(os.path.abspath(path))
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    if not os.path.isdir(path):
        raise Exception("Failed to create directory %s" % path)


class BGZipFile(object):
    """ BGZip file helper
    """

    def __init__(self, filename, force=False):
        """ Make a subprocess for bgzip
        :param filename: name of the output file
        :param force: true to overwrite if file exists
        """
        if os.path.exists(filename) and not force:
            raise Exception("File %s exists, use force=True to overwrite" % filename)

        self.write_file = open(filename, "wb")
        zip_pipe = subprocess.Popen(["bgzip", "-f"],
                                    stdin=subprocess.PIPE,
                                    stdout=self.write_file,
                                    stderr=subprocess.PIPE,
                                    shell=True)
        self.zip_pipe = zip_pipe
        self.name = filename

    def close(self):
        self.zip_pipe.stdin.flush()
        self.zip_pipe.stdin.close()
        self.zip_pipe.wait()
        self.write_file.flush()
        self.write_file.close()

    def write(self, *args, **kwargs):
        self.zip_pipe.stdin.write(*args, **kwargs)
