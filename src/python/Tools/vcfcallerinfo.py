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

import tempfile
import itertools
import subprocess
import logging
import os
import json


class CallerInfo(object):
    """ Class for collecting caller info and version
    """

    def __init__(self):
        # callers and aligners are stored in tuples of three:
        # (caller/aligner, version, parameters)
        self.callers = []
        self.aligners = []

    def __repr__(self):
        return "aligners=[" + ",".join(["/".join(xx) for xx in self.aligners]) + "] " + \
               "callers=[" + ",".join(["/".join(xx) for xx in self.callers]) + "]"

    def asDict(self):
        kvd = ["name", "version", "parameters"]
        return {"aligners": [dict(y for y in itertools.izip(kvd, x)) for x in self.aligners],
                "callers": [dict(y for y in itertools.izip(kvd, x)) for x in self.callers]}

    def addVCF(self, vcfname):
        """ Add caller versions from a VCF
        :param vcfname: VCF file name
        """
        tf = tempfile.NamedTemporaryFile(delete=False)
        tf.close()
        vfh = {}
        try:
            sp = subprocess.Popen("vcfhdr2json '%s' '%s'" % (vcfname, tf.name),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            o, e = sp.communicate()

            if sp.returncode != 0:
                raise Exception("vcfhdr2json call failed: %s / %s" % (o, e))

            vfh = json.load(open(tf.name))
        finally:
            try:
                os.unlink(tf.name)
            except:
                pass

        cp = ['unknown', 'unknown', '']
        gatk_callers = ["haplotypecaller", "unifiedgenotyper", "mutect"]
        sent_callers = ["haplotyper"]
        source_found = False

        for hf in vfh["fields"]:
            try:
                k = hf["key"]
                if k == "source":
                    try:
                        cp[0] = str(hf["values"])
                    except:
                        cp[0] = hf["value"]
                    if cp[0].startswith("Platypus_Version_"):
                        cp[1] = cp[0][len("Platypus_Version_"):]
                        cp[0] = "Platypus"
                    source_found = True
                elif k == "source_version":
                    try:
                        cp[1] = str(hf["values"])
                    except:
                        cp[1] = hf["value"]
                    source_found = True
                elif k == "cmdline":
                    try:
                        cp[2] = str(hf["values"])
                    except:
                        cp[2] = hf["value"]
                    source_found = True
                elif k == "platypusOptions":
                    try:
                        cp[2] = str(hf["values"])
                    except:
                        cp[2] = hf["value"]
                    source_found = True
                elif k == "octopus":
                    # octopus doesn't add a version
                    self.callers.append(["octopus", "unknown", str(hf["values"])])
                elif k.startswith("GATKCommandLine"):
                    caller = "GATK"
                    try:
                        caller += "-" + hf["values"]["ID"]
                    except:
                        pass
                    version = "unknown"
                    try:
                        version = hf["values"]["Version"]
                    except:
                        pass
                    options = ""
                    try:
                        options = hf["values"]["CommandLineOptions"]
                    except:
                        pass
                    if any(g in caller.lower() for g in gatk_callers):
                        self.callers.append([caller, version, options])
                elif k.startswith("SentieonCommandLine"):
                    caller = "Sentieon"
                    try:
                        caller += "-" + hf["values"]["ID"]
                    except:
                        pass
                    version = "unknown"
                    try:
                        version = hf["values"]["Version"]
                    except:
                        pass
                    options = ""
                    if any(s in caller.lower() for s in sent_callers):
                        self.callers.append([caller, version])

            except:
                pass
        if source_found:
            self.callers.append(cp)

    def addBAM(self, bamfile):
        """ Extract aligner information from a BAM file
        :param bamfile: name of BAM file
        """
        sp = subprocess.Popen("samtools view -H '%s'" % bamfile,
                              shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        o, e = sp.communicate()

        if sp.returncode != 0:
            raise Exception("Samtools call failed: %s / %s" % (o, e))

        for line in o.split("\n"):
            if not line.startswith("@PG"):
                continue
            try:
                # noinspection PyTypeChecker
                x = dict(y.split(":", 1) for y in line.split("\t")[1:])
            except:
                logging.warn("Unable to parse SAM/BAM header line: %s" % line)
                continue
            cp = ['unknown', 'unknown', '']
            try:
                cp[0] = x['PN']
            except:
                try:
                    cp[0] = x['ID']
                    if "-" in cp[0]:
                        cp[0] = cp[0].split("-")[0]
                except:
                    pass
            try:
                cp[1] = x['VN']
            except:
                pass
            try:
                cp[2] = x['CL']
            except:
                pass

            self.aligners.append(cp)
