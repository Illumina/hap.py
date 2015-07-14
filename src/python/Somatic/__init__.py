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

import abc
import pandas

from Tools.vcfextract import vcfExtract

from Strelka import extractStrelkaSNVFeatures, extractStrelkaIndelFeatures
from Mutect import extractMutectSNVFeatures, extractMutectIndelFeatures


class FeatureSet(object):
    """ VCF paired Feature set for somatic comparison """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.chr_depth = {}

    @abc.abstractmethod
    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        pass

    sets = {}

    @staticmethod
    def register(name, xclass):
        FeatureSet.sets[name] = xclass

    @staticmethod
    def make(name):
        # noinspection PyCallingNonCallable
        return FeatureSet.sets[name]()

    def setChrDepths(self, cd):
        """ set depth normalisation factors (can come from VCF or BAM) """
        self.chr_depth = cd


class GenericFeatures(FeatureSet):
    """ Collect generic variant features """

    features = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"]

    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        return GenericFeatures.collectFeatures(vcfname, tag, GenericFeatures.features)

    @staticmethod
    def processValue(t):
        _, val = t
        if type(val) is list:
            return ",".join(map(str, val))
        else:
            return val

    @staticmethod
    def collectFeatures(vcfname, tag, features, processor=None):
        if not processor:
            processor = GenericFeatures.processValue

        records = []

        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, v in enumerate(vr):
                rec[features[i]] = processor((features[i], v))
            rec["tag"] = tag
            records.append(rec)

        if records:
            df = pandas.DataFrame(records, columns=features + ["tag"])
        else:
            df = pandas.DataFrame(columns=features + ["tag"])
        return df


FeatureSet.register("generic", GenericFeatures)


class StrelkaAdmixSNVFeatures(FeatureSet):
    """ Collect SNV features from Strelka-to-admixture comparison """

    @staticmethod
    def processValue(t):
        n, val = t
        if type(val) is list:
            return ",".join(map(str, val))
        return val


    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        if tag not in ["TP", "FN"]:
            return extractStrelkaSNVFeatures(vcfname, tag, self.chr_depth)
        else:
            features = ["CHROM", "POS", "REF", "ALT", "I.editDistance", "S.2.GT"]
            return GenericFeatures.collectFeatures(vcfname, tag, features, processor=StrelkaAdmixSNVFeatures.processValue)

FeatureSet.register("admix.strelka.snv", StrelkaAdmixSNVFeatures)


class StrelkaAdmixIndelFeatures(StrelkaAdmixSNVFeatures):
    """ Collect Indel features from Strelka-to-admixture comparison """

    @staticmethod
    def processValue(t):
        n, val = t
        if n == "I.SGT":
            try:
                if val == "ref->het":
                    val = "het"
                elif val == "ref->hom":
                    val = "loh"
                else:
                    val = "unknown"
            except:
                val = "unknown"
        else:
            if type(val) is list:
                return ",".join(map(str, val))
        return val


    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        if tag not in ["TP", "FN"]:
            return extractStrelkaIndelFeatures(vcfname, tag, self.chr_depth)
        else:
            features = ["CHROM", "POS", "REF", "ALT", "I.editDistance", "S.2.GT"]
        return GenericFeatures.collectFeatures(vcfname, tag, features, processor=StrelkaAdmixIndelFeatures.processValue)

FeatureSet.register("admix.strelka.indel", StrelkaAdmixIndelFeatures)


class StrelkaHCCSNVFeatures(StrelkaAdmixSNVFeatures):
    """ Collect SNV features from Strelka-to-HCC truthset comparison """
    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        if tag not in ["TP", "FN"]:
            return extractStrelkaSNVFeatures(vcfname, tag, self.chr_depth)
        else:
            features = ["CHROM", "POS", "REF", "ALT", "QUAL",
                        "I.T_ALT_RATE", "I.DP_normal", "I.DP_tumor", "I.tag", "I.count"]
        return GenericFeatures.collectFeatures(vcfname, tag, features, processor=StrelkaAdmixSNVFeatures.processValue)

FeatureSet.register("hcc.strelka.snv", StrelkaHCCSNVFeatures)


class StrelkaHCCIndelFeatures(StrelkaAdmixIndelFeatures):
    """ Collect Indel features from Strelka-to-HCC truthset comparison """
    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        if tag not in ["TP", "FN"]:
            return extractStrelkaIndelFeatures(vcfname, tag, self.chr_depth)
        else:
            features = ["CHROM", "POS", "REF", "ALT", "QUAL", "S.1.VT",
                        "I.T_ALT_RATE", "I.DP_normal", "I.DP_tumor", "I.tag", "I.count"]
        return GenericFeatures.collectFeatures(vcfname, tag, features, processor=StrelkaAdmixIndelFeatures.processValue)

FeatureSet.register("hcc.strelka.indel", StrelkaHCCIndelFeatures)


class MutectHCCSNVFeatures(StrelkaAdmixSNVFeatures):
    """ Collect SNV features from Strelka-to-HCC truthset comparison """
    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        if tag not in ["TP", "FN"]:
            return extractMutectSNVFeatures(vcfname, tag, self.chr_depth)
        else:
            features = ["CHROM", "POS", "REF", "ALT", "QUAL",
                        "I.MapQrange", "I.somatic", "I.filtered", "S.1.VT",
                        "I.T_ALT_RATE", "I.DP_normal", "I.DP_tumor", "I.tag", "I.count"]
        return GenericFeatures.collectFeatures(vcfname, tag, features, processor=StrelkaAdmixSNVFeatures.processValue)

FeatureSet.register("hcc.mutect.snv", MutectHCCSNVFeatures)


class MutectHCCIndelFeatures(StrelkaAdmixIndelFeatures):
    """ Collect Indel features from Strelka-to-HCC truthset comparison """
    def collect(self, vcfname, tag):
        """ Return a data frame with features collected from the given VCF, tagged by given type """
        if tag not in ["TP", "FN"]:
            return extractMutectIndelFeatures(vcfname, tag, self.chr_depth)
        else:
            features = ["CHROM", "POS", "REF", "ALT", "QUAL", "S.1.VT",
                        "I.T_ALT_RATE", "I.DP_normal", "I.DP_tumor", "I.tag", "I.count"]
        return GenericFeatures.collectFeatures(vcfname, tag, features, processor=StrelkaAdmixIndelFeatures.processValue)

FeatureSet.register("hcc.mutect.indel", MutectHCCIndelFeatures)
