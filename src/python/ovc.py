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
# 10/1/2015
#
# Testing helper, will return exit code 1 if the passed bed file has
# overlapping intervals.
#
# Usage:
#
# ovc.py input.bed
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import sys

f = open(sys.argv[1])

last = -1

lines = 1
for line in f:
	l = line.split("\t")
	if len(l) > 3 and (last-1) > int(l[1]):
		print "Overlap at %s:%i (line %i)" % (l[0], int(l[1]), lines)
		exit(1)
	elif len(l) > 3:
		last = int(l[2])
	lines += 1
