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
import sys
import pandas
import time
import math

import Tools

def dataframeToMetricsTable(table_id, df):
    """ Convert a pandas dataframe to a PUMA metrics table

    This is how a metrics table looks like:

     {
         id:  'readCounts',
         label: 'Number of Read Pairs',
         type: 'Table',
         data:  [ { id : 'x', label: 'Insert size', type: 'int32', values: [ 2, 3, 4, 5, 6, 7, 8, 9] },
                  { id : 'y', label: 'Number of read pairs', type: 'int32', values: [ 88, 24, 14, 4, 2, 3, 7, 8]}
         ],
         properties: [ {key: 'represents', value : 'histogram'},
                       {key: 'xlim', value : [0, null]}
         ]
     }

    :param df: the dataframe
    :type df: pandas.DataFrame
    :return: a dictionary in metrics table format
    :rtype: dict
    """
    mdict = {"id": table_id, "label": table_id, "type": "Table", "data": [], "properties": []}


    ldict = {'id': 'types', "label": 'types', "type": "string", "values": list(df.index.values)}
    mdict["data"].append(ldict)
    for header in list(df):
        c = str(df[header].dtype)
        ccast = str
        if c in ["int32", "int64"]:
            coltype = c
            ccast = int
        elif c.startswith("float"):
            coltype = "double"
            ccast = float
        else:
            coltype = "string"
        ldict = {'id': header, "label": header, "type": coltype, "values": list(map(ccast, df[header].values))}
        mdict["data"].append(ldict)

    return replaceNaNs(mdict)


def makeMetricsObject(name):
    """ Create PUMA metrics dictionary

    This is how they look like
     {
      name: 'PAM.GCbias',
      timestamp: 'Wed Sep 05 16:23:46 2012',
      version: '1.0',
      runInfo: [ {key:  'version-of-PUMA', value: '2.1.4'} ].
      sampleInfo: [ {key: 'sample', value: 'human1'},
                    {key: 'organism', value: 'hg19'} ],
      parameters: [ {key: 'refFile', value: '/path/to/ref/file'},
                    {key: 'window', value: 10} ],
      metrics: [ { id : 'iSize', label: 'Insert size histogram', type: 'Table', data: [], properties: {} },
                 { id : 'iSizeBAM', label: 'Insert size from BAM', type: 'Table', data: [] },
                 { id : 'readPairs', label: 'Numbers of read pairs', type: 'Map', data: [] },
                 { id : 'GCb', label: 'GC related metric', type: 'Table', data: [], properties: {} }
     ]
     }
    """

    version = "%s" % Tools.version

    mdict = {'name': name,
             'timestamp': time.strftime("%a %b %d %X %Y"),
             'version': version,
             'runInfo': [{"key": "commandline", "value": " ".join(sys.argv)}],
             'metadata': {
                 "required": {
                     "id": "haplotypes",
                     'version': version,
                     "module": "%s" % os.path.basename(sys.argv[0]),
                     "description": "%s generated this JSON file via command line %s" % (
                         sys.argv[0], " ".join(sys.argv))
                 }
             },
             'sampleInfo': [],
             'parameters': [],
             'metrics': []}

    return mdict


def replaceNaNs(xobject):
    """ Replace all NaNs in the given object with string values
    :param xobject: a dictionary
    :return: object, fixed
    """

    if type(xobject) is dict:
        for k in xobject.keys():
            if type(xobject[k]) is dict or type(xobject[k]) is list or type(xobject[k]) is float:
                xobject[k] = replaceNaNs(xobject[k])
    elif type(xobject) is list:
        for k in xrange(0, len(xobject)):
            if type(xobject[k]) is dict or type(xobject[k]) is list or type(xobject[k]) is float:
                xobject[k] = replaceNaNs(xobject[k])
    elif type(xobject) is float:
        # NaN and Inf become null. Not elegant, a JavaScript-y solution would be
        # to put these in as strings
        if math.isnan(xobject) or math.isinf(xobject):
            xobject = None
    return xobject
