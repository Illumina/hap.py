#!/usr/bin/env python

# Compare two som.py csv files

import sys


def csvread(filename):
    f = open(filename)

    header = None
    rows = []
    for l in f:
        l = l.replace("\n", "")
        if not header:
            header = l.split(",")
        else:
            rows.append(l.split(","))

    data = {}
    for r in rows:
        label = r[1]
        for i in xrange(2, len(r)):
            if header[i] not in data:
                data[header[i]] = {}
            try:
                data[header[i]][label] = float(r[i])
            except:
                data[header[i]][label] = r[i]
    f.close()
    return data


def main():
    data1 = csvread(sys.argv[1])
    data2 = csvread(sys.argv[2])

    for metric in ['fp', 'ambiguous', 'recall', 'recall_lower', 'recall_upper', 'precision',
                   'precision_lower', 'precision_upper', 'tp', 'total.query', 'ambi', 'na',
                   'recall2', 'unk', 'total.truth', 'fn', 'fp.region.size', 'fp.rate']:
        if len(data1[metric]) == 0 or len(data2[metric]) == 0:
            raise Exception("Number of metrics to compare is wrong")

        for field in data1[metric].iterkeys():
            if field == "no-ALTs":
                continue
            try:
                if ("%8.3f" % data1[metric][field]) != ("%8.3f" % data2[metric][field]):
                    raise Exception("Failed: Results should be the same")
                print metric + " / " + field + "... ok (%8.3f)" % data2[metric][field]
            except:
                print metric + " / " + field + "... ERROR (%s / %s)" % (data1[metric][field], data2[metric][field])
                raise


if __name__ == '__main__':
    main()

