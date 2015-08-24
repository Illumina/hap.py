#!/usr/bin/env python

# Compare two som.py csv files

import sys
import csv


def csvread(filename):
    f = open(filename)
    header = None
    data = []
    rn = 0
    for r in csv.reader(f):
        rn += 1
        if not header:
            header = r
            continue
        rec = {}
        for i in xrange(0, len(r)):
            try:
                try:
                    rec[header[i]] = float(r[i])
                except:
                    rec[header[i]] = r[i]
            except:
                print >> sys.stderr, "%i: %i / %i" % (rn, len(header), len(r))
                print "%i: %i\n%s\n%s" % (rn, i, "\t".join(header), "\t".join(r))
                raise

        data.append(rec)
    f.close()
    return header, data


def main():
    header1, data1 = csvread(sys.argv[1])
    header2, data2 = csvread(sys.argv[2])

    if header1 != header2 or len(header1) == 0:
        raise Exception("Header mismatch, \n%s\n != \n%s\n" % (str(header1), str(header2)))

    if len(data1) != len(data2) or len(data1) == 0:
        raise Exception("Data length mismatch, %i != %i" % (len(data1), len(data2)))

    print "Comparing fields: %s" % str(header1)

    for i in xrange(0, len(data1)):
        rec1 = data1[i]
        rec2 = data2[i]
        for metric in header1:
            if type(rec1[metric]) is str and type(rec2[metric]) is str:
                if rec1[metric] != rec2[metric]:
                    raise Exception("Str mismatch in row %i: %s != %s" % (i, rec1[metric], rec2[metric]))
            elif type(rec1[metric]) is float and type(rec2[metric]) is float:
                if ("%.3g" % rec1[metric]) != ("%.3g" % rec2[metric]):
                    raise Exception("Float mismatch in row %i: %g != %g" % (i, rec1[metric], rec2[metric]))
            else:
                raise Exception("Type mismatch in row %i: %s != %s" % (i, str(rec1[metric]), str(rec2[metric])))
    print "ok (%i records)" % len(data1)


if __name__ == '__main__':
    main()
