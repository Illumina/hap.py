#!/usr/bin/env python

# Compare two summary csv files

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
        label = r[0]
        if "hap.py" in label:
            # ignore version line
            continue
        for i in xrange(1, len(r)):
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

    for metric in ["TRUTH.TOTAL", "QUERY.TOTAL",
                   "TRUTH.TOTAL.TiTv_ratio", "QUERY.TOTAL.TiTv_ratio",
                   "TRUTH.TOTAL.het_hom_ratio", "QUERY.TOTAL.het_hom_ratio",
                   "METRIC.Recall", "METRIC.Precision", "METRIC.Frac_NA"]:
        for field in ["Locations.SNP", "Locations.SNP.het",
                      "Locations.INDEL", "Locations.INDEL.het"]:
            print metric + " / " + field
            print data1[metric][field]
            print data2[metric][field]
            if metric.endswith("_ratio") and data1[metric][field] == "" and data2[metric][field] == "":
                # allow empty ratio match
                continue
            if ("%.3g" % data1[metric][field]) != ("%.3g" % data2[metric][field]):
                raise Exception("Failed: Results should be the same")

if __name__ == '__main__':
    main()
