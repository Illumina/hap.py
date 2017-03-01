#!/usr/bin/env python

# Compare two extended csv files

import sys
import csv

def csvread(filename):
    f = open(filename)
    freader = csv.DictReader(f)
    data = {}
    label_columns = ["Type", 
                     "Subtype", 
                     "Subset",
                     "Filter", 
                     "Genotype",
                     "QQ.Field",
                     "QQ"]

    for l in freader:
        record = dict(l)
        label = "_".join([record[k] for k in label_columns])
        if "hap.py" in label:
            # ignore version line
            continue

        for k in record.keys():
            try:
                record[k] = int(record[k])
            except:
                try:
                    record[k] = float(record[k])
                except:
                    if record[k] == "." or record[k] == "":
                        record[k] = float("NaN")
        data[label] = record

    return data


def main():
    data1 = csvread(sys.argv[1])
    data2 = csvread(sys.argv[2])

    all_keys_1 = set()
    for l in data1.values():
        all_keys_1 |= set(l.keys())
    all_keys_2 = set()
    for l in data2.values():
        all_keys_2 |= set(l.keys())

    if len(all_keys_1 - all_keys_2) or len(all_keys_2 - all_keys_1):
        print >> sys.stderr, "ERROR -- missing metrics:"
        print >> sys.stderr, "-------------------------\n"
        print >> sys.stderr, "In 1 but not 2: %s" % str(all_keys_1 - all_keys_2)
        print >> sys.stderr, "In 2 but not 1: %s" % str(all_keys_1 - all_keys_2)
        print >> sys.stderr, "-------------------------"
        sys.exit(1)

    all_labels_1 = set(data1.keys())
    all_labels_2 = set(data2.keys())

    if len(all_labels_1 - all_labels_2) or len(all_labels_2 - all_labels_1):
        print >> sys.stderr, "ERROR -- missing rows:"
        print >> sys.stderr, "-------------------------\n"
        print >> sys.stderr, "In 1 but not 2: %s" % str(all_labels_1 - all_labels_2)
        print >> sys.stderr, "In 2 but not 1: %s" % str(all_labels_1 - all_labels_2)
        print >> sys.stderr, "-------------------------"
        sys.exit(1)

    print "Comparing %i labels and %i metrics..." % (len(all_labels_1), len(all_keys_1))

    different_metrics = []

    for label in all_labels_1:
        for metric in all_keys_1:
            v1 = data1[label][metric]
            v2 = data2[label][metric]

            if metric.endswith("_ratio"):
                if v1 == ".":
                    v1 = ""
                if v2 == ".":
                    v2 = ""

            if type(v1) != type(v2):
                different_metrics.append((label,
                                          metric,
                                          ("%s" % type(v1)),
                                          ("%s" % type(v2)),
                                          "type mismatch"))
            elif type(v1) is float:
                if ("%.3g" % v1) != ("%.3g" % v2):
                    different_metrics.append((label,
                                              metric,
                                              ("%.3g" % v1),
                                              ("%.3g" % v2),
                                              str(v1 - v2)))
            elif type(v1) is int:
                if v1 != v2:
                    different_metrics.append((label,
                                              metric,
                                              str(v1),
                                              str(v2),
                                              str(v1 - v2)))
            else:
                if v1 != v2:
                    different_metrics.append((label,
                                              metric,
                                              str(v1),
                                              str(v2),
                                              "Non-numeric mismatch"))

    if different_metrics:
        print >> sys.stderr, "ERROR -- Metric differences detected:"
        print >> sys.stderr, "-------------------------------------\n"
        for m in different_metrics:
            print >> sys.stderr, "%s / %s: %s != %s difference: %s" % m
        print >> sys.stderr, "-------------------------------------"
        sys.exit(1)


if __name__ == '__main__':
    main()
