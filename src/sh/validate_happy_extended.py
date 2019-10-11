#!/usr/bin/env python

# Compare som.py stats.csv and hap.py extended.csv files

import sys
import argparse
import csv
import pprint as pp
import re
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser("Validate hap.py stats")
    parser.add_argument("--sompy_stats", required=True, help="Path to som.py stats.csv")
    parser.add_argument("--happy_extended", required=True, help="Path to hap.py extended.csv")
    args = parser.parse_args()
    return args

def eval_equal(metric_name, count_a, count_b):
    a = int(count_a)
    b = int(count_b)
    e = "PASS" if a == b else "FAIL"
    logging.info("%s: %d vs %d - %s" % (metric_name, a, b, e))
    return e

def parse_sompy_stats(path):
    sompy_stats = csv.DictReader(open(path))
    result = dict()
    for s in sompy_stats:
        subset = s["type"]
        if re.match("indels.", subset):
            # store results per af_bin
            m = re.findall("indels.(\d+\.\d+)-(\d+\.\d+)", subset)[0]
            af_low = m[0][:4]
            af_high = m[1][:4]
            af_bin = "[%s,%s]" % (af_low, af_high) if af_low == "1.00" else "[%s,%s)" % (af_low, af_high)
            result["INDEL." + af_bin] = s

        if re.match("SNVs.", subset):
            # store results per af_bin
            m = re.findall("SNVs.(\d+\.\d+)-(\d+\.\d+)", subset)[0]
            af_low = m[0][:4]
            af_high = m[1][:4]
            af_bin = "[%s,%s]" % (af_low, af_high) if af_low == "1.00" else "[%s,%s)" % (af_low, af_high)
            result["SNP." + af_bin] = s

    return result

if __name__ == '__main__':
    args = parse_args()

    sompy_stats = parse_sompy_stats(path=args.sompy_stats)
    happy_extended = csv.DictReader(open(args.happy_extended))

    outcomes = dict(ALL=set(), PASS=set())
    for h in happy_extended:
        k = h["Type"] + "." + h["Subset"]
        try:
            s = sompy_stats[k]
        except KeyError:
            s = {"total.truth": 0, "tp": 0, "fn": 0, "total.query": 0, "fp": 0, "unk": 0}

        outcomes[h["Filter"]].add(eval_equal(metric_name="%s %s TRUTH.TOTAL" % (k, h["Filter"]), count_a=s["total.truth"], count_b=h["TRUTH.TOTAL"]))
        outcomes[h["Filter"]].add(eval_equal(metric_name="%s %s TRUTH.TP" % (k, h["Filter"]), count_a=s["tp"], count_b=h["TRUTH.TP"]))
        outcomes[h["Filter"]].add(eval_equal(metric_name="%s %s TRUTH.FN" % (k, h["Filter"]), count_a=s["fn"], count_b=h["TRUTH.FN"]))
        outcomes[h["Filter"]].add(eval_equal(metric_name="%s %s QUERY.TOTAL" % (k, h["Filter"]), count_a=s["total.query"], count_b=h["QUERY.TOTAL"]))
        outcomes[h["Filter"]].add(eval_equal(metric_name="%s %s QUERY.FP" % (k, h["Filter"]), count_a=s["fp"], count_b=h["QUERY.FP"]))
        outcomes[h["Filter"]].add(eval_equal(metric_name="%s %s QUERY.UNK" % (k, h["Filter"]), count_a=int(s["unk"])+int(s["ambi"]), count_b=h["QUERY.UNK"]))

    failed_vfilters = [x for x in outcomes if "FAIL" in outcomes[x]]
    if len(failed_vfilters) == 2:
        logging.info("Failed filters: %s" % failed_vfilters)
        sys.exit(1)
    else:
        logging.info("DONE")
