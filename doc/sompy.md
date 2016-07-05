Somatic Variant Benchmarking Tool
=================================

## Overview

This tool compares variants by location and alleles (using bcftools isec).

The input are a truth and a query file. Optionally, false-positive regions
can be specified.

Counting works as follows:

*  Variants that are in both files count as true positives (tp). 
*  Variants that are only in the query become false positives or unknowns
    -  if FP regions are specified, query-only files in these regions become fp,
       and other calls become unknowns (unk)
    -  if no FP regions are specified, all query-only variants become FPs
*  Variants only in the truth file become false-negatives.

We compute the following metrics:

*  Relative recall: `recall = tp/(tp + fn)`
*  Absolute recall: `recall2 = tp/total(truth)`
*  Precision: `precision = tp/(tp + fp)`
*  Not assessed: `na = unk/total(query)`

Confidence intervals for relative recall and precision are calculated based on
the Jeffreys method for binomial proportions. The lower and upper bound for both
metrics is reported as `recall_lower`, `recall_upper` and `precision_lower`,
`precision_upper`, respectively.

```
${HAPPY}/bin/som.py example/sompy/PG_admix_truth_snvs.vcf.gz \ 
                    example/sompy/strelka_admix_snvs.vcf.gz \
                    -f example/sompy/FP_admix.bed.gz \
                    -o test
[...]

      type  total.truth  total.query     tp   fp    fn  unk  ambi    recall  recall_lower  recall_upper   recall2  precision  precision_lower  precision_upper   na  ambiguous  fp.region.size    fp.rate
1     SNVs        16235        15934  15131  803  1104    0     0  0.931999      0.928049      0.935794  0.931999   0.949605         0.946126         0.952920  0.0        0.0        32984914  24.344462
5  records        16235        15936  15131  805  1104    0     0  0.931999      0.928049      0.935794  0.931999   0.949485         0.946003         0.952804  0.0        0.0        32984914  24.405096

ls test.*
test.stats.csv
```

It is possible to stratify the unknown / FP variants using *ambiguity regions*.
These regions can be given as a bed file with at least four columns, e.g.:

```
chr1    0   10000   *   all_n
```

This would define an ambiguity region named all_n. All variants falling into 
these regions will be counted as ambi in the table above. The command line
option `--explain_ambiguous` will show how many calls were captured by each 
region type.

## Feature Tables

A feature table is a CSV file which contains all variant locations in truth and query together with
their status in a "tag" column (TP/FP/FN/UNK/AMBI). We implement feature tables for some common
variant callers and truth set types.

For more information, see [/src/python/Somatic/\_\_init\_\_.py](/src/python/Somatic/\_\_init\_\_.py).

## ROC Curves

Som.py can create ROC curve data based on feature tables. This is implemented using 
the `--roc` switch. Different variant callers produce different features, and therefore
might require different feature tables to be computed -- the `--roc` switch therefore
overrides the setting in `--feature-table`.

## Stratification by allele frequency

If the feature table contains columns for an expected and an observed allele frequency,
we can stratify precision and recall into AF bins. This feature requires a feature 
table to be selected using `--feature-table` or `--roc`. The columns from which to obtain
the AF values for TP/FN can be selected using `--af-truth`, and the columns for FP/UNK/AMBI,
are selected using `--af-query`.

## All Options

```
usage: Somatic Comparison [-h] -o OUTPUT [-l LOCATION] [-R REGIONS_BEDFILE]
                          [-T TARGETS_BEDFILE] [-f FP] [-a AMBI] [--ambi-fp]
                          [--no-ambi-fp] [--count-unk] [--no-count-unk] [-e]
                          [-r REF] [--scratch-prefix SCRATCH_PREFIX]
                          [--keep-scratch] [--continue] [-P]
                          [--feature-table {
                              hcc.strelka.snv,
                              hcc.mutect.snv,
                              hcc.varscan2.indel,
                              admix.strelka.snv,
                              generic,
                              hcc.strelka.indel,
                              admix.strelka.indel,
                              hcc.varscan2.snv,
                              hcc.mutect.indel}]
                          [--bam BAMS] [--normalize-truth] [--normalize-query]
                          [-N] [--fixchr-truth] [--fixchr-query]
                          [--fix-chr-truth] [--fix-chr-query]
                          [--no-fixchr-truth] [--no-fixchr-query]
                          [--no-order-check]
                          [--roc {
                              varscan2.indel,
                              strelka.snv.qss,
                              mutect.snv,
                              strelka.snv,
                              strelka.indel.evs,
                              varscan2.snv,
                              mutect.indel,
                              strelka.indel,
                              strelka.snv.vqsr}]
                          [--bin-afs] [--af-binsize AF_STRAT_BINSIZE]
                          [--af-truth AF_STRAT_TRUTH]
                          [--af-query AF_STRAT_QUERY] [-FN]
                          [--fp-region-size FPR_SIZE] [--ci-level CI_LEVEL]
                          [--logfile LOGFILE] [--verbose | --quiet]
                          truth query

positional arguments:
  truth                 Truth VCF file
  query                 Query VCF file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file prefix for statistics and feature table
                        (when selected)
  -l LOCATION, --location LOCATION
                        Location for bcftools view (e.g. chr1)
  -R REGIONS_BEDFILE, --restrict-regions REGIONS_BEDFILE
                        Restrict analysis to given (sparse) regions (using -R
                        in bcftools).
  -T TARGETS_BEDFILE, --target-regions TARGETS_BEDFILE
                        Restrict analysis to given (dense) regions (using -T
                        in bcftools).
  -f FP, --false-positives FP
                        False-positive region bed file to distinguish UNK from
                        FP
  -a AMBI, --ambiguous AMBI
                        Ambiguous region bed file(s) to distinguish from FP
                        (e.g. variant only observed in some replicates)
  --ambi-fp             Use FP calls from ambiguous region files also.
  --no-ambi-fp          Do not use FP calls from ambiguous region files also.
  --count-unk           Assume the truth set covers the whole genome and only
                        count FPs in regions specified by the truth VCF or
                        ambiguous/false-positive regions.
  --no-count-unk        Do not use FP calls from ambiguous region files also.
  -e, --explain_ambiguous
                        print a table giving the number of ambiguous events
                        per category
  -r REF, --reference REF
                        Specify a reference file.
  --scratch-prefix SCRATCH_PREFIX
                        Filename prefix for scratch report output.
  --keep-scratch        Filename prefix for scratch report output.
  --continue            Continue from scratch space (i.e. use VCFs in there if
                        they already exist).
  -P, --include-nonpass
                        Use to include failing variants in comparison.
  --feature-table {hcc.strelka.snv,hcc.mutect.snv,hcc.varscan2.indel,
      admix.strelka.snv,generic,hcc.strelka.indel,admix.strelka.indel,
      hcc.varscan2.snv,hcc.mutect.indel}
                        Select a feature table to output.
  --bam BAMS            pass one or more BAM files for feature table
                        extraction
  --normalize-truth     Enable running of bcftools norm on the truth file.
  --normalize-query     Enable running of bcftools norm on the query file.
  -N, --normalize-all   Enable running of bcftools norm on both truth and
                        query file.
  --fixchr-truth        Add chr prefix to truth file (default: true).
  --fixchr-query        Add chr prefix to query file (default: true).
  --fix-chr-truth       Same as --fixchr-truth.
  --fix-chr-query       Same as --fixchr-query.
  --no-fixchr-truth     Disable chr replacement for truth (default: false).
  --no-fixchr-query     Add chr prefix to query file (default: false).
  --no-order-check      Disable checking the order of TP features (dev
                        feature).
  --roc {varscan2.indel,strelka.snv.qss,mutect.snv,strelka.snv,strelka.indel.evs,
      varscan2.snv,mutect.indel,strelka.indel,strelka.snv.vqsr}
                        Create a ROC-style table. This is caller specific -
                        this will override the --feature-table switch!
  --bin-afs             Stratify into different AF buckets. This needs to have
                        features availablefor getting the AF both in truth and
                        query variants.
  --af-binsize AF_STRAT_BINSIZE
                        Bin size for AF binning (should be < 1). Multiple bin
                        sizes can be specified using a comma, e.g.
                        0.1,0.2,0.5,0.2 will split at 0.1, 0.3, 0.8 and 1.0.
  --af-truth AF_STRAT_TRUTH
                        Feature name to use for retrieving AF for truth
                        variants (TP and FN)
  --af-query AF_STRAT_QUERY
                        Feature name to use for retrieving AF for query
                        variants (FP/UNK/AMBI)
  -FN, --count-filtered-fn
                        Count filtered vs. absent FN numbers. This requires
                        the -P switch (to use all variants) and either the
                        --feature-table or --roc switch.
  --fp-region-size FPR_SIZE
                        How to obtain the normalisation constant for FP rate.
                        By default, this will use the FP region bed file size
                        when using --count-unk and the size of all reference
                        contigs that overlap with the location specified in -l
                        otherwise. This can be overridden with: 1) a number of
                        nucleotides, or 2) "auto" to use the lengths of all
                        contigs that have calls. The resulting value is used
                        as fp.region.size.
  --ci-level CI_LEVEL   Confidence level for precision/recall confidence
                        intervals (default: 0.95)
  --logfile LOGFILE     Write logging information into file rather than to
                        stderr
  --verbose             Raise logging level from warning to info.
  --quiet               Set logging level to output errors only.
```
