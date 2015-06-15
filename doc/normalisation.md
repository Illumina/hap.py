Variant Normalisation
=====================

This document describes how we normalize and preprocess VCF variant records.

A simple way to debug normalisation is via the multimerge tool that is  compiled
as part of hap.py. This is a customised variant of the bcftools merge command
which includes a set of normalisation steps.

Here is how to run multimerge with the full set of operations that will also be
performed by hap.py (assuming that the environment variable `HG` points to a
reference sequence that matches the input VCFs).

```bash
multimerge input.vcf.gz  -r $HG -o output.vcf.gz --process-full=1
```

Multimerge also accepts multiple input VCF files (e.g. a truth and a query file)
and the results might differ due to the nature of the normalisation process
described below. It is possible to only use specific steps of this process, run
`multimerge -h` for a list of command line switches.

## Step 0 (optional): Realign and split into allelic primitives

This operation will realign complex alleles against the reference sequence and
output primitive alleles. It is part of hap.py's "partial credit" mode, in which
we are able to match complex alleles partially.

## Step 1: Remove unused alleles

This removes alleles not seen in any genotype:

```
chr1    100000  .   A    C,T  ...   GT  0/1
```

... becomes

```
chr1    100000  .   A    C  ...   GT  0/1
```

## Step 2: Split het-alts

This splits het-alt calls into separate rows.

```
chr1    100000  .   A    C,T  ...   GT  1/2
```

... becomes

```
chr1    100000  .   A    C  ...   GT  0/1
chr1    100000  .   A    T  ...   GT  0/1
```

## Step 3: Left-shift Indel Alleles

We re-implemented the algorithm from here: [http://genome.sph.umich.edu/wiki/Variant_Normalization](http://genome.sph.umich.edu/wiki/Variant_Normalization).

## Step 4: Aggregate Calls at the same location

There are multiple levels for doing this (`multimerge --merge-by-location X`
where `X` is 0, 1, or 2). Level 2 is used inside hap.py.

Level 0 disables this step.

Level 1 merges calls across samples:

```
chr1    100000  .   A    C  ...   GT  ./.  0/1
chr1    100000  .   A    T  ...   GT  0/1  ./.
```

... becomes

```
chr1    100000  .   A    C,T  ...   GT  0/2  0/1
```

Level 2 combines het calls to het-alt:

```
chr1    100000  .   A    C  ...   GT  0/1
chr1    100000  .   A    T  ...   GT  0/1
```

... becomes

```
chr1    100000  .   A    C,T  ...   GT  1/2
```

## Step 5: Identify matching alleles

```
chr1    100000  .   A    C,T,T  ...   GT  1/2  1/3
```

... becomes

```
chr1    100000  .   A    C,T  ...   GT  1/2  1/2
```

