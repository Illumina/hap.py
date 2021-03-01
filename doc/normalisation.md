Variant Normalisation
=====================

This document describes how we normalize and preprocess VCF variant records.

## Hap.py variant pre-processing

The variant pre-processing step in hap.py can be run separately via the pre.py
command line script, which wraps bcftools and hap.py's internal preprocess tool
into a single step:

```bash
${HAPPY_BIN}/pre.py
usage: VCF preprocessor [-h] [--location LOCATIONS] [--pass-only]
                        [--filters-only FILTERS_ONLY] [-R REGIONS_BEDFILE]
                        [-T TARGETS_BEDFILE] [-L] [-D] [--bcftools-norm]
                        [--fixchr] [--no-fixchr] [--bcf] [-v] [-r REF]
                        [-w WINDOW] [--threads THREADS] [--force-interactive]
                        [--logfile LOGFILE] [--verbose | --quiet] [--filter-nonref]
                        [--convert-gvcf-to-vcf]
                        input output
```

The main positional arguments give the names of the input and output files:

```bash
${HAPPY_BIN}/pre.py in.vcf out.vcf.gz
```

This will create a VCF index and make sure there are no major problems with
the input VCF.

We can specify a reference fasta file (the default is an auto-discovered version
of hg19):

```
  -r REF, --reference REF
                        Specify a reference file.
```

We can extract a subset of the VCF with the `-l`, `-R`, `-T` command line options 
as follows:

```
  --location LOCATIONS, -l LOCATIONS
                        Comma-separated list of locations [use naming after
                        preprocessing], when not specified will use whole VCF.
  -R REGIONS_BEDFILE, --restrict-regions REGIONS_BEDFILE
                        Restrict analysis to given (sparse) regions (using -R
                        in bcftools).
  -T TARGETS_BEDFILE, --target-regions TARGETS_BEDFILE
                        Restrict analysis to given (dense) regions (using -T
                        in bcftools).
```

We can also apply filters, or require all variants to pass:

```
  --pass-only           Keep only PASS variants.
  --filters-only FILTERS_ONLY
                        Specify a comma-separated list of filters to apply (by
                        default all filters are ignored / passed on.
```

We can also attempt to canonicalize variant representations by means of left-shifting
and decomposition. These operations are implemented in the preprocess tool which is
part of hap.py. 

```
  -L, --leftshift       Left-shift variants safely (off by default).
  -D, --decompose       Decompose variants into primitives. This results in
                        more granular counts (on by default).
```

BCFtools also implements a normalisation step,  which can be used as an alternative
or together with the switches above. Bcftools will not split variants into primitives,
and also does not attempt safe left-shifting (which means that this pre-processing step
may break the original VCF haplotypes in some circumstances).

```
  --bcftools-norm       Enable preprocessing through bcftools norm -c x -D
```

Normally, the reference fasta file, truth and query should have matching
chromosome names. However, e.g. Platinum Genomes doesn't have a specific
GRCH37 version with numeric chromosome names (names like `1` rather than
`chr1`). This can be worked around in pre-processing (assuming the sequences
are otherwise identical).

Normally this is detected automatically from the VCF headers / tabix indexes /
the reference FASTA file. The reason for having this set of options is that
Platinum Genomes truth VCFs for hg19 (which have chr1 ... chr22 chromosome names)
can be used also on grc37 query VCFs(which have numeric chromosome names) because
PG only provides truth calls on chr1-chr22, chrX (which have identical sequence
in these two references, only chromosome names are different).

Generally, *the truth files (BED and VCF) should have consistent chromosome names
with the reference that is used*, and hap.py can be used to add on a chr prefix
to the query VCF if necessary.


```
  --fixchr              Add/remove chr prefix (default: auto, attempt to match
                        reference).
  --no-fixchr           Add chr prefix to query file (default: auto, attempt
                        to match reference).
```

Using BCF rather than gzipped VCF speeds up file I/O. However, some VCF files
do not translate into BCF easily since BCF needs complete headers that are 
correct for every record in the file. Therefore, this is optional for the time
being. See also the vcfcheck tool which is deployed with hap.py, this tool will
test if a given VCF can be translated into BCF.

```
  --bcf                 Use BCF internally. This is the default when the input
                        file is in BCF format already. Using BCF can speed up
                        temp file access, but may fail for VCF files that have
                        broken headers or records that don't comply with the
                        header.
```

When left-shifting and pre-blocking, we will assume that variants that are further
apart than the following window size will not interact w.r.t. haplotype representations.
The default value is 10000, which should be sufficient for short reads.

```
  -w WINDOW, --window-size WINDOW
                        Preprocessing window size (variants further apart than
                        that size are not expected to interfere).
```

The presence of the <NON_REF> symbolic allele in genome VCFs can cause problems
for hap.py, especially if it is part of a genotype. As a workaround, we 
provide several options. Since variants genotyped as <NON_REF> cannot be
sensibly scored, the we provide the following option, which is safe to use
on both genome VCFs and standard VCFs:

```
  --filter-nonref       Remove any variants genotyped as <NON_REF>.                 
```

If hap.py still crashes when processing a genome VCF, we provide separate
options to perform on-the-fly conversion of a genome VCF to a standard VCF
by removing all <NON_REF> alleles and non-variant blocks. Note that this 
also removes some fields from the INFO column. These options should only
be used on genome VCFs since attempting to convert a standard VCF will 
cause all biallelic variants to be filtered out (most of them).

```
  --convert-gvcf-to-vcf Convert the input genome VCF to a standard VCF.
```

Runtime behaviour can also be controlled as follows:

```
  --threads THREADS     Number of threads to use.
  --force-interactive   Force running interactively (i.e. when JOB_ID is not
                        in the environment)
  --logfile LOGFILE     Write logging information into file rather than to
                        stderr
  --verbose             Raise logging level from warning to info.
  --quiet               Set logging level to output errors only.
```

## Normalisation details

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

