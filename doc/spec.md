Haplotype Comparison Tools
==========================

Introduction
------------

The core idea in haplotype comparison is to enumerate all possible haplotypes
that are specified by a given VCF.

Here is an example:

```
>chrQ
ACCACCACC

chrQ  4  A  T   0|1
chrQ  5  C  T   1|0

=> two haplotypes

ACCTCCACC
ACCATCACC
```

When variants are not phased, we don't know exactly which haplotype they are on.
We can then enumerate all possible haplotype pairs:

```
chrQ  4  A  T   0/1
chrQ  5  C  T   0/1

=> pair 1 (variants are on the same haplotype)
ACCACCACC | ACCTTCACC

=> pair 2 (variants are on alternate haplotypes)
ACCTCCACC | ACCATCACC
```

Now assume we are given two samples and we would like to know if they have the
same haplotypes:

```
CHROM POS  REF  ALT  SAMPLE1   SAMPLE2
chrQ  4    AC   TT   ./.       0|1
chrQ  4    A    T    0/1       ./.
chrQ  5    C    T    0/1       ./.
```

We know from above that `SAMPLE1` has two possible haplotype pairs:

1.  `ACCACCACC | ACCTTCACC` ... or
2.  `ACCTCCACC | ACCATCACC`

In `SAMPLE2`, we have exactly one pair: `ACCACCACC | ACCTTCACC`, which matches
pair 1. from above.

This type of comparison can be applied in small windows across the genome to
identify exactly matching haplotypes that are represented differently in VCF
records.

Enumeration of all possible sequences is implemented using
[reference graphs](refgraph.md).

Tools / Python wrappers
-----------------------

The main tools of this package wrap the building blocks in the next section.
Each tool will show its command line options when run with the `-h` switch.

### hap.py

This tool performs multi-threaded haplotype comparison. See
[happy.md](happy.md).

### som.py

This is a wrapper around bcftools to compare VCFs only by alleles (useful for
comparing/benchmarking somatic calls). See [sompy.md](sompy.md).

### cnx.py

This tool tries to identify variant callers and aligners by analysing the VCF
and BAM headers. The output is a JSON file.

### ftx.py

This is a VCF feature extraction script that can extract feature tables from
VCF files.

Building blocks
---------------

This section outlines all the tools that come with this package. All tools
typically require access to the reference sequence in indexed fasta format.

To create an index for your reference file, run
`samtools faidx reference.fasta`.

Also, all the tools below require input VCFs to be bgzipped and indexed,
otherwise they will fail.

### Enumerate haplotypes: `hapenum`

Input: a VCF file, a (small) region in the genome.

Outputs:

*  all possible Haplotype sequences described by the VCF (e.g. exactly two
   for phased diploid VCF files).
*  a dot file showing the [reference graph](refgraph.md)


### Enumerate haplotype pairs: `dipenum`

Input: a VCF file, a (small) region in the genome.

Output:

All feasible haplotype pairs for this region.

### Preprocess and merge: `multimerge`

Input: one or more VCF files, optionally regions in the genome.

Output:

A multi-sample VCF file that has been [preprocessed](normalisation.md).

### Split a VCF but don't break haplotype block boundaries: `blocksplit`

Input: one or more VCF files, optionally regions in the genome.

Output:

A bed file with blocks that guarantee that block boundaries do not fall between
two variants that are closer than a window length *w*.

### Turn a VCF header into JSON format: `vcfhdr2json`

Input: a VCF file

Output:

The VCF header and tabix contigs in JSON format.

### Compare two VCFs: `xcmp`

This is the core comparison engine in hap.py.

Input: two VCF files, optionally regions of the genome

Outputs:

*   An annotated VCF file showing match / mismatch information for each location

### Count variants: `quantify`

This is a helper to debug haplotype comparison in fixed blocks.

Input: a multi-sample VCF file, optionally bed annotation regions for
       stratification

Outputs:

*  a JSON file of stratified counts
*  a tab-separated table with all ROC values
*  a VCF file with annotated regions

### Make ROC tables: `roc`

Input: tab-separated table of instances classified as TP/FP/FN with a quality
       value each

Output: table of precision / recall values with varying quality threshold.

### Split a number of VCF files into blocks/superloci: `blocksplit`

Given number of VCF files, superloci contain all variant calls that
are no further apart than a given window length.

Input: one or more VCF files

Output: a bed file with all superloci

### Check a VCF file for errors: `vcfcheck`

VCF files with invalid / incomplete headers can cause problems when converting them
into the BCF format. Also, in some cases, certain types of alleles or reference
overlaps might cause problems when comparing variants.

Input: a single VCF / BCF file

Output: return code != zero if the VCF file has problems + some statistics to stdout