Hap.py User's Manual
====================

Introduction
------------

Hap.py is a tool to compare diploid genotypes at haplotype level. Rather than 
comparing VCF records row by row, hap.py will generate and match alternate 
sequences in *haplotype blocks*. A haplotype block is a small region of the 
genome (sized between 1 and around 1000 bp) that contains one or more variants.

Matching haplotype sequences rather than VCF records is more accurate. It allows
 us to the following things:

*  We can match up variant records that represent the same alt sequences in a 
   different form (see [example/GiaB](example/GiaB)).
*  We can also more accurately merge variant call sets 
   (see [ls_example.md](ls_example.md)).

The inputs to hap.py are two VCF files (a "truth" and a "query" file), and an 
optional "confident call region" bed file.

Hap.py will report counts of 

*   ***true-positives (TP)***: variants/genotypes that match in truth and query.
*   ***false-positives (FP)***: variants that have mismatching genotypes or alt 
    alleles, as well as query variant calls in regions a truth set would call 
    confident hom-ref regions.
*   ***false-negatives (FN)*** : variants present in the truth set, but missed 
    in the query.
*   ***non-assessed calls (NA)***: variants outside the truth set regions 

From these counts, we are able to calculate

```
recall = TP/(TP+FN)
precision = TP/(TP+FP)
frac_NA = NA/total(query)
```

These counts and statistics will be calculated for the following subsets of 
variants:

```
|---------------------+-------------------------------------------------------|
|         Type        |                      Description                      |
|---------------------+-------------------------------------------------------|
| Alleles.SNP/INS/DEL | Allele counts. We count occurrences of SNP, insertion |
|                     | or deletion alleles. For hom-alt SNPs, each allele is |
|                     | counted exactly once (so the counts aren't biased     |
|                     | towards hom-alt calls). These counts should be        |
|                     | comparable between different callers.                 |
|---------------------+-------------------------------------------------------|
| Locations.*         | Location counts are useful to compute recall for het, |
|                     | hom, or het-alt calls separately. However,  due to    |
|                     | the way that hap.py performs normalisation, these     |
|                     | counts are not necessarily comparable between         |
|                     | different variant callers (you can read read more     |
|                     | about this in [normalisation.md](normalisation.md))   |
|---------------------+-------------------------------------------------------|
```


Simple Usage
------------

Hap.py is set up to run inside SGE. Unless forced, it will not run interactively
(it detects this by looking for the environment variable SGE_JOB_ID).

Below, we assume that the code has been installed to the directory `${HAPPY}`.

```bash
$ ${HAPPY}/bin/hap.py  \
      example/happy/PG_NA12878_chr21.vcf.gz \
      example/happy/NA12878_chr21.vcf.gz \
      -f example/happy/PG_Conf_chr21.bed.gz \
      -o test \
      --force-interactive
$ ls test.*
test.metrics.json  test.summary.csv
```

This example compares an example run of GATK 1.6 on NA12878 agains the Platinum
Genomes reference dataset (***Note: this is a fairly old version of GATK, so 
don't rely on these particular numbers for competitive comparisons!***).

The summary CSV file contains all computed metrics:

|         Type         | TRUTH.TOTAL | QUERY.TOTAL | METRIC.Recall.HC | METRIC.Precision.HC | METRIC.Frac_NA.HC |
|----------------------|-------------|-------------|------------------|---------------------|-------------------|
| Alleles.DEL          |        6069 |        7020 |         0.907460 |            0.973996 |          0.205698 |
| Alleles.INS          |        6654 |        7179 |         0.880879 |            0.975355 |          0.186098 |
| Alleles.SNP          |       72752 |       67481 |         0.904442 |            0.998361 |          0.023547 |
| Locations.SNP.het    |       32254 |       28665 |         0.873368 |            0.997875 |          0.015175 |
| Locations.SNP.homalt |       20231 |       19270 |         0.929317 |            0.999097 |          0.023560 |


False-positives
---------------

When comparing two VCFs, genotype and allele mismatches are counted as false-
positives. Truth sets like Platinum Genomes or NIST/Genome in a Bottle also
include "confident call regions", which show places where the truth dataset does
not expect variant calls. Hap.py can use these regions to count query variant
calls that do not match truth calls and which fall into  these regions as false
positives.

Full List of Command line Options
---------------------------------

### Minimal Options

You can run hap.py with the -h switch to get help. 

The first two positional arguments are used as the input VCF files. The output 
file prefix is specified using `-o` (this should be something in the form
of directory/prefix):

```
$ ${HAPPY}/bin/hap.py truth.vcf.gz query.vcf.gz \
      -o output-prefix --force-interactive
```

### Running / Debugging issues

```
  --force-interactive   
```

Force running interactively (i.e. when JOB_ID is not in the environment)
Hap.py is set up to run inside SGE. Unless forced, it will not run interactively
(it detects this by looking for the environment variable SGE_JOB_ID). 

The reason for this is that parallelism is implemented using the multiprocessing
module in Python, which spawns processes that can take a lot of memory and also
may be difficult to kill interactively.

```
  --threads THREADS
```

The number of threads to use. This is detected automatically by default using 
Python's multiprocessing module (we recommend around 1GB of RAM per thread).

```
  --logfile LOGFILE     
```

Write logging information into file rather than to stderr.

```
  --scratch-prefix SCRATCH_PREFIX
  --keep-scratch
```

All temporary files go into a scratch folder, which normally defaults to a 
subdirectory of `/tmp`. This can be customised (e.g. when fast local storage is
available).

### Restricting to Subsets of the Genome / Input

```
  --location LOCATIONS, -l LOCATIONS 
```

Add a location to the compare list (when not given, hap.py will use chr1-22,
chrX, chrY).

```
  -P, --include-nonpass
```

Use to include failing variants in comparison. Failing variants are actually 
considered _optional_ rather than mandatory during haplotype comparison (
therefore, )

```
  -R REGIONS_BEDFILE, --restrict-regions REGIONS_BEDFILE
```

Restrict analysis to given (sparse) regions (similar to using -R in bcftools).
Sparse regions should be used when there are not many regions to look at
(the corresponding calls are retrieved via tabix lookup, so this will be
slow for many regions). Also, the regions must not overlap (otherwise, hap.py
will fail).

```
  -T TARGETS_BEDFILE, --target-regions TARGETS_BEDFILE
```

Restrict analysis to given (dense) regions (similar to using -T in bcftools).
One example use for this is to restrict the analysis to exome-only data.

### Additional Input Options

```  
  -f FP_BEDFILE, --false-positives FP_BEDFILE
```

False positive / confident call regions (.bed or .bed.gz).

```
  -r REF, --reference REF 
```

Specify the reference FASTA file to use. Hap.py detects a default reference 
sequence file at the following locations:

*  at `/opt/hap.py-data/hg19.fa` (see the Dockerfile)
*  at the location of the HGREF or the HG19 environment variable

To specify a default reference file location, you can run 

```bash
export HGREF=path-to-your-reference.fa
```

before running hap.py.

### Additional Outputs

```
  -V, --write-vcf 
```

Write an annotated VCF. This file will show the merged and normalised truth
and query calls, together with annotation that shows how they were counted.
Each record contains an INFO field with the following values:

*   `gtt1`/`gtt2` Genotypes in truth and query
*   `type`: FP/TP/FN indicates the call classification
*   `kind`: comparison outcome (missing/GT mismatch/allele mismatch/...);
*   `ctype`: comparison result for the whole block that contains this record
    E.g. `simple:mismatch` would indicate that the block was compared without 
    haplotype matching, and that the variants within it did not match.
*   `HapMatch`: tag that indicates that the surrounding haplotype block 
    matches in truth and query. This tag is only added when haplotype matching
    was run (see [spec.md](spec.md))
*   `Regions`: This is a list of matching regions (currently, only CONF is 
    supported, indicating that a variant call falls into the confident
    call regions)

```
  -B, --write-bed
```

Write a bed file with the haplotype blocks used during comparison. This file
will show the location, comparison method, and comparison outcomes for each
block.

```
  -X, --write-counts
```

Write advanced counts and metrics. This writes additional files / numbers:

*   `output-prefix.counts.json/csv`: raw variant counts. This is produced by 
    the quantify tool, which counts and stratifies variants.
*   `output-prefix.extended.csv`: This file is similar to the summary csv file,
    but it contains more rows and columns showing detailed counts for more 
    allele/location types, as well as TP/FP/FN/UNK counts.

### Input Preprocessing

Hap.py has a range of options to control pre-processing separately for truth 
and query. Most of these require the `--external-preprocessing`  switch to work.

#### BCFtools norm

Truth and query can be preprocessed using `bcftools norm -c x -D` as follows:

```
  --preprocess-truth    Preprocess truth file using bcftools.
  --bcftools-norm       Preprocess query file using bcftools.
```

### Chromosome naming

Normally, the reference fasta file, truth and query should have matching 
chromosome names. However, since e.g. Platinum Genomes doesn't have a specific
GRCH37 version with numeric chromosome names (names like `1` rather than 
`chr1`), this can be worked around in pre-processing (assuming the sequences
are otherwise identical).

```
  --fixchr-truth        Add chr prefix to truth file (default: auto).
  --fixchr-query        Add chr prefix to query file (default: auto).
  --no-fixchr-truth     Disable chr replacement for truth (default: auto).
  --no-fixchr-query     Add chr prefix to query file (default: auto).
  --no-auto-index       Disable automatic index creation for input files. The
                        index is only necessary at this stage if we want to
                        auto-detect locations. When used with -l, and when it
                        is known that there are variants at all given
                        locations this is not needed and can be switched off
                        to save time.
```

Normally this is detected automatically from the VCF headers / tabix indexes.

### Internal Variant Normalisation

```
  --no-internal-preprocessing
```

Switch off xcmp's internal VCF [leftshifting preprocessing](normalisation.md).

### Haplotype Comparison Parameters

```
  -w WINDOW, --window-size WINDOW
```

Window size for haplotype block finder. We use a sliding window, this parameter
determines the maximum distance between two variants that ensures that they
end up in the same haplotype block. For larger values, the haplotype blocks 
get larger and might capture more het variants that will cause the enumeration
threshold (below) to be reached. Also, longer haplotype blocks take more time
to compare. The default value here is 30.

```
  --enumeration-threshold MAX_ENUM
```

Enumeration threshold / maximum number of sequences to enumerate per block. 
Basically, each unphased heterozygous variant in a block doubles the number
of possible alternate sequences that could be created from the set of variants 
within a haplotype block (10 hets in a row would result in 1024 different
alternate sequences).

```
  -e HB_EXPAND, --expand-hapblocks HB_EXPAND
```

Reference-pad and expand the sequences generate haplotype blocks by this many
basepairs left and right.  This is useful for approximate block matching.

