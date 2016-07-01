# Micro-benchmark

To test hap.py and show how ROC curves work, we have a small set of test
datasets in this repository. These are contained in the [example/happy]()
folder.

The folder contains partial callsets for Platinum Genomes 2016.1 and for
GATK3, Platypus and Freebayes on NA12878. All three methods were run using
their joint-calling mode across the whole PG pedigree. We benchmark
their precision and recall against the PG callset using different methods
and produce ROCs. [example/happy/microbenchmark.sh]() contains a script to
run hap.py in different configurations on these datasets. It needs to be
run either from a hap.py build folder, or pointed to a hap.py installation
using the `HCDIR` environment variable.

The script will produce benchmarking datasets for the following configurations:

{ GATK, Platypus, Freebayes } x { xcmp, vcfeval } x { decompose, no-decompose }

Xcmp and vcfeval are two different comparison engines, and decompose /
no-decompose are two options for pre-processing. Decomposing query variants
into primitives allows for more granular counting of TPs and FPs.

A first observat

| ![microbench_GATK.SNP.png]() | ![microbench_GATK.SNP.png]() |
|------------------------------|------------------------------|