Haplotype Comparison Tools
==========================

Peter Krusche <pkrusche@illumina.com>

This is a set of programs based on [htslib](https://github.com/samtools/htslib)
to compare VCF files by specified haplotype.

Rather than comparing entries individually, we produce a graph-based reference
from the VCF entries, create all possible haplotype sequences, and compare
these by alignment / exact matching. This is more accurate in cases like this:

*Variant representation 1 (shown in purple in the image below):*

```
CHROM POS   REF  ALT             GT
chrQ  10    G    GTGTGTGCATGCT   0/1
```

*Variant representation 2 (shown in green in the image below):*

```
CHROM POS   REF  ALT             GT
chrQ  16    G    GCATGCT         0/1
chrQ  19    T    TGTGTG          0/1
```

![](doc/rep_ex.PNG)

Both representations in this example are able to produce the same alt sequences,
but we are not able to match them up with standard VCF tools. In particular,
we can see from this example that the second representation actually allows us
to create two different sets of alt sequences (because we don't know the phasing
of our variants -- the insertions could happen on different haplotypes with
representation 2).

With this tool, we can produce all haplotypes sequences by enumerating paths
through a reference graph. By finding the paths / alt alleles that are
consistent between two VCFs files we can produce accurate benchmarking
numbers for comparing a VCF to a gold standard truth set.

See [doc/spec.md](doc/spec.md) for more information.

Simple Usage
============

The main two tools are hap.py (diploid precision/recall evaluation) and som.py
(somatic precision/recall evaluation -- this ignores the GT and just checks for
presence of alleles).

Here are some small example command lines. Advanced features like confident call
 / ambiguity / FP regions are also available, see the documentation for each
 tool for these.

Below, we assume that the code has been installed to the directory `${HAPPY}`.

### hap.py

See also [doc/happy.md](doc/happy.md).

```bash
$ ${HAPPY}/bin/hap.py  \
      example/happy/PG_NA12878_chr21.vcf.gz \
      example/happy/NA12878_chr21.vcf.gz \
      -f example/happy/PG_Conf_chr21.bed.gz \
      -o test
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

These numbers tell us the SNP and indel recall of our query VCF against the
truth dataset. See [doc/happy.md](doc/happy.md) for more documentation and some
advice for their interpretation.

### som.py

See [doc/sompy.md](doc/sompy.md) for more documentation.

```
${HAPPY}/bin/som.py example/sompy/PG_admix_truth_snvs.vcf.gz \
                    example/sompy/strelka_admix_snvs.vcf.gz \
                    -f example/sompy/FP_admix.bed.gz \
                    -o test
[...]
      type  total.truth  total.query     tp     fp   fn    unk  ambi    recall   recall2  precision        na  ambiguous
1     SNVs        16235        47530  15573  14698  662  17259     0  0.959224  0.959224   0.514453  0.363118          0
3  records        16235        47779  15573  14737  662  17469     0  0.959224  0.959224   0.513791  0.365621          0

ls test.*
test.stats.csv
```

The most relevant metrics here again are recall and precision.

Installation
============

### Hardware and OS Requirements

## Hardware

Compiling and testing can be done on a standard desktop system with 8GB of RAM. Whole-genome
comparisons (e.g. comparing a gVCF file against the [Platinum Genomes truth dataset](http://www.illumina.com/platinumgenomes/))
can use up to 64GB of RAM (20GB typical, depending on the input VCF) and about 10-20 minutes
using 32 processor cores. Whole exome comparison (using an exome bed mask and the `-T` switch)
can be carried out on a desktop system.

## Linux

Hap.py is known to build and run on the following linux distributions (see also the [Dockerfile](Dockerfile)
for a list of required packages):

    Ubuntu 12.04,14.04
    CentOS 5,6,7

## OS X

Hap.py builds and passes basic tests on OS X 10.9, but full WGS analyses are not tested for this platform.

## Windows

Hap.py is not tested on Windows. The main dependency that fails compilation is htslib. Given a build
of htslib and pysam, using hap.py on Windows should be possible.

### Requirements:

Firstly, hap.py requires a human reference sequence which contains at least
chromosomes 1-22,X,Y,M. The chromosomes should be named chr1-chr22, chrX, chrY,
chrM. there is a script  in [src/sh/make_hg19.sh](src/sh/make_hg19.sh) to create
such a sequence, but you can also  specify your own. In order for the
integration tests to run successfully, it is necessary  to point hap.py to the
reference sequence using

```bash
export HGREF=<path-to-hg19.fa>
```

Note that, while the test cases are based on hg19, other reference sequences are
usable as well  once the tool is installed.

Hap.py also requires a copy of the [Boost libraries](http://www.boost.org) to
work, with version >=  1.55. If compilation should fail using the included version
of boost, you can compile a subset of boost like this:

```bash
cd ~
wget http://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
tar xjf boost_1_55_0.tar.bz2
cd boost_1_55_0
./bootstrap.sh --with-libraries=filesystem,chrono,thread,iostreams,system,regex,test,program_options
./b2 --prefix=$HOME/boost_1_55_0_install install
```

You can point Cmake to your version of boost as follows:

```bash
export BOOST_ROOT=$HOME/boost_1_55_0_install
```

The complete list of dependencies / packages to install beforehand can be found
in the [Dockerfile](Dockerfile).

### Installation Procedure

There are two fast ways to get a running installation of hap.py:

1. Use the installer script. In the simplest use case, this script can create an
   installation of hap.py from source that uses the system Python. You will need
   tohave the following packages installed: Cython, numpy, pandas, pybedtools,
   pysam, bx-python.

   The simplest installer command line is the following, it installs everything into
   ~/hap.py-install using the system version of Python:

   ```
   python install.py ~/hap.py-install
   ```

   The installer has an option `--boost-root` that allows us to use a specific installation of boost
   (see above for instructions):

   ```
   python install.py ~/hap.py-install --boost-root $HOME/boost_1_55_0_install
   ```

   To use a special version of Python, run the installer with it:

   ```
   $HOME/my-virtualenv/bin/python install.py ~/hap.py-install
   ```

   To create a virtualenv, you can use the following options:

   ```
   python install.py ~/workspace-is/hap.py-install --python=virtualenv --python-virtualenv-dir=$HOME/my-virtualenv/hc.ve
   ```

   There are various workaround / testing switches:

   * `--python-virtualenv-update`  updates an existing virtualenv
   * `--python-virtualenv-force`  overwrites the virtualenv if it exists
   * `--pip-fix-cert` works around outdated SSL certificates when using pip
   * `--no-tests` disables the unit/integration tests after installation
   * `--no-rebuild-external` don't rebuild the external dependencies (htslib, ...) unless necessary
   * `--sge-mode` require switch `--force-interactive` to run hap.py
     interactively (useful to prevent running on a head node when installing on
     systems with SGE)

2. Use [Docker](https://www.docker.com/). Clone this repository and build a
   Docker image as follows.
   ```
   $ sudo docker build .
   $ sudo docker images
   REPOSITORY     TAG            IMAGE ID            CREATED             VIRTUAL SIZE
   <...>          latest         3d03a99b3d81        1 second ago        <...>
   $ sudo docker run -ti --rm 3d03a99b3d81 bin/bash
   $/ /opt/hap.py/bin/hap.py
   ```


Compiling
=========

This section shows how to compile hap.py from source without using the installer.

List of Dependencies
--------------------

You will need these tools / libraries on your system to compile the code.

* CMake &gt; 2.8
* GCC/G++ 4.8+ for compiling
* Boost 1.55+
* Python 2, version 2.7.8 or greater
* Python packages: Pandas, Numpy, pysam, bx-python

Compiling using CMake
---------------------

1.  Get a hap.py checkout:
    ```bash
    git clone https://github.com/sequencing/hap.py
    ```
2.  Make a build folder
    ```bash
    mkdir hap.py-build
    cd hap.py-build
    ```
3.  Run CMake
    ```bash
    cmake ../hap.py
    ```
4.  Build
    ```bash
    make
    ```

If this is successful, the bin subdirectory of your build folder will contain binaries and scripts:

```bash
$ python bin/hap.py --version
Hap.py v0.2.3
```

Additional Cmake build flags
----------------------------

The source for hap.py contains a script [configure.sh](configure.sh) which shows some basic additional
configuration flags, and an automated way to pre-package CMake setups.

Here is a list of additional flags for CMake to change compile options help it find dependencies:

*  `-DCMAKE_BUILD_TYPE=Debug` -- set the build type, allowed values are `Debug` and `Release`
*  `-DCMAKE_C_COMPILER=/usr/bin/gcc` and `-DCMAKE_CXX_COMPILER=/usr/bin/g++` -- change the compiler path
*  `-DCMAKE_INSTALL_PREFIX=/usr/local` -- set an installation directory that will be used by make install.
*  `-DBOOST_ROOT=$HOME/boost_1_55_0_install` -- set the path to Boost. Run the following commands to compile and install boost:
   ```bash
   cd ~
   wget http://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
   tar xjf boost_1_55_0.tar.bz2
   cd boost_1_55_0
   ./bootstrap.sh --with-libraries=filesystem,chrono,thread,iostreams,system,regex,test,program_options
   ./b2 --prefix=$HOME/boost_1_55_0_install install
   ```
*  `-DUSE_SGE` -- enable the `--force-interactive` switch in hap.py.
