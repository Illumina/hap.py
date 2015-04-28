Haplotype Comparison Tools
==========================

Peter Krusche <pkrusche@illumina.com>

This is a set of programs based on [htslib](https://github.com/samtools/htslib)
to compare VCF files by specified haplotype.

Rather than comparing entries individually, we produce a graph-based reference 
from the VCF entries, create all possible haplotype sequences, and compare 
these by alignment / exact matching. This is more accurate in cases like this:

![](doc/rep_ex.PNG)

As an output, we produce the haplotypes / paths through the reference graph 
which are consistent between the input VCFs.

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

Hap.py is set up to run inside SGE. Unless forced, it will not run interactively
(it detects this by looking for the environment variable SGE_JOB_ID).

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
work, with version >=  1.55. If compilation should fail using your system-wide
installation of boost, you can compile a  subset of boost like this:

```bash
cd ~
wget http://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
tar xjf boost_1_55_0.tar.bz2
cd boost_1_55_0
./bootstrap.sh --with-libraries=filesystem,chrono,thread,iostreams,system,regex,test,program_options
./b2 --prefix=$HOME/boost_1_55_0_install install
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
      
   The installer has an option `--boost-root` that allows us to use a specific installation of boost:

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

Known Issues
============

*  The testing scripts aren't working under MacOS X -- needs some shell tweaking
   [src/sh/run_tests.sh](src/sh/run_tests.sh)
   