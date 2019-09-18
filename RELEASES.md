# Hap.py Release Notes / Change Log

## v0.3.12

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-355 | Ensure compatibility with genome VCF from DRAGEN pipelines                         |
| HAP-356 | quantify module throws regex error                                                 |
| HAP-357 | Deal with change in Java license terms                                             |
| HAP-359 | Hap.py crashes when ingesting variants genotyped as <NON_REF>                      |
| HAP-360 | Apply non-ref filter to truth set when --preprocess-truth is set                   |
| HAP-361 | Hap.py "blocksplit" processes still failing on DRAGEN gVCFs even with pre-filtering NON_REF genotypes 
| HAP-362 | Hap.py blocksplit sometimes crashes with "invalid next size (fast)" error          |
| HAP-363 | Use 'git describe' to obtain version number during Cmake                           |

## v0.3.11

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-341 | Update rtgtools to 3.10.1                                                          |
| HAP-342 | som.py output tables to quantify / GA4GH format                                    |
| HAP-346 | src/sh/illumina-setup.sh out of date                                               |
| HAP-352 | string formatting in som.py: last interval is always [1.0, 1.0] in *extended.csv   |

## v0.3.10

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-326 | Support Pisces VCFs in som.py                                                      |
| HAP-332 | Fix bug in install.py related to external dependencies                             |
| HAP-333 | Update rtgtools to 3.8.2 and add support for distance-based matching               |
| HAP-336 | Change Dataframe.sort to sort_value in som.py (#21)                                |

## v0.3.9

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-317 | Update htslib to 1.4.1 and always check for BCF conversion errors                  |
| HAP-318 | Update rtgtools to 3.8.1                                                           |
| HAP-319 | Slimmer docker image without hg19 reference                                        |

## v0.3.8

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-291 | Upgrade rtgtools dependency to version 3.7.1                                       |
| HAP-292 | Improve Dockerfile and update compiler requirements                                |
| HAP-295 | Preserve INFO fields when running with vcfeval comparison engine also              |
| HAP-296 | Don't fail when no reference is found in default locations, just require `-r`      |
| HAP-298 | Output run and session info in prefix.runinfo.json                                 |
| HAP-299 | Improve command line documentation for confident regions.                          |
| HAP-300 | Correct documentation for `--adjust-conf-regions`                                  |
| HAP-306 | Add distance-based matching method to assess discovery separately from genotyping  |

## v0.3.7

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-195 | Safer left-shifting and test for leftshifting-related issues added                 |
| HAP-259 | Som.py alternative GA4GH-based workflow: hap.py --engine scmp                      |
| HAP-266 | Optional new insertion handling: both surrounding bases must be covered to capture |
| HAP-267 | Add confident region overlap counts to stratification regions                      |
| HAP-268 | Fix bamstats.py (missing imports, now working)                                     |
| HAP-270 | Handle special chars in input filenames better                                     |
| HAP-274 | Issue #10 allow hap.py to be symlinked                                             |
| HAP-275 | Somatic mode for hap.py comparisons added: --somatic switch squashes ploidy        |
| HAP-276 | List-of-float INFO fields were getting rounded to int                              |
| HAP-279 | Issue #11 Output error message when preprocess fails                               |
| HAP-280 | Add Subset.Level column to stratified output, stratify into exons and genes        |
| HAP-281 | Quantify precision and recall separately at truthset boundaries                    |
| HAP-282 | Output stratified counts also when regions are empty                               |
| HAP-285 | Hapenum outputs individual variants used in haplotype construction                 |
| HAP-286 | MacOS X build failure fixed                                                        |
| HAP-287 | Report F-score in summary CSV and also in summary output                           |
| HAP-288 | Read vcfeval SDF file from same path as fasta file if present                      |

## v0.3.6

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-258 | Som.py avoid CI failure for empty subsets                                          |
| HAP-260 | Add Dockerfile based on Centos 6                                                   |
| HAP-261 | Som.py FP rate calculation only uses truth contig sizes                            |
| HAP-263 | Som.py fix AF binning (final bin now complete)                                     |
| HAP-264 | Som.py support SomaticEVS feature for Strelka somatic                              |
| HAP-265 | Hap.py haplotype-matched TPs now assigned min qual from across the superlocus      |

## v0.3.5

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-248 | qfy.py fails when -V is not on                                                     |
| HAP-249 | stratification regions: include length of stratification region in *extended.csv   |
| HAP-250 | stratification regions: allow to summarise performance by ID (4th column)          |
| HAP-254 | Fix haploid GTs on chrX for male samples.                                          |
| HAP-255 | support * alleles (and alleles that contain "*") to support octopus VCFs           |
| HAP-256 | Sentieon caller extraction                                                         |
| HAP-257 | Fix PASS recall for complex matches                                                |

## v0.3.4

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-240 | Fix division by zero in som.py CI computation                                      |
| HAP-241 | Add genotype validation to vcfcheck. Will now fail when alleles don't exist.       |
| HAP-242 | Improved counts in bamstats.py                                                     |
| HAP-244 | Support the FT format field to enable per-sample filtering                         |

## v0.3.3

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-232 | Fix help text in hap.py and pre.py                                                 |
| HAP-233 | INFO fields don't get mixed up in preprocessing anymore                            |
| HAP-234 | Forward unstructured headers through hap.py                                        |
| HAP-236 | Fix quanitify test to work with compressed outputs                                 |

## v0.3.2

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-180 | Confidence intervals for precision and recall in som.py                            |
| HAP-221 | document how stratification region overlap is calculated                           |
| HAP-222 | document all output columns in result CSV files                                    |
| HAP-223 | getFormatInt is broken / only returns DP, which leads to incorrect ROCs            |
| HAP-225 | System BOOST install can break installation process. Only use included version.    |
| HAP-226 | cnx gets wrong caller version for GATK in some cases                               |
| HAP-227 | QQ.Field should not be IQQ when using xcmp for comparison                          |
| HAP-229 | Hap.py fails when query does not contain all default chromosomes                   |
| HAP-230 | Compress large output files to save space when writing ROCs to JSON                |

## v0.3.1

| Ticket  | Description                                                                        |
|---------|------------------------------------------------------------------------------------|
| HAP-175 | Enable GQX as a ROC feature                                                        |
| HAP-184 | check for empty query/ref vcf                                                      |
| HAP-185 | som.py --fix-chr-query does not fix MT appropriately                               |
| HAP-187 | som.py automatically resolve whether the --fix-chr-query switch is required        |
| HAP-189 | refactor storage of INFO and FORMAT, extract preprocessing into separate script    |
| HAP-192 | Missing chromosomes cause fp.rate computation to fail in som.py                    |
| HAP-197 | AF binning with arbitrary boundaries in som.py                                     |
| HAP-201 | hap.py integration tests failing because of trailing white space                   |
| HAP-202 | F-score computation is not correct                                                 |
| HAP-204 | When no / wrong command line is specified, hap.py should have non-zero return code |
| HAP-205 | fix EVS feature support in som.py                                                  |
| HAP-206 | Fix header for BLT in hap.py output                                                |
| HAP-207 | Count TP, FN and FP outside confident regions as UNK                               |
| HAP-208 | Simpler deployment with vcfeval                                                    |
| HAP-209 | Update documentation and release notes                                             |
| HAP-210 | Calculate truth-QQ values as max over all query BS TPs                             |
| HAP-211 | Output BS in xcmp                                                                  |
| HAP-212 | Improve ROC filter handling / add selective filtering                              |
| HAP-213 | Update Dockerfile                                                                  |
| HAP-215 | include FDP in som.py indel features                                               |
| HAP-217 | pad symbolic deletions in pre-processing for vcfeval                               |
| HAP-218 | improve pre.py command line usage                                                  |
| HAP-219 | Don't mix variants with different filtering during preprocessing                   |
| HAP-220 | Quantify counts hom-ref records as INDEL for UNK and for N                         |

## v0.3.0

| Ticket   | Description                                                                    |
|----------|--------------------------------------------------------------------------------|
| HAP-181  | Need documentation with detailed description of hap.py output                  |
| HAP-188  | Support GA4GH stratification regions                                           |
| HAP-190  | Implement final output format for 0.3.0                                        |
| HAP-193  | Implement indel length binning                                                 |

## v0.2.9

| Ticket  | Description                                                                  |
|---------|------------------------------------------------------------------------------|
| HAP-149 | Fix vcfeval integration via GA4GH interface                                  |
| HAP-154 | ROC Support for Mutect2                                                      |
| HAP-161 | Output GT precision and FP counts for wrong GTs and wrong alleles separately |
| HAP-166 | Indel haplotype enumeration bug                                              |
| HAP-169 | Change quantify/xcmp output be GA4GH compliant                               |
| HAP-170 | Som.py ROCs not correct in some cases                                        |
| HAP-171 | ROC support using vcfeval                                                    |
| HAP-176 | One run for PASS and ALL to fix ROCs                                         |
| HAP-177 | CNX doesn't extract caller name correctly for Mutect2                        |
| HAP-179 | Symmetric FP and FN for genotype mismatch                                    |
| HAP-183 | Handling GATK gVCF \<NON\_REF\> alleles                                      |
| HAP-186 | Update documentation and release notes                                       |


## v0.2.8

| Ticket   | Description                                                                    |
|----------|--------------------------------------------------------------------------------|
| HAP-147  | Normalisation problem where het-alts cancel each other out                     |
| HAP-151  | handle INDELs of undetermined length                                           |
| HAP-152  | Implement Ti/Tv computation                                                    |
| HAP-155  | Improve support for GRC37 / non-chr-prefix references                          |
| HAP-157  | Calculate per-MB FP rate (as an "absolute" equivalent to precision)s           |
| HAP-162  | Som.py outputs INF precision when FP and TP are 0                              |
| HAP-163  | Release notes and final test                                                   |
| HAP-164  | Sync up som.py empirical scoring features                                      |
| HAP-165  | Speed improvements in quantify and upgrade to htslib 1.3                       |

## v0.2.7

| Ticket   | Description                                                                    |
|--------- | ------------------------------------------------------------------------------ |
| HAP-133  | Parse platypus version strings correctly                                       |
| HAP-134  | Calculate AFs from BAM for som.py feature table                                |
| HAP-135  | Som.py labels some variants twice in admixture comparisons                     |
| HAP-136  | Better Strelka feature extraction                                              |
| HAP-137  | GT ordering makes unit test fail                                               |
| HAP-138  | Add het/hom ratio to output table                                              |
| HAP-139  | Document hap.py in Docker usage                                                |
| HAP-140  | Misleading error messages when checkout headers in hap.py                      |
| HAP-141  | Remove / improve reference check error message                                 |
| HAP-142  | Som.py -- stratify by AF brackets                                              |
| HAP-143  | Empty VCF files cause failure                                                  |
| HAP-144  | Add filtered FN vs. absent FN to result table                                  |
| HAP-145  | Hap.py -o switch doesn't accept prefix in current directory                    |

## v0.2.6

| Ticket   | Description                                                                    |
|----------|--------------------------------------------------------------------------------|
| HAP-83   | premature exit when multimerging gvcfs                                         |
| HAP-123  | Varscan2 DP rates are not calculated correctly                                 |
| HAP-124  | Som.py --fix-chr-truth doesn't fix chr names in the ambi bed files             |
| HAP-126  | Add VCF Validation Step                                                        |
| HAP-127  | Make -X the default behaviour (write extended metrics unless told otherwise)   |
| HAP-128  | Hap.py ROCs for Freebayes and Platypus                                         |
| HAP-129  | Som.py should support FP regions from Ambi file                                |

## v0.2.5

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-119   | ROCs should use TP variant types from query, not truth                           |
| HAP-120   | Write ROC data into JSON output                                                  |
| HAP-122   | Integrate VCFEval for matching                                                   |
| HAP-121   | Installer fails on some systems                                                  |

## v0.2.4

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-103   | Germline ROC curves                                                              |
| HAP-115   | --no-internal-* options are broken                                               |
| HAP-116   | Quantify errors don't get printed                                                |
| HAP-117   | Raw count quantify calls don't respect locations passed to hap.py                |

## v0.2.3

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-106   | blocksplit / quantify segfault when reading faulty VCF                           |
| HAP-107   | Som.py tests fail in SD                                                          |
| HAP-109   | Hap.py fails for ambiguous bases in reference fasta                              |
| HAP-111   | Remove muscle dependency code and legacy stuff that depends on it                |
| HAP-112   | Add option for hap.py to disable haplotype matching                              |
| HAP-113   | Check for output folder at start                                                 |
| HAP-79    | Well-formed vcf output file                                                      |
| HAP-100   | Clean up tests, fix platform-dependent sorting issues                            |

## v0.2.2

| Ticket    | Description                                                                              |
|-----------|------------------------------------------------------------------------------------------|
| HAP-23    | Multi-sample variant graph implementation                                                |
| HAP-52    | Create runner / wrapper script to run inside Docker                                      |
| HAP-78    | Reference FASTA handling                                                                 |
| HAP-88    | Refactor Cmake files to build included Boost, move build of dependencies to build folder |
| HAP-90    | Implement parallel graph enumeration                                                     |
| HAP-92    | Quantify should annotate the output VCF according to actual status                       |
| HAP-94    | Add allele counts by type (rather than per NT) back in                                   |
| HAP-95    | Stratified feature output                                                                |
| HAP-96    | cnx aligners empty                                                                       |
| HAP-97    | Hap.py does not quote or escape file path internally and fails on paths with spaces      |
| HAP-99    | som.py should also be able to fix truth chr names                                        |
| HAP-102   | Eliminate / clear up Locations.unknown                                                   |
| HAP-104   | Fix samtools build when using modules, add eb files                                      |

## v0.2.1

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-86    | Not all vcfeval sites are matched                                                |
| HAP-87    | auto chr naming detection broken                                                 |
| HAP-89    | Treatment of filtered variants should not be as optional                         |
| HAP-93    | GT mismatches are not counted as FP when no confident regions are specified.     |
| HAP-82    | Print error (or at least warning) for bogus command-line arguments               |

## v0.2.0

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-39    | Implement Allelic-primitive variant splitting                                    |
| HAP-72    | QUAL misleading for type=missing records                                         |
| HAP-74    | Remove setup.py                                                                  |
| HAP-75    | Document the fact that bed files with tracks aren't supported.                   |
| HAP-76    | Add version information to output / simple output                                |
| HAP-77    | som.py should write puma metrics JSON file also                                  |
| HAP-80    | Handle -P switch correctly for silver variants in PG 8.0                         |
| HAP-81    | Support VQSR ROCs in som.py                                                      |
| HAP-84    | no hap.py versions in output                                                     |
| HAP-71    | Hap.py fails on BED files with track information                                 |

## v0.1.6

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-8     | Write user documentation                                                         |
| HAP-54    | Version number does not display                                                  |
| HAP-59    | Change som.py test to PG Admixture data                                          |
| HAP-62    | hap.py non-verbose fail when FP regions file not found                           |
| HAP-63    | Add more comprehensive PG hap.py test (PGv7 vs. GATK 1.6 on chr21)               |

## v0.1.5

| Ticket    | Description                                                                      |
|-----------|----------------------------------------------------------------------------------|
| HAP-22    | Split ref-match blocks into SNPs and indels                                      |
| HAP-44    | Clean up + Migrate Dockerfile to Ubuntu 14.04                                    |
| HAP-46    | Test + add GA4GH difficult indels to test set                                    |
| HAP-49    | som.py should do ROC curves                                                      |
| HAP-56    | Define value in metrics.json                                                     |
| HAP-57    | -R option can make xcmp fail                                                     |
| HAP-58    | Add option for verbosity levels or to write message to log file                  |
| HAP-54    | Version number does not display                                                  |
| HAP-62    | hap.py non-verbose fail when FP regions file not found                           |
| HAP-63    | Add more comprehensive PG hap.py test (PGv7 vs. GATK 1.6 on chr21)               |

