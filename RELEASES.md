# Hap.py Release Notes / Change Log

## v0.2.8

P3 (Major)  Bug HAP-147 Normalisation problem where het-alts cancel each other out 
P3 (Major)  Bug HAP-151 handle INDELs of undetermined length   
P3 (Major)  Task    HAP-152 Implement Ti/Tv computation
P3 (Major)  Task    HAP-155 Improve support for GRC37 / non-chr-prefix references  
P3 (Major)  Task    HAP-157 Calculate per-MB FP rate (as an "absolute" equivalent to precision)s   
P3 (Major)  Bug HAP-162 Som.py outputs INF precision when FP and TP are 0  
P3 (Major)  Task    HAP-163 Release notes and final test   
P3 (Major)  Task    HAP-164 Sync up som.py empirical scoring features
P3 (Major)  Task    HAP-165 Speed improvements in quantify and upgrade to htslib 1.3

## v0.2.7

P3 (Major)  Bug HAP-133 Parse platypus version strings correctly    
P3 (Major)  New Feature HAP-134 Calculate AFs from BAM for som.py feature table 
P3 (Major)  Bug HAP-135 Som.py labels some variants twice in admixture comparisons  
P3 (Major)  Improvement HAP-136 Better Strelka feature extraction   
P3 (Major)  Bug HAP-137 GT ordering makes unit test fail    
P3 (Major)  Task    HAP-138 Add het/hom ratio to output table   
P3 (Major)  Task    HAP-139 Document hap.py in Docker usage 
P3 (Major)  Bug HAP-140 Misleading error messages when checkout headers in hap.py   
P3 (Major)  Improvement HAP-141 Remove / improve reference check error message  
P3 (Major)  Task    HAP-142 Som.py -- stratify by AF brackets   
P3 (Major)  Bug HAP-143 Empty VCF files cause failure   
P3 (Major)  Improvement HAP-144 Add filtered FN vs. absent FN to result table   
P3 (Major)  Improvement HAP-145 Hap.py -o switch doesn't accept prefix in current directory 

## v0.2.6

P3 (Major)  Bug HAP-83  premature exit when multimerging gvcfs 
P3 (Major)  Bug HAP-123 Varscan2 DP rates are not calculated correctly 
P3 (Major)  Bug HAP-124 Som.py --fix-chr-truth doesn't fix chr names in the ambi bed files 
P3 (Major)  Improvement HAP-126 Add VCF Validation Step
P3 (Major)  Improvement HAP-127 Make -X the default behaviour (write extended metrics unless told otherwise)  
P3 (Major)  Improvement HAP-128 Hap.py ROCs for Freebayes and Platypus 
P3 (Major)  Task    HAP-129 Som.py should support FP regions from Ambi file

## v0.2.5

P3 (Major)  Bug HAP-119 ROCs should use TP variant types from query, not truth 
P3 (Major)  Improvement HAP-120 Write ROC data into JSON output
P3 (Major)  Task    HAP-122 Integrate VCFEval for matching 
P4 (Minor)  Bug HAP-121 Installer fails on some systems

## v0.2.4

P3 (Major)  New Feature HAP-103 Germline ROC curves
P3 (Major)  Bug HAP-115 --no-internal-* options are broken 
P3 (Major)  Bug HAP-116 Quantify errors don't get printed  
P3 (Major)  Bug HAP-117 Raw count quantify calls don't respect locations passed to hap.py  

## v0.2.3

P3 (Major)  Bug HAP-106 blocksplit / quantify segfault when reading faulty VCF 
P3 (Major)  Bug HAP-107 Som.py tests fail in SD 
P3 (Major)  Bug HAP-109 Hap.py fails for ambiguous bases in reference fasta 
P3 (Major)  Task    HAP-111 Remove muscle dependency code and legacy stuff that depends on it         
P3 (Major)  Improvement HAP-112 Add option for hap.py to disable haplotype matching 
P3 (Major)  Improvement HAP-113 Check for output folder at start    
P4 (Minor)  Bug HAP-79  Well-formed vcf output file 
P4 (Minor)  Improvement HAP-100 Clean up tests, fix platform-dependent sorting issues   

## v0.2.2

P3 (Major)  Task    HAP-23  Multi-sample variant graph implementation   
P3 (Major)  Task    HAP-52  Create runner / wrapper script to run inside Docker 
P3 (Major)  Improvement HAP-78  Reference FASTA handling    
P3 (Major)  Task    HAP-88  Refactor Cmake files to build included Boost, move build of dependencies to build folder    
P3 (Major)  Task    HAP-90  Implement parallel graph enumeration    
P3 (Major)  Improvement HAP-92  Quantify should annotate the output VCF according to actual status          
P3 (Major)  Improvement HAP-94  Add allele counts by type (rather than per NT) back in  
P3 (Major)  New Feature HAP-95  Stratified feature output   
P3 (Major)  Bug HAP-96  cnx aligners empty  
P3 (Major)  Bug HAP-97  Hap.py does not quote or escape file path internally and fails on paths with spaces         
P3 (Major)  Improvement HAP-99  som.py should also be able to fix truth chr names   
P3 (Major)  Improvement HAP-102 Eliminate / clear up Locations.unknown  
P3 (Major)  Task    HAP-104 Fix samtools build when using modules, add eb files 

## v0.2.1

P3 (Major)  Bug HAP-86  Not all vcfeval sites are matched  
P3 (Major)  Bug HAP-87  auto chr naming detection broken   
P3 (Major)  Task    HAP-89  Treatment of filtered variants should not be as optional   
P3 (Major)  Bug HAP-93  GT mismatches are not counted as FP when no confident regions are specified.  
P4 (Minor)  Bug HAP-82  Print error (or at least warning) for bogus command-line arguments 

## v0.2.0

P3 (Major)  New Feature HAP-39  Implement Allelic-primitive variant splitting  
P3 (Major)  Bug HAP-72  QUAL misleading for type=missing records   
P3 (Major)  Improvement HAP-74  Remove setup.py
P3 (Major)  Task    HAP-75  Document the fact that bed files with tracks aren't supported. 
P3 (Major)  New Feature HAP-76  Add version information to output / simple output  
P3 (Major)  Task    HAP-77  som.py should write puma metrics JSON file also
P3 (Major)  Bug HAP-80  Handle -P switch correctly for silver variants in PG 8.0   
P3 (Major)  New Feature HAP-81  Support VQSR ROCs in som.py
P3 (Major)  Improvement HAP-84  no hap.py versions in output   
P4 (Minor)  Improvement HAP-71  Hap.py fails on BED files with track information   

## v0.1.6

P3 (Major)  Task    HAP-8   Write user documentation    
P3 (Major)  Bug HAP-54  Version number does not display 
P3 (Major)  Task    HAP-59  Change som.py test to PG Admixture data 
P3 (Major)  Bug HAP-62  hap.py non-verbose fail when FP regions file not found  
P3 (Major)  Bug HAP-63  Add more comprehensive PG hap.py test (PGv7 vs. GATK 1.6 on chr21)  

## v0.1.5

P3 (Major)  New Feature HAP-22  Split ref-match blocks into SNPs and indels
P3 (Major)  Improvement HAP-44  Clean up + Migrate Dockerfile to Ubuntu 14.04  
P3 (Major)  Task    HAP-46  Test + add GA4GH difficult indels to test set   
P3 (Major)  Task    HAP-49  som.py should do ROC curves
P3 (Major)  Improvement HAP-56  Define value in metrics.json   
P3 (Major)  Bug HAP-57  -R option can make xcmp fail   
P3 (Major)  Improvement HAP-58  Add option for verbosity levels or to write message to log file
P3 (Major)  Bug HAP-54  Version number does not display 
P3 (Major)  Bug HAP-62  hap.py non-verbose fail when FP regions file not found  
P3 (Major)  Bug HAP-63  Add more comprehensive PG hap.py test (PGv7 vs. GATK 1.6 on chr21)  
