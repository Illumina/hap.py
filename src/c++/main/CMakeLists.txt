# multimerge is our version of bcftools merge
add_executable(multimerge multimerge.cpp)
target_link_libraries(multimerge ${HAPLOTYPES_ALL_LIBS})

# hapenum takes a VCF and outputs all possible haplotypes
add_executable(hapenum hapenum.cpp)
target_link_libraries(hapenum ${HAPLOTYPES_ALL_LIBS})

# dipenum takes a VCF and outputs paired diploid haplotypes
add_executable(dipenum dipenum.cpp)
target_link_libraries(dipenum ${HAPLOTYPES_ALL_LIBS})

# blocksplit takes a set of VCFs and regions which do not break haplotype blocks at the edges
add_executable(blocksplit blocksplit.cpp ../include/helpers/Roc.hh)
target_link_libraries(blocksplit ${HAPLOTYPES_ALL_LIBS})

# hapcmp takes a pair of VCFs and outputs comparison results based on alignments
add_executable(hapcmp hapcmp.cpp)
target_link_libraries(hapcmp ${HAPLOTYPES_ALL_LIBS})

# xcmp takes a pair of VCFs and outputs an annotated comparison VCF (by genotypes)
add_executable(xcmp xcmp.cpp)
target_link_libraries(xcmp ${HAPLOTYPES_ALL_LIBS})

# scmp takes a pair of VCFs and outputs an annotated comparison VCF (by alleles)
add_executable(scmp scmp.cpp)
target_link_libraries(scmp ${HAPLOTYPES_ALL_LIBS})

# alleles is a pre-processor for scmp. It takes all samples, writes one record per allele
# collects FORMATs into INFO and GT.
add_executable(alleles alleles.cpp)
target_link_libraries(alleles ${HAPLOTYPES_ALL_LIBS})

# quantify counts variants in VCFs
add_executable(quantify quantify.cpp)
target_link_libraries(quantify ${HAPLOTYPES_ALL_LIBS})

# vcfhdr2json turns a VCF header into JSON
add_executable(vcfhdr2json vcfhdr2json.cpp)
target_link_libraries(vcfhdr2json ${HAPLOTYPES_ALL_LIBS})

# ROC takes a table with labels and a value + calculates a ROC curve
add_executable(roc roc.cpp)
target_link_libraries(roc ${HAPLOTYPES_ALL_LIBS})

# validatevcf checks ref alleles and tests if haplotypes can be enumerated
add_executable(validatevcf validatevcf.cpp)
target_link_libraries(validatevcf ${HAPLOTYPES_ALL_LIBS})

# vcfcheck tests if something will be acceptable by hap.py
# this is different from validatevcf, which counts the number of
# superloci that will cause matching problems.
add_executable(vcfcheck vcfcheck.cpp)
target_link_libraries(vcfcheck ${HAPLOTYPES_ALL_LIBS})

# preprocess does variant decomposition and leftshifting
add_executable(preprocess preprocess.cpp)
target_link_libraries(preprocess ${HAPLOTYPES_ALL_LIBS})

# fastainfo returns contig lengths and non-n lengths in JSON format
add_executable(fastainfo fastainfo.cpp)
target_link_libraries(fastainfo ${HAPLOTYPES_ALL_LIBS})

# gVCF2bed converts a gVCF file to bed, indicating regions of pass/fail homref/variants
add_executable(gvcf2bed gvcf2bed.cpp)
target_link_libraries(gvcf2bed ${HAPLOTYPES_ALL_LIBS})
