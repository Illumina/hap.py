# multimerge is our version of bcftools merge
add_executable(multimerge multimerge.cpp)
target_link_libraries(multimerge ${HAPLOTYPES_ALL_LIBS})

# hapenum takes a VCF and outputs all possible haplotypes
add_executable(hapenum hapenum.cpp)
target_link_libraries(hapenum ${HAPLOTYPES_ALL_LIBS})

# dipenum takes a VCF and outputs paired diploid haplotypes
add_executable(dipenum dipenum.cpp)
target_link_libraries(dipenum ${HAPLOTYPES_ALL_LIBS})

# hapblockfind takes a set of VCFs and outputs potential haplotype block regions
add_executable(hapblockfind hapblockfind.cpp)
target_link_libraries(hapblockfind ${HAPLOTYPES_ALL_LIBS})

# blocksplit takes a set of VCFs and regions which do not break haplotype blocks at the edges
add_executable(blocksplit blocksplit.cpp)
target_link_libraries(blocksplit ${HAPLOTYPES_ALL_LIBS})

# hapcmp takes a pair of VCFs and outputs comparison results based on alignments
add_executable(hapcmp hapcmp.cpp)
target_link_libraries(hapcmp ${HAPLOTYPES_ALL_LIBS})

# xcmp takes a pair of VCFs and outputs an annotated comparison VCF
add_executable(xcmp xcmp.cpp)
target_link_libraries(xcmp ${HAPLOTYPES_ALL_LIBS})

# quantify counts variants in VCFs
add_executable(quantify quantify.cpp)
target_link_libraries(quantify ${HAPLOTYPES_ALL_LIBS})

# denovo counts variants in only one sample
add_executable(denovo denovo.cpp)
target_link_libraries(denovo ${HAPLOTYPES_ALL_LIBS})

# vcfhdr2json turns a VCF header into JSON
add_executable(vcfhdr2json vcfhdr2json.cpp)
target_link_libraries(vcfhdr2json ${HAPLOTYPES_ALL_LIBS})

# ROC takes a table with labels and a value + calculates a ROC curve
add_executable(roc roc.cpp)
target_link_libraries(roc ${HAPLOTYPES_ALL_LIBS})

# postvcfeval merges TP/FP/FN files into a single file that is understood by quantify
add_executable(postvcfeval postvcfeval.cpp)
target_link_libraries(postvcfeval ${HAPLOTYPES_ALL_LIBS})

# validatevcf checks ref alleles and tests if haplotypes can be enumerated
add_executable(validatevcf validatevcf.cpp)
target_link_libraries(validatevcf ${HAPLOTYPES_ALL_LIBS})
