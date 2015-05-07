#!/bin/bash
#
# Download HG19 reference chr1-22,X,Y,M
#
# Author: Peter Krusche <pkrusche@illumina.com>

##############################################################
# Create Simplified HG19 Reference
##############################################################

set -e

PREFIX="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes"

if [[ -e hg19.fa ]]; then
	echo "hg19.fa already exists, delete to redownload"
	exit 0
fi

for x in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" "M" ; do
	echo "wget ${PREFIX}/chr${x}.fa.gz -O chr${x}.fa.gz"
	wget ${PREFIX}/chr${x}.fa.gz -O chr${x}.fa.gz
	gunzip -c "chr${x}.fa.gz" >> hg19.fa
	rm "chr${x}.fa.gz"
done
