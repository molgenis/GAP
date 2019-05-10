#!/bin/bash

set -e
set -u

# identify on which cluster this script will be executed
HOST=$(hostname -s)

# source config files

source ${EBROOTGAP}/automated/etc/${HOST}.cfg

#Checking if VCF File array is older then ArchiveTime and archive file if this is the case

for i in $(find "/groups/umcg-gd/${TMP_LFS}/Concordance/array/" -type f -maxdepth 1 -mtime +31 -print)
do
echo "moving ${i} to the neverUsed folder since 31 days has passed and the Concordance has not been executed yet" >> /groups/umcg-gd/${TMP_LFS}/logs/Archive_Concordance_log.txt
 mv "${i}" "/groups/umcg-gd/${TMP_LFS}/Concordance/array/neverUsed/"
done

#Checking if VCF File NGS is older then ArchiveTime and archive file if this is the case

for i in $(find "/groups/umcg-gd/${TMP_LFS}/Concordance/ngs/" -type f -maxdepth 1 -mtime +31 -print)
do
	echo "moving ${i} to the neverUsed folder since 31 days has passed and the Concordance has not been executed yet" >> /groups/umcg-gd/${TMP_LFS}/logs/Archive_Concordance_log.txt
	mv "${i}" "/groups/umcg-gd/${TMP_LFS}/Concordance/ngs/neverUsed/"
done
