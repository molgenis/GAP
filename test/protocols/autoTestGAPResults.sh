#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=8gb ppn=1

#string tmpName
#string	Project
#string logsDir
#string intermediateDir

set -e
set -u

module load HTSlib
module load BEDTools
module load ngs-utils
module list

## Filtering of the VCF files from the pipeline

#1 Compress the VCF file


for i in  $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/run01/results/vcf/*.vcf)
do
	file=$(basename ${i})
	mkdir -p /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/
	echo  "zipping vcf file : ${file} ..."
	bgzip -c ${i} > /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/${file}.gz

done

#2 making the index of the VCF file

for i in $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/*.gz)
do
	echo "indexing the vcf file of ${i} ..."
	tabix -p vcf ${i}
done

#3 Filtering the VCF using the bedfile

for i in $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/*.gz)
do
	echo "filtering the vcf file: ${i} ..."
	file=$(basename ${i})
	sample=$(basename ${file} .FINAL.vcf.gz)
	bedtools intersect -a ${i} -b /home/umcg-molgenis/GAP/autoTestArray.bed > /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/${sample}.FINAL_FILTERED.vcf
done

## Comparing the VCF files

for i in $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/*.FINAL_FILTERED.vcf)
do
	file=$(basename ${i})
	echo ${file}
	sample=$(basename ${file} .FINAL_FILTERED.vcf)
	echo "Comparing the VCF file with the TRUE VCF for the sample sample : ${sample} ..."
	${EBROOTNGSMINUTILS}/vcf-compare_2.0.sh -1 "/home/umcg-molgenis/GAP/vcf/${sample}.FINAL_TRUE_FILTERED.vcf" -2 "${i}" -o "/groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/VCF_Compare/${sample}/"
done

## Checking if the output is correct

for i in $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/*.FINAL_FILTERED.vcf)
do
	file=$(basename ${i})
	echo ${file}
	sample=$(basename ${file} .FINAL_FILTERED.vcf)
	vcfCheckvalue="$(awk 'NR == 2 {print $4}' /groups/umcg-gsad/tmp03/projects/NIST_TRIO/temp/VCF_Compare/${sample}/vcfStats.txt)"
	echo "${sample}: vcf Check value is: ${vcfCheckvalue}"

	if [ ${vcfCheckvalue} = 100.00% ]
	then
		echo "VCF files are correct for ${sample}"
	else
		echo "there are differences in the VCF files between the test and the original data"
		echo "please fix the bug or update this test"
		exit 1
	fi
done


## Checking the PennCNV per sample files
echo "start comparing the PennCNV per sample files ..."

for i in $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/run01/results/PennCNV_reports/*.txt)
do
	file=$(basename ${i})
	sample=$(basename ${file} .txt)

	echo "Comparing the PennCNV per sample file for ${sample}"
	diffPennCNV="false"

	diff -q  <(tail -n +4 /home/umcg-molgenis/GAP/PennCNV_reports/${sample}_TRUE.txt) <(tail -n +4 ${i}) || diffPennCNV="true"
#	diff -q /home/umcg-molgenis/GAP/PennCNV_reports/${sample}_TRUE.txt ${i} || diffPennCNV="true"
	if [ "${diffPennCNV}" == "true" ]
	then
		echo "there are differences in the PennCNV files between the test and original data for sample ${sample}"
		echo "please fix the bug or update this test"
		exit 1
	else
	echo "there are no differences between the PennCNV per sample files for sample ${sample}"
	fi
done


## Checking the SD files per sample

for i in $(ls /groups/umcg-gsad/tmp03/projects/NIST_TRIO/run01/results/vcf/*.sd)
do
	file=$(basename ${i})
        sample=$(basename ${file} .vcf.sd)
	echo ${sample}
        echo "Comparing the standard deviation's per sample for ${sample}"
        diffSD="false"
        diff -q /home/umcg-molgenis/GAP/vcf/${sample}_TRUE.vcf.sd ${i} || diffPennCNV="true"
        if [ "${diffSD}" == "true" ]
        then
		echo "there are differences in the SD files between the test and original data for sample ${sample}"
                echo "please fix the bug or update this test"
                exit 1
        else
	echo "there are no differences between the SD files per sample, for sample ${sample}"
        fi
done



## Checking the PennCNV project file
echo "Comparing the Project PennCNV files ... "

diffPennCNVProjectFile="false"

diff -q /home/umcg-molgenis/GAP/NIST_TRIO_PennCNV_TRUE.txt /groups/umcg-gsad/tmp03/projects/NIST_TRIO/run01/results/NIST_TRIO_PennCNV.txt || diffPennCNVProjectFile="true"
        if [ "${diffPennCNVProjectFile}" == "true" ]
        then
		echo "there are differences in the PennCNV Project files between the test and original data for sample ${sample}"
                echo "please fix the bug or update this test"
                exit 1
        else
	echo "there are no differences between the PennCNV Project files"
        fi

## Checking the callrate file
echo "Comparing the Project callrate files ... "

diffProjectCallrateFile="false"

diff -q /home/umcg-molgenis/GAP/results/Callrates_NIST_TRIO_TRUE.txt /groups/umcg-gsad/tmp03/projects/NIST_TRIO/run01/results/Callrates_NIST_TRIO.txt || diffProjectCallrateFile="true"
        if [ "${diffProjectCallrateFile}" == "true" ]
        then
                echo "there are differences in the Project Callrate files between the test and original data for sample ${sample}"
                echo "please fix the bug or update this test"
                exit 1
        else
	echo "there are no differences between the Project Callrate files"
        fi
