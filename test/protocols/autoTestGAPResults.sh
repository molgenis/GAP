#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=8gb ppn=1

#string tmpName
#string	Project
#string logsDir
#string intermediateDir
#string HTSlibVersion
#string bedToolsVersion
#string ngsUtilsVersion

set -e
set -u

module load "${HTSlibVersion}"
module load "${bedToolsVersion}"
module load "${ngsUtilsVersion}"
module list

## Filtering of the VCF files from the pipeline

#1 Compress  and indexing the VCFfile

project="/groups/umcg-gsad/tmp01/projects/GAP/NIST_TRIO"

for i in  $(ls "${project}/run01/results/vcf/"*".vcf")
do
	file=$(basename "${i}")
	mkdir -p "${project}/temp/"
	echo  "zipping vcf file : ${file} ..."
	bgzip -c "${i}" > "${project}/temp/${file}.gz"
	echo "indexing the vcf file of ${file}.gz ..."
		tabix -p vcf "${project}/temp/${file}.gz"

done

#2 Filtering the VCF using the bedfile

for i in $(ls "${project}/temp/"*".gz")
do
	echo "filtering the vcf file: ${i} ..."
	file=$(basename "${i}")
	sample=$(basename "${file}" ".FINAL.vcf.gz")
	bedtools intersect -a "${i}" -b "/home/umcg-molgenis/GAP/autoTestArray.bed" > "${project}/temp/${sample}.FINAL_FILTERED.vcf"
done

## Comparing the VCF files

for i in $(ls "${project}/temp/"*".FINAL_FILTERED.vcf")
do
	file=$(basename "${i}")
	echo "${file}"
	sample=$(basename "${file}" ".FINAL_FILTERED.vcf")
	echo "Comparing the VCF file with the TRUE VCF for the sample sample : ${sample} ..."
	export TERM=xterm-256color
	"${EBROOTNGSMINUTILS}/vcf-compare_2.0.sh" -1 "/home/umcg-molgenis/GAP/vcf/${sample}.FINAL_TRUE_FILTERED.vcf" -2 "${i}" -o "${project}/temp/VCF_Compare/${sample}/"

	## Checking if the output is correct

	vcfCheckvalue="$(awk 'NR == 2 {print $4}' ${project}/temp/VCF_Compare/${sample}/vcfStats.txt)"
	echo "${sample}: vcf Check value is: ${vcfCheckvalue}"

	if [ "${vcfCheckvalue}" = 100.00% ]
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

for i in $(ls "${project}/run01/results/PennCNV_reports/"*".txt")
do
	file=$(basename "${i}")
	sample=$(basename "${file}" ".txt")

	echo "Comparing the PennCNV per sample file for ${sample}"
	diffPennCNV="false"

	diff -q  <(tail -n +4 "/home/umcg-molgenis/GAP/PennCNV_reports/${sample}_TRUE.txt") <(tail -n +4 "${i}") || diffPennCNV="true"
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

for i in $(ls "${project}/run01/results/vcf/"*".sd")
do
	file=$(basename "${i}")
        sample=$(basename "${file}" ".vcf.sd")
	echo "${sample}"
        echo "Comparing the standard deviation's per sample for ${sample}"
        diffSD="false"
        diff -q "/home/umcg-molgenis/GAP/vcf/${sample}_TRUE.vcf.sd" "${i}" || diffPennCNV="true"
        if [ "${diffSD}" == "true" ]
        then
		echo "there are differences in the SD files between the test and original data for sample ${sample}"
                echo "please fix the bug or update this test"
                exit 1
        else
	echo "there are no differences between the SD files per sample, for sample ${sample}"
        fi
done


## Checking the callrate file
echo "Comparing the Project callrate files ... "

diffProjectCallrateFile="false"

diff -q "/home/umcg-molgenis/GAP/Callrates_NIST_TRIO_TRUE.txt" "${project}/run01/results/Callrates_NIST_TRIO.txt" || diffProjectCallrateFile="true"
        if [ "${diffProjectCallrateFile}" == "true" ]
        then
                echo "there are differences in the Project Callrate files between the test and original data for sample ${sample}"
                echo "please fix the bug or update this test"
                exit 1
        else
	echo "there are no differences between the Project Callrate files"
        fi
