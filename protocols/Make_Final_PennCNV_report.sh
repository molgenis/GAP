#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#list Sample_ID
#list SentrixBarcode_A
#list SentrixPosition_A
#string intermediateDir
#string Project
#string logsDir
#string diagnosticOutputFolder
#string resultDir
#string FinalReportDir
#string PennCNV_reportDir

set -e
set -u


mkdir -p "${FinalReportDir}/"
makeTmpDir "${FinalReportDir}/"
tmpFinalReportDir="${MC_tmpFile}"

samples=()
count=0

barcodelist=()

n_elements=${Sample_ID[@]}
max_index=${#Sample_ID[@]}-1
for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
	barcodelist+=("${Sample_ID[samplenumber]}:${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}")
done


columnCount=9
for i in ${barcodelist[@]}
do
	echo "processing $i"
	echo "countboven: $count"
	sampleName=$(echo ${i} | awk 'BEGIN {FS=":"}{print $1}')
	barcodeCombined=$(echo ${i} | awk 'BEGIN {FS=":"}{print $2}')
	sampleName2=$(echo ${barcodelist[1]} | awk 'BEGIN {FS=":"}{print $1}')
	barcodeCombined2=$(echo ${barcodelist[1]} | awk 'BEGIN {FS=":"}{print $2}')

	if [ $count == 0 ]
	then
		echo "count=0"
				echo -en "Name\tChr\tPosition\t${sampleName}.GType\t${sampleName}.Log R Ratio\t${sampleName}.B Allele Freq\t${sampleName2}.GType\t${sampleName2}.Log R Ratio\t${sampleName2}.B Allele Freq" > ${tmpFinalReportDir}/header.txt
		echo "join -1 2 -2 2 -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,2.5,2.6,2.7 ${intermediateDir}/${barcodeCombined} ${intermediateDir}/${barcodeCombined2} > ${tmpFinalReportDir}/bron.txt"
		join -1 2 -2 2 -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,2.5,2.6,2.7  "${PennCNV_reportDir}/${barcodeCombined}.gtc.txt" "${PennCNV_reportDir}/${barcodeCombined2}.gtc.txt" > "${tmpFinalReportDir}/bron.txt"
	elif [[ $count -gt 1 ]]
	then
		echo "count>1"

		format=''
		for (( col = 1 ; col <= ${columnCount:-0} ; col++ ))
		do
			format="${format}1.${col},"
		done
		format="${format}2.5,2.6,2.7"
		echo ${format}
		echo ${i}
		echo "join -1 1 -2 2 -t $'\t' -o ${format}  ${tmpFinalReportDir}/bron.txt   ${PennCNV_reportDir}/${barcodeCombined}.gtc.txt >> ${tmpFinalReportDir}/output.txt"

		join -1 1 -2 2 -t $'\t' -o "${format}" "${tmpFinalReportDir}/bron.txt" "${PennCNV_reportDir}/${barcodeCombined}.gtc.txt" >> "${tmpFinalReportDir}/output.txt"
		columnCount=$columnCount+3
		mv "${tmpFinalReportDir}/output.txt" "${tmpFinalReportDir}/bron.txt"

		echo -en "\t${sampleName}.GType\t${sampleName}.Log R Ratio\t${sampleName}.B Allele Freq" >> ${tmpFinalReportDir}/header.txt

	fi
	echo "countonder: $count"
	count=$((count+1))
done

(cat "${tmpFinalReportDir}/header.txt"; printf "\n"; cat "${tmpFinalReportDir}/bron.txt") > "${tmpFinalReportDir}/${Project}_PennCNV.txt"

#Move PennCNV report to Intermediate dir
echo "mv ${tmpFinalReportDir}/${Project}_PennCNV.txt ${FinalReportDir}/${Project}_PennCNV.txt"
mv "${tmpFinalReportDir}/${Project}_PennCNV.txt" "${FinalReportDir}/${Project}_PennCNV.txt"


#Copy results to resultDir

rsync -a "${FinalReportDir}/${Project}_PennCNV.txt" "${resultDir}"
rsync -a "${FinalReportDir}/${Project}_PennCNV.txt" "${diagnosticOutputFolder}/${Project}/"

#Touch file for DARWIN so they know pipeline is finished and can start proceeding step to put data in SNP Module...
touch "${diagnosticOutputFolder}/${Project}/${Project}.finished"
