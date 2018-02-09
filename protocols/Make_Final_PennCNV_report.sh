#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#list Sample_ID
#list SentrixBarcode_A
#list SentrixPosition_A
#string intermediateDir
#string Project

set -e
set -u

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
                echo -en "Name\tChr\tPosition\t${sampleName}.GType\t${sampleName}.Log R Ratio\t${sampleName}.B Allele Freq\t${sampleName2}.GType\t${sampleName2}.Log R Ratio\t${sampleName2}.B Allele Freq" > ${intermediateDir}/header.txt
		echo "join -1 2 -2 2 -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,2.5,2.6,2.7 ${intermediateDir}/${barcodeCombined} ${intermediateDir}/${barcodeCombined2} > ${intermediateDir}/bron.txt"
		join -1 2 -2 2 -t $'\t' -o 1.2,1.3,1.4,1.5,1.6,1.7,2.5,2.6,2.7  "${intermediateDir}/${barcodeCombined}.gtc.txt" "${intermediateDir}/${barcodeCombined2}.gtc.txt" > "${intermediateDir}/bron.txt"
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
		echo "join -1 1 -2 2 -t $'\t' -o ${format}  ${intermediateDir}/bron.txt   ${intermediateDir}/${barcodeCombined}.gtc.txt >> ${intermediateDir}/output.txt"

		join -1 1 -2 2 -t $'\t' -o "${format}" "${intermediateDir}/bron.txt" "${intermediateDir}/${barcodeCombined}.gtc.txt" >> "${intermediateDir}/output.txt"
		columnCount=$columnCount+3
		mv "${intermediateDir}/output.txt" "${intermediateDir}/bron.txt"

		echo -en "\t${sampleName}.GType\t${sampleName}.Log R Ratio\t${sampleName}.B Allele Freq" >> ${intermediateDir}/header.txt

	fi
	echo "countonder: $count"
	count=$((count+1))
done

(cat "${intermediateDir}/header.txt"; printf "\n"; cat "${intermediateDir}/bron.txt") > "${intermediateDir}/${Project}_PennCNV.txt"

