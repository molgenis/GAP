#!/bin/bash

set -e
set -u

callsFile="${1}"
probsFile="${2}"
outputFolder="${3}"
fileName="${callsFile%.calls}"
chrNr="${fileName#*_}"

echo "$callsFile $probsFile output=$outputFolder $fileName $chrNr"
#exit 0

if [[ -z "$callsFile" ]]
then
	echo "Error: set input"
	exit
fi

if [[ -z "$outputFolder" ]]
then
	echo "Error: set output"
	exit
fi

	echo "${chrNr}"

	input="${probsFile}"
	output="chr_${chrNr}"

	sampleFile="${fileName}.sample"

	if [ -e $input ]
	then

		awk -v chr=${chrNr} -v sampleFile="${sampleFile}" '

		NR == 1 {
			print "ID_1 ID_2 missing" > sampleFile
			print "0 0 0" > sampleFile
			sampleCount = 0
			for(i=5;i<=NF;i+=1){
				print $i,$i,0 > sampleFile
				sampleCount += 1
			}
		}

		NR > 1 {

			ORS="";
			print chr, $1, $2, substr($3,1,1), substr($3,2,1);
			for(i=0;i<sampleCount;i+=1){
				x = 4 + i * 4
				print $x,$(x+2),$(x+1)
			}
			ORS="";
			print "\n"
		}
		' < ${input} > ${output}.gen
fi

#echo "mv ${tmpOutputFolder}/chr_${chrNr}.gen ${outputFolder}/"
#echo "mv ${tmpOutputFolder}/chr_${chrNr}.sample ${outputFolder}/"

#mv "${tmpOutputFolder}/chr_${chrNr}.gen" "${outputFolder}/"
#mv "${tmpOutputFolder}/chr_${chrNr}.sample" "${outputFolder}/"
