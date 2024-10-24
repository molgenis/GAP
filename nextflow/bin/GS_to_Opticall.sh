#!/bin/bash

set -e
set -u

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Usage:
	$(basename $0) OPTIONS
Options:
	-h	Show this help.
	Required:
	-i inputfile   (final report file)
	-r recalculateControls
	-c includeControls
	-o outputDir
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

#Get command line options
while getopts "i:o:c:r:h" opt; do
	case "${opt}" in
		i) input=$OPTARG ;;
		r) recalculate=$OPTARG ;;
		c) includeControls=$OPTARG ;;
		o) output=$OPTARG ;;
		h) showHelp;;
	esac
done

if [[ -z "${input:-}" ]]
then
	echo -e '\nERROR: -i input GS final report file not specified.\n'
	showHelp
	exit 1
fi
if [[ -z "${output:-}" ]]
then
	echo -e "\nERROR: -o output dir not specified.\n"
	showHelp
	exit 1
fi
if [[ -z "${recalculate:-}" ]]
then
	recalculate="no"
fi
if [[ -z "${includeControls:-}" ]]
then
	includeControls="no"
fi


# Flag to extract header
flag=1
shift $(( OPTIND - 1 ))

SAMPLESHEET_SEP="\t"
declare -a sampleSheetColumnNames=()
declare -A sampleSheetColumnOffsets=()

headerNumber=$(( $(head -30 "${input}" |grep -n "\[Data\]" | grep -Eo '^[^:]+')+1))
echo "${headerNumber}"

#get header, columnNames, and columnNumbers
IFS=$'\t' sampleSheetColumnNames=($(head -30 "${input}" | awk -v headerNumber="${headerNumber}" '{if(NR==headerNumber){print $0}}'))
for (( offset = 0 ; offset < ${#sampleSheetColumnNames[@]} ; offset++ ))
do
	columnName="${sampleSheetColumnNames[${offset}]}"
	sampleSheetColumnOffsets["${columnName}"]="${offset}"
done

# map ColumnHeader to correct variable
my_snp=$(( ${sampleSheetColumnOffsets['SNP Name']} +1 ))
alleles=$(( ${sampleSheetColumnOffsets['SNP']} +1 ))
coor=$(( ${sampleSheetColumnOffsets['MapInfo']} +1 ))
#coor=$(( ${sampleSheetColumnOffsets['Position']} +1 ))
chr=$(( ${sampleSheetColumnOffsets['Chr']} +1 ))
intA=$(( ${sampleSheetColumnOffsets['X']} +1 ))
intB=$(( ${sampleSheetColumnOffsets['Y']} +1 ))

headerNumber=$(( $(head -30 "${input}" |grep -n "\[Data\]" | grep -Eo '^[^:]+')+1))

# Split file per sample id, in this case: second column (first 10 columns is a header)
cd "${output}"
rm -f *"tmp"*
rm -f "info"
awk -v headerNumber="${headerNumber}" '{if(NR>headerNumber){print >> $2".tmp"}}' "${input}"
# Extract info of all snps in the first condition + intensities of first sample
	for a in *tmp;
	do
	sid="${a%.tmp}"
	if [[ "${flag}" -eq 1 ]]
	then
		less "${a}" | awk -F "\t" -v name="${sid}" -v chr="${chr}" -v snp="${my_snp}" -v al="${alleles}" -v c="${coor}" -v A="${intA}" -v B="${intB}" 'BEGIN { print "Chr" "\t" "SNP" "\t" "Coor" "\t" "Alleles" "\t" name "A" "\t" name "B" }{ print $chr "\t" $snp "\t" $c "\t" substr($al,2,1) substr($al,4,1) "\t" $A "\t" $B}' > "info"
		flag=2
	# keep intensities per sample
	else
		echo "${sid}"
		less "${a}" | awk -F "\t" -v name="${sid}" -v chr="${chr}" -v snp="${my_snp}" -v al="${alleles}" -v c="${coor}" -v A="${intA}" -v B="${intB}" 'BEGIN { print name "A" "\t" name "B" }{ print $A "\t" $B}' > "${a}.tmp2"
	fi
	done
	
	## If there was no recalculation done we need to include the .tmp2 files from the GDIO controls
	if [[ "${includeControls}" == 'yes' ]]
	then
		if [[ "${recalculate}" == 'no' ]]
		then
			echo "no need for a recalculation, start to unpack the ones from: /apps/data/PGx/gdio_controls.tar.gz in ${output}"
			tar xvzf '/apps/data/PGx/gdio_controls.tar.gz' -C "${output}"
		elif [[ "${recalculate}" == 'yes' ]]
		then
			echo "everything in the samplesheet will be recalculated"
		else
			echo "recalculate is nothing? [${recalculate}]"
			exit 1
		fi
	else
		echo "no Controls included"
	fi
		
	
	ls -1 *.tmp2 | split -l 500 -d - tmp2list

	for list in tmp2list*
	do 
	cat << EOH > "${output}/job_${list}.sh"
#!/bin/bash
#SBATCH --job-name=list_${list}
#SBATCH --output=${output}/job_${list}.out
#SBATCH --error=${output}/job_${list}.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 8gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=60L

set -o  pipefail
set -eu

paste \$(cat ${output}/$list) > ${output}/merge${list##lists} 

touch ${output}/${list}.list.finished

EOH

	echo "$list done"
	sbatch "${output}/job_${list}.sh"

	done	
	
	sizeOfList=$(ls -1 tmp2list* | wc -l)
	sizeOfList=$((sizeOfList +1))
	count=0
	touch "${output}/fake.list.finished"
	finished=$(ls -1 ${output}/*.list.finished | wc -l)
	while [[ "${sizeOfList}" -ne  "${finished}" ]]
	do
		finished=$(ls -1 ${output}/*.list.finished | wc -l)
		sleep 120
		count=$((count+2))
		if [[ "${count}" -eq '30' ]]
		then
			echo "not finished in time"
			exit 1
		fi
		
	done
	echo "done"
	
# Paste all intensities and paste snp info with the merged intensities
paste "info" mergetmp2list* > "all_final.tmp"
# Split by chromosome
awk '{if(NR>1){print >> $1".tmp3"}}' "all_final.tmp"

# Add header and remove chromosome column
for z in *"tmp3";
	do
	sid2="${z%.tmp3}"
	head -1 "all_final.tmp" >> "chr_tmp_${sid2}"
	cat "${z}" >> "chr_tmp_${sid2}"
	cut -f2- "chr_tmp_${sid2}" > "chr_${sid2}"
	done
# Remove temp files.

cd "${output}"

#for i in {0..22};do
#	rm _tmp_${i}*
#done

rm *tmp* 
rm *list*
rm info
cd -
