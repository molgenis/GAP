#!/bin/bash

set -e
set -u

## USAGE ##

# bash GS_to_Opticall.sh \
# -i inputfile  \ (finalreport file)
# -o outputfolder

#Get command line options
while getopts ":i:o:" opt; do
  case "$opt" in
    i) input=$OPTARG ;;
    o) output=$OPTARG ;;
  esac
done
# Flag to extract header
flag=1
shift $(( OPTIND - 1 ))

SAMPLESHEET_SEP="\t"
declare -a sampleSheetColumnNames=()
declare -A sampleSheetColumnOffsets=()

headerNumber=$(( $(head -30 "${input}" |grep -n "\[Data\]" | grep -Eo '^[^:]+')+1))
echo ${headerNumber}

#get header, columnNames, and columnNumbers
IFS=$'\t' sampleSheetColumnNames=($(head -30 "${input}" | awk -v headerNumber="$headerNumber" '{if(NR==headerNumber){print $0}}'))
for (( offset = 0 ; offset < ${#sampleSheetColumnNames[@]:-0} ; offset++ ))
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
rm -f *tmp*
rm -f info
awk -v headerNumber="$headerNumber" '{if(NR>headerNumber){print >> $2".tmp"}}' "${input}"
# Extract info of all snps in the first condition + intensities of first sample
for a in *tmp;
    do
    sid=${a%.tmp}
    if [ "$flag" -eq 1 ]
        then
            less "$a" | awk -F "\t" -v name=$sid -v chr=$chr -v snp=$my_snp -v al=$alleles -v c=$coor -v A=$intA -v B=$intB 'BEGIN { print "Chr" "\t" "SNP" "\t" "Coor" "\t" "Alleles" "\t" name "A" "\t" name "B" }{ print $chr "\t" $snp "\t" $c "\t" substr($al,2,1) substr($al,4,1) "\t" $A "\t" $B}' > info
            flag=2
# keep intensities per sample
    else
        echo $sid
        less "$a" | awk -F "\t" -v name=$sid -v chr=$chr -v snp=$my_snp -v al=$alleles -v c=$coor -v A=$intA -v B=$intB 'BEGIN { print name "A" "\t" name "B" }{ print $A "\t" $B}' > "$a".tmp2
    fi
done
# Paste all intensities and paste snp info with the merged intensities
paste *tmp2 > final.tmp
paste info final.tmp > all_final.tmp
# Split by chromosome
awk '{if(NR>1){print >> $1".tmp3"}}' all_final.tmp

# Add header and remove chromosome column
for z in *tmp3;
    do
    sid2=${z%.tmp3}
    head -1 all_final.tmp >> "chr_tmp_"$sid2
    cat $z >> "chr_tmp_"$sid2
    cut -f2- "chr_tmp_"$sid2 > "chr_"$sid2
done
# Remove temp files.

cd ${output}

#for i in {0..22};do
#	rm _tmp_${i}*
#done

rm *tmp*
rm info
cd -
