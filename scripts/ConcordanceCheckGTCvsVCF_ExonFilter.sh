#!/bin/bash
#SBATCH --job-name=ConcordanceCheck
#SBATCH --output=ConcordanceCheck.out
#SBATCH --error=ConcordanceCheck.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 15gb
#SBATCH --nodes 1
#SBATCH --open-mode=append

set -e
set -u


PARSED_OPTIONS=$(getopt -n "$0"  -o w:a:n:o:t: --long "workDir:arrayVcfDir:ngsVcfDir:outputDir:tempDir:"  -- "$@")

#
# Bad arguments, something has gone wrong with the getopt command.
#
if [ $? -ne 0 ]; then
        usage
        echo "FATAL: Wrong arguments."
        exit 1
fi


eval set -- "$PARSED_OPTIONS"

#
# Now goes through all the options with a case and using shift to analyse 1 argument at a time.
# $1 identifies the first argument, and when we use shift we discard the first argument, so $2 becomes $1 and goes again through the case.
#
while true; do
    case "$1" in
    -w|--workDir)
                case "$2" in
                *) workDir=$2 ; shift 2 ;;
            esac ;;
    -a|--arrayVcfDir)
                case "$2" in
                *) arrayVcfDir=$2 ; shift 2 ;;
            esac ;;
    -n|--ngsVcfDir)
                case "$2" in
                *) ngsVcfDir=$2 ; shift 2 ;;
            esac ;;
    -o|--outputDir)
                case "$2" in
                *) outputDir=$2 ; shift 2 ;;
            esac ;;
    -t|--tempDir)
                case "$2" in
                *) tempDir=$2 ; shift 2 ;;
            esac ;;
         --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

empty=""
#
# Check required options were provided.
#
if [[ -z "${workDir-}" ]]; then
        workDir="/groups/umcg-gd/tmp05/Concordance/"
fi
if [[ -z "${arrayVcfDir-}" ]]; then
        arrayVcfDir="/groups/umcg-gd/tmp05/Concordance/array/"
fi
if [[ -z "${ngsVcfDir-}" ]]; then
        ngsVcfDir="/groups/umcg-gd/tmp05/Concordance/ngs/"
fi
if [[ -z "${outputDir-}" ]]; then
        outputDir="/groups/umcg-gd/tmp05/Concordance/output/"
fi
if [[ -z "${tempDir-}" ]]; then
        tempDir="/groups/umcg-gd/tmp05/Concordance/temp/"
fi

module load HTSlib
module load CompareGenotypeCalls/1.8.1-Java-1.8.0_74
module load BEDTools/2.25.0-foss-2015b

module list

for vcfFile in $(ls "${ngsVcfDir}"*"final.vcf")
do

    bedType="$(grep -m 1 -o -P 'intervals=\[[^\]]*.bed\]' "${vcfFile}" | cut -d [ -f2 | cut -d ] -f1)"
    echo "bedType: ${bedType}"
    bedDir="$(dirname ${bedType})"
    echo "bedDir: ${bedDir}"
    bedFile="${bedDir}/captured.merged.bed"
    echo "${bedFile}"

    ngsVcfId="$(basename "${vcfFile}" .final.vcf)"

    ##remove indel-calls from ngs-vcf
    declare -a patientList=($(grep '#CHROM' "${vcfFile}" | awk '{for(i=10;i<=NF;++i)print $i}'))
    declare -a dnaList=($(grep '#CHROM' "${vcfFile}" | awk '{for(i=10;i<=NF;++i)print $i}' | awk 'BEGIN {FS="_"}{if (NR>0){print substr($3,4)}}'))

    echo "${patientList[@]}"
    echo "${dnaList[@]}"

    ##find with DNA number the right NGS and array vcf file
    for patientNo in ${patientList[@]}
    do
        dnaNo="$(echo "${patientNo}" | awk 'BEGIN {FS="_"}{print substr($3,4)}')"
        declare -a arrayFile=($(ls -1 "${arrayVcfDir}/DNA-${dnaNo}_"*".FINAL.vcf"))

        echo "DNAno: ${dnaNo}"

    #################
    ## chech for duplo NGS vcf files.
    #################

    if [ -e "${ngsVcfDir}/archive/"*"${dnaNo}"* ]
    then
        echo "NGS sample duplo with DNAno:${dnaNo}, concordance is already calculated"
        mv "${vcfFile}" "${ngsVcfDir}/archive"
        continue
    fi        

    #################
    ## mail met notificatie als er meer dan 1 file is met hetzelfde DNA nummer
    #################

    if [[ "${#arrayFile[@]:-0}" != 1 ]]
    then
        echo "more than 1 file or file with DNAno:"${dnaNo}" does not exist "
        messageGCC="Dear GCC helpdesk,\n\nThere is more than 1 array sample with id:"${dnaNo}". \nPlease check if there is some think wrong with the concordance check.\nKind regards\nGCC"
        #printf '%b\n' "${messageGCC}" | mail -s "Concordance check error" 'helpdesk@somemailadres.com'
        continue  # exit 1
    fi

    mkdir -p "${tempDir}/${ngsVcfId}/"
    echo "${ngsVcfId}"

    grep '^#' "${vcfFile}" > "${tempDir}/${ngsVcfId}/${ngsVcfId}.FINAL.vcf"
    grep -v '^#' "${vcfFile}" | awk '{if (length($4)<2 && length($5)<2 ){print $0}}' >> "${tempDir}/${ngsVcfId}/${ngsVcfId}.FINAL.vcf"

    bgzip -c "${tempDir}/${ngsVcfId}/${ngsVcfId}.FINAL.vcf" > "${tempDir}/${ngsVcfId}/${ngsVcfId}.FINAL.vcf.gz"
    tabix -p vcf "${tempDir}/${ngsVcfId}/${ngsVcfId}.FINAL.vcf.gz"

        arrayId="$(basename "${arrayFile}" .FINAL.vcf)"
        echo "arrayID: ${arrayId}"
        touch "${tempDir}/${ngsVcfId}/${arrayId}.sampleId.txt"
        echo -e "data1Id\tdata2Id\n${arrayId}\t${patientNo}" >> "${tempDir}/${ngsVcfId}/${arrayId}.sampleId.txt"

        bedtools intersect -a "${arrayFile}" -b "${bedFile}" -header  > "${tempDir}/${ngsVcfId}/${arrayId}.FINAL.ExonFiltered.vcf"

        bgzip -c "${tempDir}/${ngsVcfId}/${arrayId}.FINAL.ExonFiltered.vcf" > "${tempDir}/${ngsVcfId}/${arrayId}.FINAL.ExonFiltered.vcf.gz"
        tabix -p vcf "${tempDir}/${ngsVcfId}/${arrayId}.FINAL.ExonFiltered.vcf.gz"

        echo "vcf-file:${vcfFile}"
        echo "ngs id :${ngsVcfId}"
        echo "vcf-file-no:${patientNo}"
        echo "arrayfile: ${arrayFile}"
        echo "arrayid: ${arrayId}"

        java -XX:ParallelGCThreads=1 -Djava.io.tmpdir="${tempDir}" -Xmx9g -jar ${EBROOTCOMPAREGENOTYPECALLS}/CompareGenotypeCalls.jar \
        -d1 "${tempDir}/${ngsVcfId}/${arrayId}.FINAL.ExonFiltered.vcf.gz" \
        -D1 VCF \
        -d2 "${tempDir}/${ngsVcfId}/${ngsVcfId}.FINAL.vcf.gz" \
        -D2 VCF \
        -ac \
        --sampleMap "${tempDir}/${ngsVcfId}/${arrayId}.sampleId.txt" \
        -o "${outputDir}/${arrayId}" \
        -sva

        mv "${vcfFile}" "${ngsVcfDir}/archive"
        mv "${arrayVcfDir}/${arrayId}.FINAL.vcf" "${arrayVcfDir}/archive"
        mv "${tempDir}/${ngsVcfId}/" "${tempDir}/archive/"

    done
done


