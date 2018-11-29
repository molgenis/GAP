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


PARSED_OPTIONS=$(getopt -n "$0"  -o w:a:n:o:t: --long "workdir:arrayvcfdir:ngsvcfdir:outputdir:tempdir:prmdir:"  -- "$@")

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
    -w|--workdir)
                case "$2" in
                *) workdir=$2 ; shift 2 ;;
            esac ;;
    -a|--arrayvcfdir)
                case "$2" in
                *) arrayvcfdir=$2 ; shift 2 ;;
            esac ;;
    -n|--ngsvcfdir)
                case "$2" in
                *) ngsvcfdir=$2 ; shift 2 ;;
            esac ;;
    -o|--outputdir)
                case "$2" in
                *) outputdir=$2 ; shift 2 ;;
            esac ;;
    -t|--tempdir)
                case "$2" in
                *) tempdir=$2 ; shift 2 ;;
            esac ;;
    -p|--prmdir)
                case "$2" in
                *) tempdir=$2 ; shift 2 ;;
            esac ;;
         --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done


empty=""
#
# Check required options were provided.
#
if [[ -z "${workdir-}" ]]; then
        workdir="/groups/umcg-gd/tmp05/Concordance/"
fi
if [[ -z "${arrayvcfdir-}" ]]; then
        arrayvcfdir="/groups/umcg-gd/tmp05/Concordance/array/"
fi
if [[ -z "${ngsvcfdir-}" ]]; then
        ngsvcfdir="/groups/umcg-gd/tmp05/Concordance/ngs/"
fi
if [[ -z "${outputdir-}" ]]; then
        outputdir="/groups/umcg-gd/tmp05/Concordance/output/"
fi
if [[ -z "${tempdir-}" ]]; then
        tempdir="/groups/umcg-gd/tmp05/Concordance/temp/"
fi
if [[ -z "${prmdir-}" ]]; then
        tempdir="/groups/umcg-gd/prm03/Concordance/"
fi

module load HTSlib
module load CompareGenotypeCalls/1.8.1-Java-1.8.0_74
module load BEDTools/2.25.0-foss-2015b

module list

for vcffile in $(ls "${ngsvcfdir}"*"final.vcf")
do
    bedtype=($(grep -m 1 -o -P 'intervals=\[[^\]]*.bed\]' ${vcffile} | cut -d [ -f2 | cut -d ] -f1))
    echo "bedtype: ${bedtype}"
    beddir=($(dirname ${bedtype}))
    echo "beddir: ${beddir}"
    bedfile="${beddir}/captured.merged.bed"
    echo "${bedfile}"

    ngsvcfid=$(basename "${vcffile}" .final.vcf)

    grep '^#' ${vcffile} > ${ngsvcfid}.FINAL.vcf
    grep -v '^#' ${vcffile} | awk '{if (length($4)<2 && length($5)<2 ){print $0}}' >> ${ngsvcfid}.FINAL.vcf

    bgzip -c ${ngsvcfid}.FINAL.vcf > ${ngsvcfid}.FINAL.vcf.gz
    tabix -p vcf ${ngsvcfid}.FINAL.vcf.gz

    patientlist=()
    dnalist=()

    ##remove indel-calls from ngs-vcf
    patientlist+=($(grep '#CHROM' ${vcffile} | awk '{for(i=10;i<=NF;++i)print $i}'))
    dnalist+=($(grep '#CHROM' ${vcffile} | awk '{for(i=10;i<=NF;++i)print $i}' | awk ' {FS="_"}{if (NR>1){print substr($3,4)}}')) 

        ##find with DNA number the right NGS and array vcf file
    for patientno in ${patientlist[@]}
    do
	dnano=$(echo ${patientno} | awk 'BEGIN {FS="_"}{print substr($3,4)}' )
        declare -a arrayfile=($(ls -1 "${arrayvcfdir}/DNA-${dnano}_"*".FINAL.vcf" 2> /dev/null ))
        echo "DNAno: ${dnano}"
        echo "arrayfile: ${arrayfile}"
            ##################
            ### mail met notificatie als er meer dan 1 file is met hetzelfde DNA nummer
            ##################

	if [[ ${#arrayfile[@]:-0} != 1 ]]
        then
            echo "more than 1 file or file ${dnancdo} does not exist "
            messageGCC="Dear GCC helpdesk,\n\nThere is more than 1 array sample with id:${dnano}. \nPlease check if there is some think wrong with the concordance check.\nKind regards\nGCC"
            #printf '%b\n' "${messageGCC}" | mail -s "Concordance check error" 'helpdesk.gcc.groningen@gmail.com'
            next # exit 1
        fi

        arrayid=($(basename ${arrayfile} .FINAL.vcf ))
        echo "arrayID: ${arrayid}"
        touch ${tempdir}/${arrayid}.sampleId.txt
        echo -e "data1Id\tdata2Id\n${arrayid}\t${patientno}" >> ${tempdir}/${arrayid}.sampleId.txt

        bedtools intersect -a "${arrayfile}" -b "${bedfile}" -header  >"${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf"

        bgzip -c "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf" > "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf.gz"
        tabix -p vcf "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf.gz"

        echo "vcf-file:${vcffile}"
        echo "ngs id :${ngsvcfid}"
        echo "vcf-file-no:${patientno}"
        echo "arrayfile: ${arrayfile}"
        echo "arrayid: ${arrayid}"

        java -XX:ParallelGCThreads=1 -Djava.io.tmpdir="${tempdir}" -Xmx9g -jar ${EBROOTCOMPAREGENOTYPECALLS}/CompareGenotypeCalls.jar \
        -d1 "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf.gz" \
        -D1 VCF \
        -d2 "${ngsvcfid}.FINAL.vcf.gz" \
        -D2 VCF \
        -ac \
        --sampleMap "${tempdir}/${arrayid}.sampleId.txt" \
        -o "${outputdir}/${arrayid}" \
        -sva

	mv "${vcffile}" "${ngsvcfdir}/archive"
        mv "${ngsvcfid}.FINAL.vcf" "${ngsvcfdir}/archive"
        mv "${ngsvcfid}.FINAL.vcf.gz" "${ngsvcfdir}/archive"
        mv "${ngsvcfid}.FINAL.vcf.gz.tbi" "${ngsvcfdir}/archive"
        mv "${arrayvcfdir}/${arrayid}.vcf" "${arrayvcfdir}/archive"
        mv "${arrayvcfdir}/${arrayid}.FINAL.vcf" "${arrayvcfdir}/archive"
        mv "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf" "${arrayvcfdir}/archive"
        mv "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf.gz" "${arrayvcfdir}/archive"
        mv "${arrayvcfdir}/${arrayid}.FINAL.ExonFiltered.vcf.gz.tbi" "${arrayvcfdir}/archive"
        mv "${tempdir}/${arrayid}.sampleId.txt" "${tempdir}/archive"

	#rsync final results to PRM
#	rsync -av  "" "umcg-gd-dm@boxy.hpc.rug.nl:${}"
	# ${DATA_MANAGER}@${sourceServerFQDN}:${SCR_ROOT_DIR}/Samplesheets/${_run}.${SAMPLESHEET_EXT}
	done
done


