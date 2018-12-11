#!/bin/bash
#SBATCH --job-name=Archive_Concordance_results
#SBATCH --output=Archive_Concordance_results.out
#SBATCH --error=Archive_Concordance_results.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 15gb
#SBATCH --nodes 1
#SBATCH --open-mode=append

set -e
set -u


PARSED_OPTIONS=$(getopt -n "$0"  -o t:p:a:d:s: --long "tmpDir:prmDir:archiveTime:dataManager:sourceServer"  -- "$@")

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
    -t|--tmpDir)
                case "$2" in
                *) tmpDir=$2 ; shift 2 ;;
            esac ;;
    -p|--prmDir)
                case "$2" in
                *) prmDir=$2 ; shift 2 ;;
            esac ;;
    -a|--archiveTime)
                case "$2" in
                *) archiveTime=$2 ; shift 2 ;;
            esac ;;
    -d|--dataManager)
                case "$2" in
                *) dataManager=$2 ; shift 2 ;;
            esac ;;
    -s|--sourceServer)
                case "$2" in
                *) sourceServer=$2 ; shift 2 ;;
       esac ;;
         --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

empty=""
#
# Check required options were provided.
#
if [[ -z "${tmpDir-}" ]]; then
        tmpDir="/groups/umcg-gd/tmp05/Concordance/"
fi
if [[ -z "${prmDir-}" ]]; then
        prmDir="/groups/umcg-gd/prm03/Concordance/"
fi
if [[ -z "${archiveTime-}" ]]; then
        archiveTime="31"
fi
if [[ -z "${dataManager-}" ]]; then
        dataManager="umcg-gap-dm"
fi
if [[ -z "${sourceServer-}" ]]; then
        sourceServer="zinc-finger.gcc.rug.nl"
fi

#Copying output to PRM
for resultsFile in $(ssh ${dataManager}@${sourceServer} ls ${tmpDir}/output/'*'sample)
do
	file="$(basename ${resultsFile})"
	sampleID="$(echo "${file}" | awk 'BEGIN {FS="."}{print $1}')"
	echo ${resultsFile}
	echo ${file}
	if [[ ! -f "${prmDir}/output/${sampleID}.sample" || ! -f "${prmDir}/output/${sampleID}.variants" ]]
	then
		#rsyncing results to PRM
		echo "rsync -av ${dataManager}@${sourceServer}:${resultsFile} ${prmDir}/output/"
		rsync -av "${dataManager}@${sourceServer}:${resultsFile}" "${prmDir}/output/"

		echo "rsync -av ${dataManager}@${sourceServer}:${tmpDir}/output/${sampleID}.variants ${prmDir}/output/"
		rsync -av "${dataManager}@${sourceServer}:${tmpDir}/output/${sampleID}.variants" "${prmDir}/output/"
	else
		echo "there are no new results to copy to prm"
	fi
done
