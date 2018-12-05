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


PARSED_OPTIONS=$(getopt -n "$0"  -o t:a: --long "tmpDir:archiveTime"  -- "$@")

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
    -a|--archiveTime)
                case "$2" in
                *) archiveTime=$2 ; shift 2 ;;
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

#Checking if VCF File array is older then ArchiveTime and archive file if this is the case

for i in $(find "${tmpDir}/array/" -type f -maxdepth 1 -mtime +"${archiveTime}" -print)
do
echo "${i}"
mv "${i}" "${tmpDir}/array/archive/"
done

#Checking if VCF File Array is older then ArchiveTime and archive file if this is the case

for i in $(find "${tmpDir}/ngs/" -type f -maxdepth 1 -mtime +"${archiveTime}" -print)
do
echo "${i}"
mv "${i}" "${tmpDir}/ngs/archive/"
done
