#MOLGENIS walltime=02:00:00 mem=4gb
#list SentrixBarcode_A,SentrixPosition_A
#string intermediateDir
#string resultDir
#string computeVersion
#string Project
#string projectJobsDir
#string projectRawTmpDataDir
#string genScripts
#string pipeline
#string runID
#string logsDir
#string perlVersion
#string group
#string gapVersion
#string workDir
#string workflowpath

umask 0007

module load ${computeVersion}
module list


#Create ProjectDirs
mkdir -p -m 2770 "${intermediateDir}"
mkdir -p -m 2770 "${resultDir}"
mkdir -p -m 2770 "${projectJobsDir}"
mkdir -p -m 2770 "${projectRawTmpDataDir}"

#Create Symlinks

rocketPoint=$(pwd)
host=$(hostname -s)

cd "${projectRawTmpDataDir}"

max_index=${#SentrixPosition_A[@]}-1

if [ "${pipeline}" == 'diagnostics' ] 
then
    for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
    do
        ln -sf "../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc" \
        "${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc"

        ln -sf "../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.md5" \
        "${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.md5"
    done
else
    for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
    do
        mkdir -p "${SentrixBarcode_A[samplenumber]}"
        ln -sf "../../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc" \
            "${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc"

        ln -sf "../../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.md5" \
            "${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.md5"
    done
fi

#Copying samplesheet to project jobs,results folder

cp "${genScripts}/${Project}.csv" "${projectJobsDir}/${Project}.csv"
cp "${genScripts}/${Project}.csv" "${resultDir}/${Project}.csv"

#
# Execute MOLGENIS/compute to create job scripts to analyse this project.
#

cd "${rocketPoint}"

perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${host}.csv" > "${rocketPoint}/parameters_host_converted.csv"
perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${group}.csv" > "${rocketPoint}/parameters_group_converted.csv"
perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${pipeline}.csv" > "${rocketPoint}/parameters_converted.csv"

sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_group_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${genScripts}/${Project}.csv" \
-p "${EBROOTGAP}/chromosomes.homo_sapiens.csv" \
-rundir "${projectJobsDir}" \
-w "${workflowpath}" \
--header "${EBROOTGAP}/templates/slurm/header.ftl" \
--submit "${EBROOTGAP}/templates/slurm/submit.ftl" \
--footer "${EBROOTGAP}/templates/slurm/footer.ftl" \
-b slurm \
-g \
-weave \
-runid "${runID}" \
-o "gapVersion=${gapVersion};\
runID=${runID}"


sampleSize=$(cat "${genScripts}/${Project}.csv" |  wc -l)

if [ "${pipeline}" == 'research' ] && [ "${sampleSize}" -gt 1000 ]
then
    echo "Samplesize is ${sampleSize}"
    ml "${perlVersion}"
    perl ${EBROOTGAP}/scripts/RemoveDuplicatesCompute.pl "${projectJobsDir}/"*"_mergeFinalReports_0.sh"
 fi
