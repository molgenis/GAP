#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string intermediateDir
#string resultDir
#string projectJobsDir
#string Project
#string logsDir
#string runID
#string finalReport
#string samplesheet
#string optiCallDir
#string genSampleDir

set -e
set -u

if [ ! -d "${logsDir}/${Project}/" ]
then
	mkdir -p "${logsDir}/${Project}/"
fi

#Copying  outputfiles to resultsDir

printf "Copied ${finalReport} file to project results directory.."
rsync -a "${finalReport}" "${resultDir}"
printf ".. finished (1/4)\n"
printf "Copied ${samplesheet} file to project results directory.."
rsync -a "${samplesheet}" "${resultDir}"
printf ".. finished (2/4)\n"
printf "Copied ${optiCallDir} dir to project results directory.."
rsync -a "${optiCallDir}" "${resultDir}"
printf ".. finished (3/4)\n"
printf "Copied ${genSampleDir} dir to project results directory.."
rsync -a "${genSampleDir}" "${resultDir}"
printf ".. finished (4/4)\n"

# Touch log file for GAP_Automated for starting copying project data to PRM

if [ -f "${logsDir}//${Project}/${runID}.pipeline.started" ]
then
	mv "${logsDir}/${Project}/${runID}.pipeline".{started,finished}
else
	touch "${logsDir}/${Project}/${runID}.pipeline.finished"
fi
rm -f "${logsDir}/${Project}/${runID}.pipeline.failed"
echo "${logsDir}/${Project}/${runID}.pipeline.finished is created"

touch "${logsDir}/${Project}/pipeline.finished"
