#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string intermediateDir
#string Project
#string logsDir
#string runID

set -e
set -u

# Touch log file for GAP_Automated for starting copying project data to PRM

if [ -f "${logsDir}//${Project}/${runID}.pipeline.started" ]
then
	mv "${logsDir}/${Project}/${runID}.pipeline".{started,finished}
else
	touch "${logsDir}/${Project}/${runID}.pipeline.finished"
fi
rm -f "${logsDir}/${Project}/${runID}.pipeline.failed"
echo "${logsDir}/${Project}/${runID}.pipeline.finished is created"

if [ ! -d "${logsDir}/${Project}/" ]
then
	mkdir -p "${logsDir}/${Project}/"
fi

touch pipeline.finished
