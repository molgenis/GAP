#!/bin/bash

set -e
set -u
set -o pipefail

#
# Disable some shellcheck warnings:
#  * SC2004: $/${} is unnecessary on arithmetic variables.
#            But for consistency we prefer to always use ${} anyway.
#  * SC2015: Note that A && B || C is not if-then-else. C may run when A is true.
#            We know and use this construct regularly to create "transactions"
#            where C is only executed when both A and B have succeeded.
#  * SC2148: the shebang is declared in the header.ftl not in the protocols
#
export SHELLCHECK_OPTS="-e SC2004 -e SC2015 -e SC2148"

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script for sanity checking of Bash code of this repo using ShellCheck.
Usage:
	$(basename "${0}") OPTIONS
Options:
	-h	Show this help.
	-v	Enables verbose output.
===============================================================================================================
EOH
	exit 0
}


#
# Parse commandline options
#
declare format='gcc'  # default
while getopts ":hv" opt
do
	case "${opt}" in
		h)
			showHelp
			;;
		v)
			verbose='1'
			;;
		\?)
			printf '%s\n' "FATAL: Invalid option -${OPTARG}. Try $(basename "${0}") -h for help."
			exit 1
			;;
		:)
			printf '%s\n' "FATAL: Option -${OPTARG} requires an argument. Try $(basename "${0}") -h for help."
			exit 1
			;;
		*)
			printf '%s\n' "FATAL Unhandled option. Try $(basename "${0}") -h for help."
			exit 1
			;;
	esac
done

#
# Check if ShellCheck is installed.
#
which shellcheck 2>&1 >/dev/null \
	|| {
		printf '%s\n' 'FATAL: cannot find shellcheck; make sure it is installed and can be found in ${PATH}.'
		exit 1
	}

MYDIR="$(cd -P "$(dirname "${0}")" && pwd)"
#
# Run ShellCheck for all Bash scripts in the protocols/ subdir.
#
if [[ "${CIRCLECI:-}" == true ]]
then
	#
	# Exclude SC2154 (warning for variables that are referenced but not assigned),
	# because we cannot easily resolve variables filled in using templating with Molgenis Compute.
	#
	export SHELLCHECK_OPTS="${SHELLCHECK_OPTS} -e SC2154"
	#
	# Exclude SC2312 (warning for masking return values of command in a subshell when using process substitution) temporarily,
	# because we did not find a reliable fix yet....
	#
	#export SHELLCHECK_OPTS="${SHELLCHECK_OPTS} -e SC2312"
fi
shellcheck -a -x -o all -f "${format}" "${MYDIR}"/../protocols/*.sh | sed "s|${MYDIR}/../||g"