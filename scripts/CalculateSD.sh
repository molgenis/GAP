#!/bin/bash

set -e
set -u

input_dir="/groups/umcg-gap/tmp05/tmp/GlobalScreeningArray-24+v1.0_000530/"*
output_dir="/groups/umcg-gd/tmp05/Concordance/array/"

for input_file in ${input_dir}
do
        log_ratios=$(awk '{if ($3 != "X" && $3 != "Y" && $3 != "XY" ) print $6}' "${input_file}")
        sd=$(echo "${log_ratios[@]}" | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

        echo "$(basename ${input_file}): ${sd}"

	if [[ "${sd}" < 0.2 ]]
	then
		echo "mv ${input_file} ${output_dir}"
		mv "${input_file}" "${output_dir}"
	fi
done
