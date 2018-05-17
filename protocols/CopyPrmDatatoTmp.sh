#list SentrixBarcode_A
#list SentrixPosition_A
#string GTCprmDataDir
#string GTCtmpDataDir


max_index=${#SentrixPosition_A[@]}-1

for i in ${SentrixBarcode_A[@]}
do
	mkdir -vp "${tmpDataDir}/${i}"
	for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
	do
		GTC_FILE="${GTCprmDataDir}/${i}/${i}_${SentrixPosition_A[samplenumber]}.gtc"
		rsync --verbose --recursive --links --no-perms --times --group --no-owner --devices --specials --checksum \
		"${GTC_FILE}"* \
		"${GTCtmpDataDir}/${i}/"
	done
done
