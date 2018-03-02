#list SentrixBarcode_A
#list SentrixPosition_A
#string prmDataDir
#string tmpDataDir


max_index=${#SentrixPosition_A[@]}-1

for i in ${SentrixBarcode_A[@]}
do
	mkdir -vp "${tmpDataDir}/${i}"
	for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
	do
		GTC_FILE="${prmDataDir}/${i}/${i}_${SentrixPosition_A[samplenumber]}.gtc"
		rsync --verbose --recursive --links --no-perms --times --group --no-owner --devices --specials --checksum \
		"${GTC_FILE}"* \
		"${tmpDataDir}/${i}/"
	done
done
