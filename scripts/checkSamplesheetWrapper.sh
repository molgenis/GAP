set -eu

for i in $(ls /groups/umcg-gap/scr01/Samplesheets/new/*.csv)
do
	mac2unix ${i}
	mv ${i} /groups/umcg-gap/scr01/Samplesheets/
done
