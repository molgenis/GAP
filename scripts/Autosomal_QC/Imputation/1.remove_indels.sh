#################################
### remove InDels ##
### 10/09/2019
### v1.0
### author: ealopera
#################################

ml plink

mkdir -p 1.InDels_removed
for chr in {1..22} "XY" "X"
do 

awk '$5=="D" ||$5=="I"' chr_${chr}.bim|cut -f2>>InDels.remove

plink --bfile chr_$chr \
--exclude InDels.remove \
--make-bed \
--out 1.InDels_removed/chr_$chr 

done





