#MOLGENIS walltime=00:30:00 mem=40gb ppn=1

#string Project

#string genSampleDir
#string AutosomeQCDir
#string plinkVersion
#string gapVersion
#string output80
#string output95
#string outputMH
#string repout
#string RPlusVersion
#string MAFref
#string logsDir
#string intermediateDir
#strinf samplesheet


module load "${plinkVersion}"
module load "${RPlusVersion}"
module load "${gapVersion}"
module list


output80="${AutosomeQCDir}/1_CR80" 
output95="${AutosomeQCDir}/2_CR95"
outputMH="${AutosomeQCDir}/3_MAF_HWE"
Het="${AutosomeQCDir}/4_Het"
Relatedness="${AutosomeQCDir}/5_Relatedness"



mkdir -p "${AutosomeQCDir}/"
mkdir -p "${output80}/"
mkdir -p "${output95}/"
mkdir -p "${outputMH}/"
mkdir -p "${repout}/"
mkdir -p "${Het}/"
mkdir -p "${Het}/proc/"
mkdir -p "${Relatedness}/"
mkdir -p "${Relatedness}/proc/"




for chr in {1..22} "XY"
do

 ### create plink files and call_rate stats for individuales and SNPs
  plink --data ${genSampleDir}/chr_${chr} \
        --make-bed  \
        --missing \
        --out ${AutosomeQCDir}/chr_${chr}

  ##create list of SNPs  and vars to exclude on the criteria callrate<=80
  awk '$5>0.20 {print $2}' ${AutosomeQCDir}/chr_${chr}.lmiss > ${output80}/chr_${chr}.extr80_var
  awk '$6>0.20 {print $1, $2}' ${AutosomeQCDir}/chr_${chr}.imiss > ${output80}/chr_${chr}.extr80_sam
  
done

##creates list of individuals and snps excluded for all teh autosomes
cat ${output80}/chr_*.extr80_sam |sort -u > ${output80}/extr80.samples
cat ${output80}/chr_*.extr80_var > ${output80}/extr80.vars



for chr in {1..22} "XY"
do

  # exclude individuals with callrate<=80 (creates excluded individuals_file) creates data set with  SNP_ and individual_ callrate>80
  plink --bfile ${AutosomeQCDir}/chr_${chr}  \
        --make-bed \
        --remove ${output80}/extr80.samples \
        --exclude ${output80}/chr_${chr}.extr80_var \
      --out ${output80}/chr_${chr}


  #calculates callrate stats from the previously filtered datafile
  plink --bfile ${output80}/chr_${chr} \
         --missing \
         --out ${output80}/chr_${chr}


  ##create list of SNPs snd samples to exclude on the criteria callrate<=95
  awk '$5>0.05 {print $2}' ${output80}/chr_${chr}.lmiss > ${output95}/chr_${chr}.extr95_var
  awk '$6>0.05 {print $1, $2}' ${output80}/chr_${chr}.imiss > ${output95}/chr_${chr}.extr95_sam

done

cat ${output95}/chr_*.extr95_sam|sort -u > ${output95}/extr95.samples
cat ${output95}/chr_*.extr95_var > ${output95}/extr95.vars


for chr in {1..22} "XY"
do

  ## exclude individuals with callrate<=95 (creates excluded individuals_file) creates data set with  SNP_ and individual_ callrate>95
  plink --bfile ${output80}/chr_${chr}  \
        --make-bed \
        --remove ${output95}/extr95.samples \
        --exclude ${output95}/chr_${chr}.extr95_var \
        --out ${output95}/chr_${chr}

  ## calculate MAF and HWE 
  plink --bfile ${output95}/chr_${chr} \
        --freq \
        --hardy \
        --out ${outputMH}/chr_${chr}


done


##creates merged files of included individuals and SNPs to be used of further analysis

cat ${AutosomeQCDir}/chr_*.fam|sort -u|awk '{print$2}' > ${AutosomeQCDir}/full.ind
cat ${AutosomeQCDir}/chr_*.bim|awk '{print$2}' > ${AutosomeQCDir}/full.snps

cat ${output80}/chr_*.fam|sort -u|awk '{print$2}' > ${output80}/incl80.samples
cat ${output80}/chr_*.bim|awk '{print$2}' > ${output80}/incl80.vars

cat ${output95}/chr_*.fam|sort -u|awk '{print$2}' > ${output95}/incl95.samples
cat ${output95}/chr_*.bim|awk '{print$2}' > ${output95}/incl95.vars


#### heterozygosity analysis by sample

for chr in {1..22} "XY" 
do
   ##first, separate a list of SNPs to exclude by LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
   plink --bfile ${output95}/chr_"${chr}" --indep 50 5 5 --out ${Het}/proc/chr_"${chr}" ;
   ##then exclude this list from the working files
   plink --bfile ${output95}/chr_"${chr}" --extract ${Het}/proc/chr_"${chr}".prune.in --make-bed --out ${Het}/proc/chr_"${chr}".temp;
done
#create list to merge
find ${Het}/proc/ -name "*.bim" > ${Het}/proc/allchr.list;
sed -i 's/.bim//g' ${Het}/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${Het}/proc/allchr.list --out ${Het}/proc/full_autosomal_het.temp
#perform heterozigocity
plink --het --bfile ${Het}/proc/full_autosomal_het.temp --out ${Het}/autosomal

#erase temp files
rm ${Het}/proc/*temp*

## create file with samples to exclude (het>4sd) and heterozygosity density plot
Rscript ${EBROOTGAP}/scripts/Het_autosomeQC.R -i ${Het} \
  -o ${repout}

## Create QCed files corrected by heterozygosity 
for chr in {1..22} "XY"
do
  plink --bfile ${output95}/chr_${chr}  \
        --make-bed \
        --remove ${Het}/Excluded.het\
        --out ${Het}/chr_${chr}
done

#### create relatedness(identity by descent) report

for chr in {1..22} "XY" 
do
   ##first, separate a list of SNPs to exclude by MAF and LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
   plink --bfile ${Het}/chr_"${chr}"  --maf 0.01 --indep 50 5 1.3 --out ${Relatedness}/proc/chr_"${chr}" ;
   ##then exclude this list from the working files
   plink --bfile ${Het}/chr_"${chr}" --extract ${Relatedness}/proc/chr_"${chr}".prune.in --make-bed --out ${Relatedness}/proc/chr_"${chr}".temp;
done

#create list to emrge files
find ${Relatedness}/proc/ -name "*.bim" > ${Relatedness}/proc/allchr.list;
sed -i 's/.bim//g' ${Relatedness}/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${Relatedness}/proc/allchr.list --out ${Relatedness}/proc/full_autosomal_rel.temp
#perform identity by descent computation
plink --genome --bfile ${Relatedness}/proc/full_autosomal_rel.temp --min 0.05 --out ${Relatedness}/autosomal_rel

#remove temp files
rm ${Relatedness}/proc/*temp*


##Call the Rscript to plot
Rscript genotypeQC.R -i ${AutosomeQCDir} \
 -o ${repout} \
 -n "test01" \
 --sampleinfo ${samplesheet} \
 -r ${MAFref}










