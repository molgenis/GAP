#MOLGENIS walltime=00:40:00 mem=40gb ppn=1

#string Project

#string genSampleDir
#string AutosomeQCDir
#string plinkVersion
#string gapVersion
#string output80
#string CR_high
#string outputMH
#string PCA
#string ref1000G
#string repout
#string RPlusVersion
#string MAFref
#string logsDir
#string intermediateDir
#string samplesheet
#X_QCDir
#X_output80
#X_CR_high
#X_outputMH
#X_repout


module load "${plinkVersion}"
module load "${RPlusVersion}"
module load "${gapVersion}"
module list


output80="${AutosomeQCDir}/1_CR80" 
CR_high="${AutosomeQCDir}/2_CR_high"
outputMH="${AutosomeQCDir}/3_MAF_HWE"
Het="${AutosomeQCDir}/4_Het"
Relatedness="${AutosomeQCDir}/5_Relatedness"
PCA="${AutosomeQCDir}/6_PCA"
ref1000G="/apps/data/1000G/phase3"


mkdir -p "${AutosomeQCDir}/"
mkdir -p "${output80}/"
mkdir -p "${output95}/"
mkdir -p "${outputMH}/"
mkdir -p "${repout}/"
mkdir -p "${Het}/"
mkdir -p "${Het}/proc/"
mkdir -p "${Relatedness}/"
mkdir -p "${Relatedness}/proc/"
mkdir -p "${AutosomeQCDir}/6_PCA"
mkdir -p "${AutosomeQCDir}/6_PCA/proc"
mkdir -p "${AutosomeQCDir}/6_PCA/proc2"
 
#################################################################################################################
#########################----------------------Main-----------------#############################################
 
#####################################################################################################
###############--------Call rate and duplicate SNP QC-----###########################################
for chr in {1..22} "XY"
do
### create plink files and call_rate stats for individuals and SNPs
   plink --data ${genSampleDir}/chr_${chr} \
         --make-bed  \
         --missing \
         --out ${AutosomeQCDir}/chr_${chr}
 
##create list of samples to exclude on the criteria callrate<=80
   awk '$6>0.20 {print $1, $2}' ${AutosomeQCDir}/chr_${chr}.imiss > ${output80}/chr_${chr}.extr80_sam
done
 
## create files with duplicate SNPs to exlude 
Rscript position_duplicates.R -i ${AutosomeQCDir} 

##creates list of individuals and dup snps excluded for all the autosomes
cat ${output80}/chr_*.extr80_sam |sort -u > ${output80}/extr80.samples
cat ${AutosomeQCDir}/chr_*.excl.duplicates > ${AutosomeQCDir}/extr.dups
 
for chr in {1..22} "XY"
 do
 # exclude individuals with callrate<=80 (creates excluded individuals_file) creates data set with individual_ callrate>80
   plink --bfile ${AutosomeQCDir}/chr_${chr}  \
         --make-bed \
         --remove ${output80}/extr80.samples \
         --exclude ${AutosomeQCDir}/extr.dups \
         --out ${output80}/chr_${chr}
 
 #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${output80}/chr_${chr} \
          --missing \
          --out ${output80}/chr_${chr}

 awk '$5>0.20 {print $2}' ${output80}/chr_${chr}.lmiss > ${output80}/chr_${chr}.extr80_var
 done
##file of snps to exclude
cat ${output80}/chr_*.extr80_var > ${output80}/extr80.vars
 
 for chr in {1..22} "XY"
 do
   # exclude  SNPs>80
   plink --bfile ${output80}/chr_${chr}  \
         --make-bed \
         --exclude ${output80}/extr80.vars\
         --out ${output80}/chr_${chr}.2

   #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${output80}/chr_${chr}.2 \
          --missing \
          --out ${output80}/chr_${chr}.2

 ##create list of SNPs snd samples to exclude on the criteria callrate<=high
 awk '$6>0.01 {print $1, $2}' ${output80}/chr_${chr}.2.imiss > ${CR_high}/chr_${chr}.extrhigh_sam
 ##information for the heterozygosity analysis
 awk '$6<0.01 {print $1, $2,$6}' ${output80}/chr_${chr}.imiss > ${CR_high}/chr_${chr}.incl_CR_sam
 done
cat ${CR_high}/chr_*.extrhigh_sam|sort -u > ${CR_high}/extrhigh.samples

 for chr in {1..22} "XY"
 do
   ## exclude individuals with callrate<=high (creates excluded individuals_file) creates data set with  individual_ callrate>99
   plink --bfile ${output80}/chr_${chr}  \
         --make-bed \
         --remove ${CR_high}/extrhigh.samples \
         --out ${CR_high}/chr_${chr}
 
  #calculates callrate stats from the previously filtered datafile
   plink  --bfile ${CR_high}/chr_${chr} \
          --missing \
          --out ${CR_high}/chr_${chr}
          
 awk '$5>0.01 {print $2}' ${output80}/chr_${chr}.lmiss > ${CR_high}/chr_${chr}.extrhigh_var
 done
cat ${CR_high}/chr_*.extrhigh_var > ${CR_high}/extrhigh.vars
  
 for chr in {1..22} "XY"
 do
   ## exclude SNPs callrate>99
   plink --bfile ${CR_high}/chr_${chr}  \
         --make-bed \
         --exclude ${CR_high}/extrhigh.vars \
         --out ${CR_high}/chr_${chr}.2

   ## calculate MAF and HWE 
   plink --bfile ${CR_high}/chr_${chr}.2 \
         --freq \
         --hardy \
         --out ${outputMH}/chr_${chr}
 done
 
##creates merged files of included individuals and SNPs to be used of further analysis
 
cat ${AutosomeQCDir}/chr_*.fam|sort -u|awk '{print$2}' > ${AutosomeQCDir}/full.ind
cat ${AutosomeQCDir}/chr_*.bim|awk '{print$2}' > ${AutosomeQCDir}/untouched.snps
sort ${AutosomeQCDir}/extr.dups ${AutosomeQCDir}/untouched.snps|uniq -u >${AutosomeQCDir}/full.snps 

cat ${output80}/chr_*.2.fam|sort -u|awk '{print$2}' > ${output80}/incl80.samples
cat ${output80}/chr_*.2.bim|awk '{print$2}' > ${output80}/incl80.vars
 
cat ${CR_high}/chr_*.2.fam|sort -u|awk '{print$2}' > ${CR_high}/inclhigh.samples
cat ${CR_high}/chr_*.2.bim|awk '{print$2}' > ${CR_high}/inclhigh.vars
 
 
##################################################################################################
################-------------MAF and HWE filtering--------########################################

###Eliminate outlier markers from  H-WE 
for chr in {1..22} "XY"
  do
 
  awk '$9<0.000001 {print $2}' ${outputMH}/chr_${chr}.hwe|sort > ${outputMH}/highhw_${chr}.temp
  awk '$5==0 {print $2}' ${outputMH}/chr_${chr}.frq|sort > ${outputMH}/zeroMAF_${chr}.temp
  cat ${outputMH}/highhw_${chr}.temp ${outputMH}/zeroMAF_${chr}.temp  > ${outputMH}/extr_${chr}hw
  rm *.temp
  done
cat ${outputMH}/extr_*hw> ${outputMH}/excl_HW.snps

for chr in {1..22} "XY"
 do
   ## extract markers with MAF>0.01 and WHE< 1x 10 exp(-6)
   plink --bfile ${CR_high}/chr_${chr}.2 \
         --make-bed \
         --exclude ${outputMH}/excl_HW.snps\
         --out ${outputMH}/chr_${chr}
 done
cat ${outputMH}/chr_*.bim|awk '{print$2}' > ${outputMH}/incl_HW.snps

##################################################################################################
################-------------Heterozygosity for samples--------###################################

#create list to merge files
find ${outputMH}/ -name "*.bim" > ${Het}/proc/allchr.list;
sed -i 's/.bim//g' ${Het}/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${Het}/proc/allchr.list \
       --out ${Het}/proc/full_autosomal_het.temp

##first, separate a list of SNPs to exclude by LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
plink --bfile ${Het}/proc/full_autosomal_het.temp \
      --indep 50 5 2.5 --out ${Het}/proc/full_extract.temp ;
#then exclude this list from the working files
plink --bfile ${Het}/proc/full_autosomal_het.temp \
      --extract ${Het}/proc/full_extract.temp.prune.in \
      --make-bed --out ${Het}/proc/pruned_autosomal_het.temp;
 
#perform heterozigocity
plink --het \
       --bfile ${Het}/proc/pruned_autosomal_het.temp \
       --homozyg \
       --out ${Het}/autosomal
#erase temp files 
rm ${Het}/proc/*temp*
#retrieve Call rate information from samples
cat ${CR_high}/chr_*.incl_CR_sam> ${Het}/CR.samples

## create file with samples to exclude (het>4sd) and heterozygosity density plot
 Rscript Het_autosomeQC.R -i ${Het} \
   -o ${repout}
 
## Create QCed files corrected by heterozygosity 
for chr in {1..22} "XY"
 do
   plink --bfile ${outputMH}/chr_"${chr}" \
         --make-bed \
         --remove ${Het}/Excluded.het \
         --out ${Het}/chr_${chr}
 done
 
#######################################################################################################
############################--------X chromosome QC and sex check--------##############################
###path variables




###Make working directories
mkdir -p "${X_QCDir}/"
mkdir -p "${X_output80}/"
mkdir -p "${X_CR_high}/"
mkdir -p "${X_outputMH}/"
mkdir -p "${X_repout}/"

###############--------Call rate and duplicate SNP QC for X-----#########

cat ${Het}/Excluded.het  ${CR_high}/extrhigh.samples ${output80}/extr80.samples > \
${X_QCDir}/excludebeforeX.samples

### create plink files and call_rate stats for individuals and SNPs
plink --data ${genSampleDir}/chr_X \
     --make-bed  \
     --remove ${X_QCDir}/excludebeforeX.samples \
     --missing \
     --out ${X_QCDir}/chr_X

## generate list of duplicated SNPs (selecting the one with best call rate).
Rscript position_duplicates.R -i ${X_QCDir} 

##creates list of individuals and dup snps excluded for the X chromosome
cat ${X_QCDir}/chr_X.excl.duplicates > ${X_QCDir}/extr.dups

plink --bfile ${X_QCDir}/chr_X  \
     --make-bed \
     --exclude ${X_QCDir}/extr.dups \
     --out ${X_output80}/chr_X

awk '$5>0.20 {print $2}' ${X_QCDir}/chr_X.lmiss > ${X_output80}/chr_X.extr80_var
cat ${X_output80}/chr_*.extr80_var > ${X_output80}/extr80.vars

# exclude SNPs with a call rate < 80%
plink --bfile ${X_output80}/chr_X  \
     --make-bed \
     --exclude ${X_output80}/extr80.vars\
     --out ${X_output80}/chr_X.2

#calculates callrate stats from the previously filtered datafile (removed SNPS with call rate < 80%)
plink  --bfile ${X_output80}/chr_X.2 \
      --missing \
      --out ${X_output80}/chr_X

#### NO SAMPLES are to be removed based on X chromosome call rate. 

awk '$5>0.01 {print $2}' ${X_output80}/chr_X.lmiss > ${X_CR_high}/extrhigh.vars

## exclude SNPs callrate>99
plink --bfile ${X_output80}/chr_X.2 \
     --make-bed \
     --exclude ${X_CR_high}/extrhigh.vars \
     --out ${X_CR_high}/chr_X

## calculate MAF and HWE 
plink --bfile ${X_CR_high}/chr_X \
     --freq \
     --hardy \
     --out ${X_outputMH}/chr_X

#### Files with variants throughout the CR filtering steps. 
cat ${X_QCDir}/chr_*.bim|awk '{print$2}' > ${X_QCDir}/untouched.snps
sort ${X_QCDir}/extr.dups ${X_QCDir}/untouched.snps|uniq -u >${X_QCDir}/full.snps 

cat ${X_output80}/chr_*.2.bim|awk '{print$2}' > ${X_output80}/incl80.vars
cat ${X_CR_high}/chr_*.bim|awk '{print$2}' > ${X_CR_high}/inclhigh.vars


######################### Impute sex and sexcheck ###################################
plink --bfile  ${X_CR_high}/chr_X  --impute-sex --make-bed --out ${X_QCDir}/imputed

###call the Rscript for plotting sex concordance
Rscript sexCheck_genotypeQC.R -i ${X_QCDir}/imputed.sexcheck \
                              -p $samplesheet \
                              -o $X_repout
                              
grep -E 'Non concordant|Failed'  ${X_QCDir}/plots/all.samples.concordance.txt| awk \
'{print $4, $4}' > ${X_QCDir}/sex.exclude

####################Filter out SNPs based on HW and MAF for chr X#############################

###Eliminate outlier markers from  H-WE 
  awk '$9<0.000001 {print $2}' ${X_outputMH}/chr_X.hwe > ${X_outputMH}/highhw_X.temp
  awk '$5==0 {print $2}' ${X_outputMH}/chr_X.frq > ${X_outputMH}/zeroMAF_X.temp
  #merge both HW and MAF SNPs
  cat ${X_outputMH}/highhw_X.temp ${X_outputMH}/zeroMAF_X.temp > ${X_outputMH}/extr_Xhw
   ## extract markers with MAF>0.01 and WHE< 1x 10 exp(-6)
   plink --bfile ${X_CR_high}/chr_X \
         --make-bed \
         --remove ${X_QCDir}/sex.exclude \
         --exclude ${X_outputMH}/extr_Xhw \
         --out ${X_outputMH}/chr_X

cat ${X_outputMH}/chr_X.bim|awk '{print$2}' > ${X_outputMH}/incl_HW.snps

 
##################################################################################################
###############---------Relatedness and identity by descent (IBD)------###########################

#create list to merge files
find ${Het}/ -name "*.bim" > ${Relatedness}/proc/allchr.list;
sed -i 's/.bim//g' ${Relatedness}/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${Relatedness}/proc/allchr.list \
       --out ${Relatedness}/proc/full_autosomal_rel.temp
#separate the HLA SNPs
awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' \
${Relatedness}/proc/full_autosomal_rel.temp.bim > ${Relatedness}/proc/HLAexclude.txt

## apply filters to reduce SNP number
 plink --bfile ${Relatedness}/proc/full_autosomal_rel.temp \
       --exclude ${Relatedness}/proc/HLAexclude.txt \
       --make-bed \
       --out ${Relatedness}/proc/full_data
##first, separate a list of SNPs to exclude by MAF and LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
 plink --bfile ${Relatedness}/proc/full_data \
       --maf 0.1 --indep 50 5 1.1 \
       --out ${Relatedness}/proc/full_extract;
##then exclude this list from the working files
 plink --bfile ${Relatedness}/proc/full_data \
       --extract ${Relatedness}/proc/full_extract.prune.in \
       --make-bed \
       --out ${Relatedness}/proc/pruned_full.temp;

#perform identity by descent computation (excluding, people with  pi-hat less than 0.05, this means people completely unrelated)
plink --genome \
      --bfile ${Relatedness}/proc/pruned_full.temp \
      --min 0.05 \
      --out ${Relatedness}/autosomal_rel
#remove temp files
rm ${Relatedness}/proc/*temp*
#create duplicate samples or tweeen file
awk '$9>0.99 {print $2, $3, $7, $8, $9}' ${Relatedness}/autosomal_rel.genome \
>${Relatedness}/equal.samples

##################################################################################################
###########################-----------------PCA analysis---------------###########################

### filter reference Databases by the common SNPs
 plink --bfile ${ref1000G}/1000G_all \
       --extract  ${outputMH}/incl_HW.snps \
       --make-bed \
       --out ${PCA}/proc/thg

 plink --bfile ${gonlref}/gonl_SNV_INDELs_ab \
       --extract  ${outputMH}/incl_HW.snps \
       --make-bed \
       --out ${PCA}/proc/goref
       
#list of common snps
awk '{print $2}' ${PCA}/proc/thg.bim > ${PCA}/proc2/common.snps
#create list to merge
find ${Het}/ -name "*.bim" > ${PCA}/proc2/allfiles.list;
sed -i 's/.bim//g' ${PCA}/proc2/allfiles.list;
###merge all chromosomes
plink --merge-list ${PCA}/proc2/allfiles.list \
      --out ${PCA}/proc2/allchr_join

##make list of individuals to merge 
echo "${PCA}/proc2/allchr_join" >> ${PCA}/proc2/allfiles.list
echo "${PCA}/proc/thg" >> ${PCA}/proc2/allfiles.list
echo "${PCA}/proc/goref" >> ${PCA}/proc2/allfiles.list

###merge all data
 plink --merge-list ${PCA}/proc2/allfiles.list \
       --out ${PCA}/proc2/full_data
###remove SNPs with more tha 3 alleles if there are
 if [ -e ${PCA}/proc2/full_data.missnp ];
  then
      plink --bfile ${PCA}/proc/thg \
            --exclude ${PCA}/proc2/full_data.missnp \
            --make-bed \
            --out ${PCA}/proc2/allg.temp
   
      plink --bfile ${PCA}/proc/goref \
            --exclude ${PCA}/proc2/full_data.missnp \
            --make-bed \
            --out ${PCA}/proc2/allgonl.temp

      find ${Het}/ -name "*.bim" > ${PCA}/proc2/allfiles.list;
      sed -i 's/.bim//g' ${PCA}/proc2/allfiles.list;
      echo "${PCA}/proc2/allchr_join" >> ${PCA}/proc2/allfiles.list
      echo "${PCA}/proc2/allg.temp" >> ${PCA}/proc2/allfiles.list
      echo "${PCA}/proc2/allgonl.temp" >> ${PCA}/proc2/allfiles.list

      plink --merge-list ${PCA}/proc2/allfiles.list \
            --out ${PCA}/proc2/full_data
 fi

#separate the HLA SNPs
awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' \
${PCA}/proc2/full_data.bim > ${PCA}/proc2/HLAexclude.txt
## apply filters to reduce SNP number
 plink --bfile ${PCA}/proc2/full_data \
      --extract ${PCA}/proc2/common.snps \
      --maf 0.1 \
      --geno 0.01 \
      --exclude ${PCA}/proc2/HLAexclude.txt \
      --make-bed \
      --out ${PCA}/proc2/full_data

###prune all data
##first, separate a list of SNPs to exclude by LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
plink --bfile ${PCA}/proc2/full_data \
      --indep 150 20 1.4 \
      --out ${PCA}/proc2/splitsnp ;
##then exclude this list and the HLA from the working files
plink --bfile ${PCA}/proc2/full_data \
      --extract ${PCA}/proc2/splitsnp.prune.in \
      --make-bed \
      --out ${PCA}/PCA_data1
##PCA analysis
plink --bfile ${PCA}/PCA_data1 --pca --out ${PCA}/PCA_1000G

###first PCA plot
Rscript PCA1.R -i ${PCA} \
 -f "PCA_1000G.eigenvec" \
 -o ${repout} \
 -r "1st."

### 2.nd PCA 
##extract selected samples by PCA
plink --bfile ${PCA}/PCA_data1 \
      --keep ${PCA}/1st.eur.samples \
      --make-bed \
      --out ${PCA}/PCA_data2

#2nd. PCA
plink --bfile ${PCA}/PCA_data2 --pca --out ${PCA}/second_PCA

#Rscript PCA1_v1.R -i ${PCA} \
Rscript PCA1.R -i ${PCA} \
 -f "second_PCA.eigenvec" \
 -o ${repout} \
 -r "2nd."

cat ${PCA}/1st.excl.samples ${PCA}/2nd.excl.samples > ${PCA}/excl_PCA.samples

for chr in {1..22} "XY"
 do
 ##extract selected samples by secondPCA
 plink --bfile ${Het}/chr_${chr} \
       --keep ${PCA}/2nd.incl.samples \
       --make-bed \
       --freq \
       --out ${PCA}/chr_${chr}
 done

#######################################################
#################plot results###########################
Rscript genotypeQC_v8.R -i ${AutosomeQCDir} \
 -o ${repout} \
 --sampleinfo ${samplesheet} \
 -r ${MAFref}
 
 
 #END

