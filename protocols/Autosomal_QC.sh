#MOLGENIS walltime=00:40:00 mem=40gb ppn=1

#string Project

#string genSampleDir
#string AutosomeQCDir
#string plinkVersion
#string gapVersion
#string output80
#string output95
#string outputMH
#string PCA
#string ref1000G
#string repout
#string RPlusVersion
#string MAFref
#string logsDir
#string intermediateDir
#string samplesheet

 

module load "${plinkVersion}"
module load "${RPlusVersion}"
module load "${gapVersion}"
module list


output80="${AutosomeQCDir}/1_CR80" 
output95="${AutosomeQCDir}/2_CR95"
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
mkdir -p "${AutosomeQCDir}/6_PCA/proc2"

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
 
 
 #### heterozygosity analysis by sample
 
 for chr in {1..22} "XY" 
 do
    ##first, separate a list of SNPs to exclude by LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
    plink --bfile ${CR_high}/chr_"${chr}".2 \
          --indep 50 5 5 --out ${Het}/proc/chr_"${chr}" ;
    ##then exclude this list from the working files
    plink --bfile ${CR_high}/chr_"${chr}".2 \
          --extract ${Het}/proc/chr_"${chr}".prune.in \
          --make-bed --out ${Het}/proc/chr_"${chr}".temp;
 done
 #create list to merge
 find ${Het}/proc/ -name "*.bim" > ${Het}/proc/allchr.list;
 sed -i 's/.bim//g' ${Het}/proc/allchr.list;
 #merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
 plink --merge-list ${Het}/proc/allchr.list \
       --out ${Het}/proc/full_autosomal_het.temp
 #perform heterozigocity
 plink --het \
       --bfile ${Het}/proc/full_autosomal_het.temp \
       --homozyg \
       --out ${Het}/autosomal
 
 #erase temp files
 rm ${Het}/proc/*temp*
  #retrieve CR information from amples
 cat ${CR_high}/chr_*.incl_CR_sam> ${Het}/CR.samples

 ## create file with samples to exclude (het>4sd) and heterozygosity density plot
 Rscript Het_autosomeQC.R -i ${Het} \
   -o ${repout}
 
 ## Create QCed files corrected by heterozygosity 
 for chr in {1..22} "XY"
 do
   plink --bfile ${CR_high}/chr_${chr}.2  \
         --make-bed \
         --remove ${Het}/Excluded.het \
         --out ${Het}/chr_${chr}
 done
 
 #### create relatedness(identity by descent) report
 
 for chr in {1..22} "XY" 
 do
    ##first, separate a list of SNPs to exclude by MAF and LD (--indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
    plink --bfile ${Het}/chr_"${chr}" \
          --maf 0.01 --indep 50 5 1.3 \
          --out ${Relatedness}/proc/chr_"${chr}" ;
    ##then exclude this list from the working files
    plink --bfile ${Het}/chr_"${chr}" \
          --extract ${Relatedness}/proc/chr_"${chr}".prune.in \
          --make-bed \
          --out ${Relatedness}/proc/chr_"${chr}".temp;
 done
 
 #create list to emrge files
 find ${Relatedness}/proc/ -name "*.bim" > ${Relatedness}/proc/allchr.list;
 sed -i 's/.bim//g' ${Relatedness}/proc/allchr.list;
#merge all chromosomes into one single genotype ...set of files (.fam, .bim, .bed)
plink --merge-list ${Relatedness}/proc/allchr.list \
      --out ${Relatedness}/proc/full_autosomal_rel.temp
#perform identity by descent computation (excluding, people with  pi-hat less than 0.05, this means people completely unrelated)
plink --genome \
      --bfile ${Relatedness}/proc/full_autosomal_rel.temp \
      --min 0.05 \
      --out ${Relatedness}/autosomal_rel
 #remove temp files
rm ${Relatedness}/proc/*temp*
 


####PCA analysis



### filter reference DB by the common SNPs
plink --bfile ${ref1000G}/1000G_all \
      --extract  ${CR_high}/inclhigh.vars \
      --make-bed \
      --out ${thg}/thg

awk '{print $2}' ${thg}/thg.bim > ${PCA}/proc2/common.snps


#create list to merge
 find ${Het}/ -name "*.bim" > ${PCA}/proc2/allfiles.list;
 sed -i 's/.bim//g' ${PCA}/proc2/allfiles.list;
 echo "${thg}/thg" >> ${PCA}/proc2/allfiles.list

###merge all data
plink --merge-list ${PCA}/proc2/allfiles.list \
      --out ${PCA}/proc2/full_data

if [ -e ${PCA}/proc2/full_data.missnp ];
  then
      plink --bfile ${thg}/thg \
            --exclude ${PCA}/proc2/full_data.missnp \
            --make-bed \
            --out ${PCA}/proc2/allg.temp
      
      find ${Het}/ -name "*.bim" > ${PCA}/proc2/allfiles.list;
      sed -i 's/.bim//g' ${PCA}/proc2/allfiles.list;
      echo "${PCA}/proc2/allg.temp" >> ${PCA}/proc2/allfiles.list


      plink --merge-list ${PCA}/proc2/allfiles.list \
            --out ${PCA}/proc2/full_data
fi

#separate the HLA SNPs
awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' ${PCA}/proc2/full_data.bim > ${PCA}/proc2/HLAexclude.txt

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
      --indep 300 10 1.1 \
      --out ${PCA}/proc2/splitsnp ;
##then exclude this list and the HLA from the working files
plink --bfile ${PCA}/proc2/full_data \
      --extract ${PCA}/proc2/splitsnp.prune.in \
      --make-bed \
      --out ${PCA}/G_filt_PCA2

##PCA analysis
plink --bfile ${PCA}/G_filt_PCA2 --pca --out ${PCA}/PCA_1000G



##Call the Rscript to plot
Rscript genotypeQC.R -i ${AutosomeQCDir} \
 -o ${repout} \
 -n "test01" \
 --sampleinfo ${samplesheet} \
 -r ${MAFref}










