#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G

#############################################
### convert from oxforfile (gen-sample) to plink files
### date=01-03-2022
### version: 1
### author: EL
#############################################
### New
## 01-03-2022
## removed chromosome loop and adequate to a launcher script to papralellize by chromosome
#############################################
### inform node
echo hostname
### Load packages
module load PLINK/1.9-beta6-20190617
### inform deifnition of variables 
echo "wkdir =" ${1}  ## the place used to process the data and make output
echo "inputdir=" ${2}  ## directory with the folders fo each batch
echo "out=" ${3} ### name of the output folder
echo "logs=" ${1}/${3}/logs ## a logfiles directory
echo "chr=" ${4} ### cromosome processed

################################################################################
### Main
################################################################################
mkdir -p ${1}/${3}/ #make results directory
mkdir -p ${1}/${3}/logs/ #make logs directory
mkdir -p ${1}/${3}/meta/ # make directory for metafiles (used independently of chromosomes)

      plink --data ${2}/chr_${4} \
            --make-bed  \
            --out ${1}/${3}/chr_${4}
      echo "${1}/${3}/chr_${4}" >> ${1}/${3}/meta/allb_${4}.list









