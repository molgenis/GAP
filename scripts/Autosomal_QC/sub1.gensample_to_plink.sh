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
echo "chr=" ${3} ### cromosome processed

################################################################################
### Main
################################################################################
       if [[ ${8} == "1" ]]; then
      
      plink --data ${2}/chr_${3} \
            --allow-no-sex \
            --make-bed  \
            --out ${1}/chr_${3}
      else
            plink --data ${2}/chr_${3} \
            --allow-no-sex \
            --set-hh-missing \
            --make-bed  \
            --out ${1}/chr_${3}
      fi







