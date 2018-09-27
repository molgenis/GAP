outputFolder='/groups/umcg-gap/tmp04/projects/test/new/s2_genfiles/'
#mkdir -p "/groups/umcg-gap/tmp04/rawdata/array/GTC/test.small/genFilesOut"
for chr in {1..22} "X" "XY" "Y" "MT"
do

echo "#!/bin/bash
#SBATCH --job-name=convert_to_Gen_chr${chr}
#SBATCH --output=convert_to_Gen_chr${chr}.out
#SBATCH --error=convert_to_Gen_chr${chr}.err
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task 8
#SBATCH --mem 20gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=30L

inputFolder='/groups/umcg-gap/tmp04/projects/test/new/s1_optical/'
outputFolder='/groups/umcg-gap/tmp04/projects/test/new/s2_genfiles/'
chr=${chr}

if [[ -z \"\$inputFolder\" ]]
then
    echo \"Error: set input\"
    exit
fi


if [[ -z \"\$outputFolder\" ]]
then
    echo \"Error: set output\"
    exit
fi

    echo \"\${chr}\"

    input=\"\${inputFolder}/chr_\${chr}.probs\"
    output=\"\${outputFolder}/chr_\${chr}\"

    sampleFile=\"\${output}.sample\"
    
    if [ -e \$input ]
    then

        awk -v chr=\${chr} -v sampleFile=\"\${sampleFile}\" '

                NR == 1 {
                        print \"ID_1 ID_2 missing\" > sampleFile
                        print \"0 0 0\" > sampleFile
                        sampleCount = 0
                        for(i=5;i<=NF;i+=1){
                                print \$i,\$i,0 > sampleFile
                                sampleCount += 1
                        }
                }

                NR > 1 {

                        ORS=\" \";
                        print chr, \$1, \$2, substr(\$3,1,1), substr(\$3,2,1);
                        for(i=0;i<sampleCount;i+=1){
                                x = 4 + i * 4
                                print \$x,\$(x+2),\$(x+1)
                        }
                        ORS=\"\";
                        print \"\\n\"
                }
        ' < \${input} > \${output}.gen
fi

" >${outputFolder}/convertToGen_chr${chr}.sh
done

echo "Done"
