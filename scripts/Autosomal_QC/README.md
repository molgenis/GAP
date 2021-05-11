# Quality Control for GSA array on the UMCG HPC
Created by: Raul Aguirre-Gamboa, Esteban Lopera-Maya. \
Contact information: e.a.lopera.maya@umcg.nl. \

# Pipeline structure
Important note: this is a self-contained and specific pipeline to work on the UMCG HPC. This means, that at the time it is not fully automated to work in completeness in every kind of data and every kind of server, and many steps might be specifically designed for the UGLI data alone. You are welcome, however, to take any individual step and/or code snippet to addapt to your own conditions. \
\
This pipeline is designed to work in three main steps as indicated by the numeration in front of the name of the main scripts. Numbered scripts should be called by the user independently one ofter the other. All the scripts with the prefix "sub_" in the front of the name are automatically by the numbered scripts. \
A feedback loop should be done manually by the user, after processing also manually the output of step 2 (script 2.) to remove familial errors (as these cannot be removed automatically), the data resulting from this should go trhough step 2 a second time. Special adittional "sub_" scripts contain also the suffix "second_it". These should replace their counterparts in the second iteration of step 2. \

#  Location and order of the quality control steps




