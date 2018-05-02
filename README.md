# GAP 
   Short for Genotyping Array Pipeline.
   Consist of the following workflow:

```
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                    Iscan writes IDAT files to GATTACA {01,02}machines                ⎜
   ⎜                                                                                      ⎜
   ⎝______________________________________________________________________________________⎠
                                         v
                                         v
                                         v
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                    DARWIN conversion of IDAT files to GTC files                      ⎜
   ⎜                    takes place on GATTACA {01,02}machines                            ⎜
   ⎝______________________________________________________________________________________⎠
                                         v
                                         v  > > > > > > > > > > GAP_Automated CopyRawDataToPRM [stores .IDAT and .GTC files on permanent storage system]
                                         v  > > > > > > > > > > GAP_Automated Start Pipeline [automatically start Pipeline when new run has finished]
                                         v
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                         Run GAP Pipeline                                             ⎜
   ⎜                         takes place on [zinc-finger of leucine-zipper]               ⎜
   ⎝______________________________________________________________________________________⎠
                                         v
                                         v  > > > > > > > > > > GAP_Automated CopyProjectDataToPRM [stores data produced by pipeline on permanent storage system]
                                         v
   ⎛¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯⎞
   ⎜                  DARWIN split pipeline output file per sample,                       ⎜
   ⎜                  calculate  standard deviation and store all data in array database  ⎜
   ⎝______________________________________________________________________________________⎠
```

#### GAP Pipeline 

The GAP pipeline consist of 3 steps:

1 Create_Callrate_file
```
This step creates a file containing Callrate information per sample.
The Callrate gives an idea of how many percent of the SNPs on the array are performing well.
If this number is below 0.97 we know the data is of pour quality.

Fileformat:

Sample ID '/t' Call Rate '/t' Gender
```
2 Make_Final_PennCNV_report
```
This step creates a file containing per SNP information about the log ratio and the B allel frequency of the specific snp.

Fileformat:
Name '/t' Chromosome '/t' Position '/t' Sample1.GType '/t' Sample1.LogRratio '/t' Sample1.B Allele Freq '/t' Sample2.GType '/t' Sample2.LogRratio '/t' Sample2.B Allele Freq

The log R ratio and B allele Frequency  per SNP are used by Nexus (commercial software) to call CNV's
```
3 CopyToResults Dir
```
This step copies the results to the ${projectname}/${resultsDir}
```

GAP_Automated steps which are not implemented yet are:
```
GAP_Automated CopyProjectDataToPRM
```
```
GAP_Automated CopyRawdataToPRM
```

GAP_pipeline steps yet to come:
```
Create array input for concordance check NGS_Data 
Add extra location to put output files from a project so DARWIN can use this as input
 ```
