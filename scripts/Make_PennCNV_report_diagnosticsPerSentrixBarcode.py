from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype
import sys
import os
import glob
import argparse

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file") #"/apps/data/GSAarray/GSAMD-24v1-0_20011747_A5.bpm"
parser.add_argument("gtc_directory", help="Directory containing GTC files per glass number") # "/groups/umcg-gsad//tmp04/projects//GSA-24+v1.0-MD_000617/run01/rawdata/array/"
parser.add_argument("output_directory", help="Directory where output has to be written") #"/groups/umcg-gsad//tmp04//tmp//GSA-24+v1.0-MD_000617/run01//PennCNV_reports/"
parser.add_argument("Sample_ID", help="Sample_ID, array_ID:SentrixBarcode_SentrixPosition") #DNA-119212_OND-545213_EXP-18091_Female:203273700072_R01C01
args = parser.parse_args()

split_Sample_ID = args.Sample_ID.split(':')
array_ID = split_Sample_ID[0]
sentrix_ID = split_Sample_ID[1]

for gtc_file in glob.glob(os.path.join(args.gtc_directory, sentrix_ID+'.gtc')):
	print "gtc_file is" + gtc_file
	save_path = args.output_directory
	print "save_path:" + save_path
	gtc_output = os.path.join(args.output_directory, array_ID + ".txt")
	print "gtc_output is" + gtc_output
	
	output = open(gtc_output, "w")
	print "output is" + str(output)
	if gtc_file.lower().endswith(".gtc"):
		manifest = BeadPoolManifest(args.manifest)
		names = manifest.names
		chrom = manifest.chroms
		map_info = manifest.map_infos 
		map_info2 = str(map_info)
		logratio = GenotypeCalls( gtc_file ).get_logr_ratios()
		BAF = GenotypeCalls( gtc_file ).get_ballele_freqs()
		output.write("[Data]" + "\n")
		output.write("SNP Name" + "\t" + "Sample ID" + "\t" + "Chr" + "\t" + "Position" + "\t" + "Log R Ratio" + "\t" + "B Allele Freq" + "\n")
		for (names, chrom, map_info, logratio, BAF) in zip(names, chrom, map_info, logratio, BAF):
			output.write(names + "\t" + array_ID + "\t" + chrom + "\t" + str(map_info) + "\t" + str(logratio) + "\t" + str(BAF) + "\n")
	output.close()
