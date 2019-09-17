from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype
from datetime import date
import sys
import os
import glob
import argparse

parser = argparse.ArgumentParser("Generate a final report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files per glass number")
parser.add_argument("output_directory", help="Directory where output has to be written")
parser.add_argument("Sample_ID", help="Sample_ID, array_ID:SentrixBarcode_SentrixPosition")
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
		delim="\t"
	    	output.write("[Header]\n")
		output.write(delim.join(["GSGT Version", "1.9.4"]) + "\n")
		output.write(delim.join(["Processing Date", date.today().strftime("%d/%m/%Y")]) + "\n")
		output.write(delim.join(["Content", os.path.splitext(manifest.manifest_name)[0]]) + "\n")
    		output.write(delim.join(["Num SNPs", str(len(names))]) + "\n")
		output.write(delim.join(["Total SNPs", str(len(names))]) + "\n")
		output.write(delim.join(["Num Samples", str(1)]) + "\n")
		output.write(delim.join(["Total Samples", str(1)]) + "\n")
		output.write("[Data]" + "\n")
		output.write("SNP Name" + "\t" + "Sample ID" + "\t" + "Chr" + "\t" + "Position" + "\t" + "Log R Ratio" + "\t" + "B Allele Freq" + "\n")
		for (names, chrom, map_info, logratio, BAF) in zip(names, chrom, map_info, logratio, BAF):
			output.write(names + "\t" + array_ID + "\t" + chrom + "\t" + str(map_info) + "\t" + str(logratio) + "\t" + str(BAF) + "\n")
	output.close()
