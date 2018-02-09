from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype

import sys

import os

import glob

import argparse

parser = argparse.ArgumentParser("Generate a Callrate report from a directory of GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_directory", help="Directory containing GTC files")
parser.add_argument("output_file", help="Directory where output has to be written")

args = parser.parse_args()

output_file2 = args.output_file
output_file3 = open(output_file2, "w")
print "output is" + str(output_file3)

for gtc_file in glob.glob(os.path.join(args.gtc_directory, '*.gtc')):

	print "processing" + gtc_file

        manifest = BeadPoolManifest(args.manifest)

        sample_id = os.path.basename(gtc_file)[:-4]

        call_rate = str(GenotypeCalls(gtc_file).get_call_rate())

        gender = GenotypeCalls(gtc_file).get_gender()

	output_file3.write(sample_id + "\t" + call_rate + "\t" + gender + "\n")



output_file3.close()
