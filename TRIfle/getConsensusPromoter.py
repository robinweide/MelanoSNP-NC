#!/usr/bin/python
import sys
import argparse
import re

__author__ = 'Robin van der Weide'

# Command-line thingies
parser = argparse.ArgumentParser(
    description='Makes fasta per region of INPUTBED')
parser.add_argument('-i', '--inputBED', help='BED-file with regions (e.g. promoters)', required=True)
parser.add_argument('-a', '--cases', help='List of VCFs', required=True)
parser.add_argument('-b', '--controls', help='List of VCFs', required=True)
parser.add_argument('-r', '--reference', help='Path to reference.fa', required=True)
args = vars(parser.parse_args())

regions = open(args['inputBED'], 'r')
refFile = (args['reference'])
caseFiles = open(args['cases'], 'r')
controlFiles = open(args['controls'], 'r')

for region in regions:
	Chr,start,stop,strand,gene,transcript = re.split(r'\t+', region.rstrip())
	coord = str(Chr + ":" + start + "-" + stop)
	for vcf in caseFiles:
		os.system("samtools faidx " + refFile + " " + coord + " | vcf-consensus " + vcf + " >> " + transcript + "_case.fasta")
	for vcf in controlFiles:
			os.system("samtools faidx " + refFile + " " + coord + " | vcf-consensus " + vcf + " >> " + transcript + "_control.fasta")