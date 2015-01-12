#!/usr/bin/python
import argparse
import sys
import os
import re
import csv
import numpy
from pprint import pprint
__author__ = 'Robin van der Weide'

#Command-line thingies
parser = argparse.ArgumentParser(
    description='This is TRIfle: integrated non-coding variant prioritisation. Please give ONE input-file.')
parser.add_argument('-i', '--input', help='TRIfle-getRawScores.py output', required=True)
args = vars(parser.parse_args())


#Main script
info = []
commonfile = open('./commonCandidates','w+')
rcfile = open('./rareCodingCandidates','w+')
rncfile = open('./rareNcCandidates','w+')

ifile = open(args['input'], 'r')
for row in ifile:
    if row.startswith('#'):
        continue
    else:
        fields = re.split(r'\t+', row)
        #coord = re.split(r'\|+', fields)
        id = fields[0].split('|' )
        coord = id[0].split(':')
        outputLine = str(coord[0]) + str(":") + str(coord[1]) + str("-") + str(int(coord[1])+1) + "\t" + id[1] + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\t" + fields[4] + "\t" + str(int(float(fields[5]))) + "\t" + str(int(float(fields[6])))+ "\n"
        if float(fields[1]) == float(0.0):
            if float(fields[2]) == float(0.0):
                cr = "common"
            else:
                cr = "rare"
        else:
            cr = "rare"
        if cr == "common":
            if float(fields[3]) > float(10.0):
                commonfile.write(outputLine)
        elif cr == "rare":
            if float(fields[3]) > float(10):
                if float(fields[1]) == float(0.0):
                #dit is een coding rare
                    rcfile.write(outputLine)
                    #print(float(fields[1]))
                if float(fields[2]) == float(0.0):
                #dit is een nc rare
                    rncfile.write(outputLine)
        #print(cr)
ifile.close()

#to make a direct vcf of this output: cat rareCodingCandidates | \
# sed 's/:/  /g' | sed 's/chr//g' | sed 's/-/       /g' | sed 's$/$ $g' | \
# awk '{print $1"\t"$2"\t.\t"$4"\t"$5"\t30\tRANK=XXX;CASE="$11";CONTROL="$10}' > rareCodingCandidates.vcf