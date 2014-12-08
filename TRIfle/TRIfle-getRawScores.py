#!/usr/bin/python
import argparse
import sys
import re
import csv
import numpy
from pprint import pprint
__author__ = 'Robin van der Weide'

#Command-line thingies
parser = argparse.ArgumentParser(
    description='This is TRIfle: integrated non-coding variant prioritisation. Please give ONE input-file.')
parser.add_argument('-f', '--Funseq2', help='input: list of Output.vcf', required=False)
parser.add_argument('-c', '--CADD', help='input: list of scores.tsv', required=False)
parser.add_argument('-d', '--DANN', help='input: list of scores.tsv', required=False)
parser.add_argument('-s', '--SuRFR', help='input: list of output.tsv', required=False)
#parser.add_argument('-i', '--input', help='input: list of variant calls', required=True)
#parser.add_argument('-cc','--Case-Control', help='Case(1) or Control(0)', required=True, default='1')
args = vars(parser.parse_args())
if not any(args.values()):
    parser.error("Error: EITHER -f, -s or -c must be given. Use -h for help.")


#FUNCTIONS
def funseq2scorecount(Ffile, scorebookFN, scorebookFC, countbookF):
    score = 0
    coord = ""
    count = 0
    scoreN = 0.0
    scoreS = 0.0
    for row in Ffile:
        if row.startswith('#'):
            continue
        else:
            fields = re.split(r'\t+', row)
            coord = str(fields[0]) + str(":") + str(fields[1] + str("|") + str(fields[4]))
            index = row.find("CDSS", 40)
            if not index == -1:
                inforow = row[index:]
                CDSSfield = re.split(r'[;]', inforow)[0]
                score = float(re.split(r'[=]+', CDSSfield)[-1].rstrip())
                if coord in scorebookFC:
                    continue
                else:
                    scorebookFC[coord] = float(score)
            else:
                index = row.find("NCDS", 40)
                if not index == -1:
                    inforow = row[index:]
                    CDSSfield = re.split(r'[;]', inforow)[0]
                    score = float(re.split(r'[=]+', CDSSfield)[-1].rstrip())
                    if coord in scorebookFN:
                        continue
                    else:
                        scorebookFN[coord] = float(score)
            if coord in countbookF:
                countbookF[coord] += 1
            else:
                countbookF[coord] = float(1)
    return(scorebookFN, scorebookFC, countbookF)
def cadd2scorecount(Ffile, scorebookC, countbookC):
    score = 0.0
    coord = ""
    for row in Ffile:
        if row.startswith('#'):
            continue
        else:
            fields = re.split(r'\t+', row)
            coord = str(str("chr") + fields[0] + str(":") + str(fields[1]) + str("|") + str(fields[3]))
            score = float(fields[5])
            if coord in scorebookC:
                if float(scorebookC[coord]) < float(score):
                    scorebookC[coord] = float(score)
            else:
                    scorebookC[coord] = score
            if coord in countbookC:
                countbookC[coord] += 1
            else:
                countbookC[coord] = float(1)
    return(scorebookC, countbookC)


#General declarations
countbookF = {}
scorebookFN = {}
scorebookFC = {}
scorebookC = {}
countbookC = {}
scorebookD = {}
countbookD = {}
variantCoordList = []

#Main script


if args['Funseq2'] is not None:
    Ffile = open(args['Funseq2'], 'r')
    for row in Ffile:
        frow = open(row.rstrip(), 'r')
        scorebookFN, scorebookFC, countbookF = funseq2scorecount(frow, scorebookFN, scorebookFC, countbookF)
        frow.close()
    Ffile.close()
if args['CADD'] is not None:
    Cfile = open(args['CADD'], 'r')
    for row in Cfile:
        crow = open(row.rstrip(), 'r')
        scorebookC, countbookC = cadd2scorecount(crow, scorebookC, countbookC)
        crow.close()
    Cfile.close()
if args['SuRFR'] is not None:
    Sfile = open(args['SuRFR'], 'r')
if args['DANN'] is not None:
    Cfile = open(args['DANN'], 'r')
    for row in Cfile:
        crow = open(row.rstrip(), 'r')
        scorebookD, countbookD = cadd2scorecount(crow, scorebookD, countbookD)
        crow.close()
    Cfile.close()



# Adding zeros to missing scores
for key,value in scorebookD.items():
    if key in scorebookFC:
        continue
    else:
        scorebookFC[key] = float(0)
    if key in scorebookFN:
        continue
    else:
        scorebookFN[key] = float(0)
    if key in scorebookC:
        continue
    else:
        scorebookC[key] = float(0)
for key,value in scorebookC.items():
    if key in scorebookFC:
        continue
    else:
        scorebookFC[key] = float(0)
    if key in scorebookFN:
        continue
    else:
        scorebookFN[key] = float(0)
    if key in scorebookD:
        continue
    else:
        scorebookD[key] = float(0)
for key,value in scorebookFC.items():
    if key in scorebookC:
        continue
    else:
        scorebookC[key] = float(0)
    if key in scorebookFN:
        continue
    else:
        scorebookFN[key] = float(0)
    if key in scorebookD:
        continue
    else:
        scorebookD[key] = float(0)
for key,value in scorebookFN.items():
    if key in scorebookC:
        continue
    else:
        scorebookC[key] = float(0)
    if key in scorebookFC:
        continue
    else:
        scorebookFC[key] = float(0)
    if key in scorebookD:
        continue
    else:
        scorebookD[key] = float(0)
#same with counts
for key,value in countbookD.items():
    if key in countbookF:
        continue
    else:
        countbookF[key] = float(value)
    if key in countbookC:
        continue
    else:
        countbookC[key] = float(value)
for key,value in countbookC.items():
    if key in countbookF:
        continue
    else:
        countbookF[key] = float(value)
    if key in countbookD:
        continue
    else:
        countbookD[key] = float(value)
for key,value in countbookF.items():
    if key in countbookC:
        continue
    else:
        countbookC[key] = float(value)
    if key in countbookD:
        continue
    else:
        countbookD[key] = float(value)

print(len(countbookC))
print(len(countbookF))
print(len(countbookD))
print(len(scorebookFN))
print(len(scorebookFC))
print(len(scorebookC))
print(len(scorebookD))

#print tab-delim: coord:FC:FN:C:Count
#header = ["#Coord|mut","Funseq2(Nc)","Funseq2(C)","CADD(Phred)","DANN","Frequency"]
#print('\t'.join(map(str,header)))
#for coord,score in scorebookD.items():
#    row = [coord, score, scorebookFC[coord],scorebookC[coord],scorebookD[coord],countbookC[coord]]
#    print('\t'.join(map(str,row)))