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
parser.add_argument('-cd', '--controlDANN', help='input: list of control scores.tsv', required=False)
#parser.add_argument('-i', '--input', help='input: list of variant calls', required=True)
#parser.add_argument('-cc','--Case-Control', help='Case(1) or Control(0)', required=True, default='1')
args = vars(parser.parse_args())
if not any(args.values()):
    parser.error("Error: EITHER -f, -s or -c must be given. Use -h for help.")


#FUNCTIONS
def funseq2scorecount(Ffile, scorebookFN, scorebookFC, countbookF, variantCoordList):
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
            variantCoordList.append(coord)
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
    return(scorebookFN, scorebookFC, countbookF, variantCoordList)
def cadd2scorecount(Ffile, scorebookC, countbookC, variantCoordList):
    score = 0.0
    coord = ""
    for row in Ffile:
        if row.startswith('#'):
            continue
        else:
            fields = re.split(r'\t+', row)
            coord = str(str("chr") + fields[0] + str(":") + str(fields[1]) + str("|") + str(fields[3]))
            variantCoordList.append(coord)
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
    return(scorebookC, countbookC, variantCoordList)


#General declarations
countbookF = {}
scorebookFN = {}
scorebookFC = {}
scorebookC = {}
countbookC = {}
scorebookD = {}
countbookD = {}
cscorebookD = {}
ccountbookD = {}
variantCoordList = []
cvariantCoordList = []
controlcounter = float(0)
casecounter = float(0)

#Main script


if args['Funseq2'] is not None:
    Ffile = open(args['Funseq2'], 'r')
    for row in Ffile:
        frow = open(row.rstrip(), 'r')
        scorebookFN, scorebookFC, countbookF, variantCoordList = funseq2scorecount(frow, scorebookFN, scorebookFC, countbookF, variantCoordList)
        frow.close()
    Ffile.close()
if args['CADD'] is not None:
    Cfile = open(args['CADD'], 'r')
    for row in Cfile:
        crow = open(row.rstrip(), 'r')
        scorebookC, countbookC, variantCoordList = cadd2scorecount(crow, scorebookC, countbookC,variantCoordList)
        crow.close()
    Cfile.close()
if args['SuRFR'] is not None:
    Sfile = open(args['SuRFR'], 'r')
if args['DANN'] is not None:
    Cfile = open(args['DANN'], 'r')
    for row in Cfile:
        crow = open(row.rstrip(), 'r')
        scorebookD, countbookD, variantCoordList = cadd2scorecount(crow, scorebookD, countbookD, variantCoordList)
        casecounter +=1
        crow.close()
    Cfile.close()
if args['controlDANN'] is not None:
    Cfile = open(args['controlDANN'], 'r')
    for row in Cfile:
        crow = open(row.rstrip(), 'r')
        cscorebookD, ccountbookD, cvariantCoordList = cadd2scorecount(crow, cscorebookD, ccountbookD, cvariantCoordList)
        controlcounter +=1
        crow.close()
    Cfile.close()


coords = list(set(variantCoordList))
coords.sort()

#print tab-delim: coord:FC:FN:C:Count
header = ["#Coord|mut","Funseq2(Nc)","Funseq2(C)","CADD(Phred)","DANN","FrequencyCase","FrequencyControl","PercCase","PercControl"]
print('\t'.join(map(str,header)))

for coord in coords:
    row = [coord, scorebookFN.get(coord, float(0)),scorebookFC.get(coord, float(0)), scorebookC.get(coord, float(0)),scorebookD.get(coord, float(0)), countbookD.get(coord, float(0)), ccountbookD.get(coord, float(0)), float(100 * float(countbookD.get(coord, float(0)))/float(casecounter)), float(100 * float(ccountbookD.get(coord, float(0)))/float(controlcounter))]
    print('\t'.join(map(str,row)))
