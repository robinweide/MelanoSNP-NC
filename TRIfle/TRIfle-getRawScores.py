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
variantCoordList = []

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
        crow.close()
    Cfile.close()



# Adding zeros to missing scores


coords = list(set(variantCoordList))
coords.sort()
for p in coords:
    print(p)


#if scorebookFN.keys() == scorebookFC.keys():
#    print("FN=FC")
#if scorebookFN.keys() == scorebookC.keys():
#    print("FN=C")
#if scorebookFN.keys() == scorebookD.keys():
#    print("FN=D")
#if scorebookFC.keys() == scorebookFN.keys():
#    print("FC=FN")
#if scorebookFC.keys() == scorebookC.keys():
#    print("FC=C")
#if scorebookFC.keys() == scorebookD.keys():
#    print("FC=D")
#if scorebookD.keys() == scorebookFC.keys():
#    print("D=FC")
#if scorebookD.keys() == scorebookC.keys():
#    print("D=C")
#if scorebookD.keys() == scorebookFN.keys():
#    print("D=FN")
#if scorebookC.keys() == scorebookFC.keys():
#    print("C=FC")
#if scorebookC.keys() == scorebookD.keys():
#    print("C=D")
#if scorebookC.keys() == scorebookFN.keys():
#    print("DC=FN")


#print tab-delim: coord:FC:FN:C:Count
#header = ["#Coord|mut","Funseq2(Nc)","Funseq2(C)","CADD(Phred)","DANN","Frequency"]
#print('\t'.join(map(str,header)))



for coord in coords:
    row = [scorebookFC.get(coord, float(0))]
    print('\t'.join(map(str,row)))
#for coord,score in scorebookD.items():
#    row = [coord, score, scorebookFC[coord],scorebookC[coord],scorebookD[coord],countbookC[coord]]
#    print('\t'.join(map(str,row)))