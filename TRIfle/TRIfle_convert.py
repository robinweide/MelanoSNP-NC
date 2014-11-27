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
parser.add_argument('-s', '--SuRFR', help='input: list of output.tsv', required=False)
#parser.add_argument('-cc','--Case-Control', help='Case(1) or Control(0)', required=True, default='1')
args = vars(parser.parse_args())
if not any(args.values()):
    parser.error("Error: EITHER -f, -s or -c must be given. Use -h for help.")
if len(sys.argv) > 3:
    parser.error("Error: only ONE of the options (-f, -s or -c) must be given. Use -h for help.")

#FUNCTIONS

def saveDict(fn,dict_rap):
    f=open(fn, "wb")
    w = csv.writer(f)
    for key, val in dict_rap.items():
        w.writerow([key, val])
    f.close()
def readDict(fn):
    f=open(fn,'rb')
    dict_rap={}

    for key, val in csv.reader(f):
        dict_rap[key]=eval(val)
    f.close()
    return(dict_rap)
def funseq2scorecount(Ffile, scorebookN, scorebookC, countbook):
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
                coord = str(coord + "*C")
                if coord in scorebookC:
                    if float(scorebookC[coord]) < float(score):
                        scorebookC[coord] = float(score)
                else:
                    scorebookC[coord] = float(score)
            else:
                index = row.find("NCDS", 40)
                if not index == -1:
                    inforow = row[index:]
                    CDSSfield = re.split(r'[;]', inforow)[0]
                    score = float(re.split(r'[=]+', CDSSfield)[-1].rstrip())
                    coord = str(coord + "*N")
                    if coord in scorebookN:
                        if float(scorebookN[coord]) < float(score):
                            scorebookN[coord] = float(score)
                    else:
                        scorebookN[coord] = float(score)
            if coord in countbook:
                countbook[coord] += 1
            else:
                countbook[coord] = float(1)
    return(scorebookN, scorebookC, countbook)
def funseqzscore(scorebookN, scorebookC):
    mean=numpy.mean(scorebookN.values())
    standardDeviation= numpy.std(scorebookN.values())
    a={}
    for key,value in scorebookN.items():
        z=(value-mean)/standardDeviation
        ZscorebookN[key]=z
    mean=numpy.mean(scorebookC.values())
    standardDeviation= numpy.std(scorebookC.values())
    a={}
    for key,value in scorebookC.items():
        z=(value-mean)/standardDeviation
        ZscorebookC[key]=z
    scorebook = ZscorebookN.copy()
    scorebook.update(ZscorebookC)
    return scorebook

#General declarations
scorebook = {}
scorebookF = {}
countbookF = {}
countbook = {}
scorebookN = {}
scorebookC = {}
ZscorebookN = {}
ZscorebookC = {}

#Main script
if args['Funseq2'] is not None:
    Ffile = open(args['Funseq2'], 'r')
    for row in Ffile:
        frow = open(row.rstrip(), 'r')
        scorebookN, scorebookC, countbookF = funseq2scorecount(frow, scorebookN, scorebookC, countbookF)
        frow.close()
    Ffile.close()
elif args['CADD'] is not None:
    Cfile = open(args['CADD'], 'r')
elif args['SuRFR'] is not None:
    Sfile = open(args['SuRFR'], 'r')

# CHECK FOR MISSING VARIANTS BY COMPARING THE COORDS IN THE FILES WITH THE



print("Scorebook:")
pprint(scorebookF)
print("Countbook:")
pprint(countbookF)
print("LET OP: ALS ALLE DRIE ALLE VARIANTEN ALS OUTPUT GEVEN, DAN KUNNEN DE COUNTBOOKS ALLEEN OP FUNSEQ GEBASEERD WORDEN!!!")
print("LET OP: MERGE DE COORDS VAN FUNSEQ-C EN -N! ZODAT ELKE COORD (CHRLPOS|NUC) IS")
