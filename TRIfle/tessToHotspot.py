#!/usr/bin/python
import argparse
import itertools
__author__ = 'Robin van der Weide'

#Command-line thingies
parser = argparse.ArgumentParser(
    description='Takes tess bulk-file and prints a list of consequences (coord | tfbs | consequence).\n Remember: it prints each consequence Ncases times.')
parser.add_argument('-i', '--input', help='TESS bulk-file', required=True)
parser.add_argument('-t', '--treshold', help='Delta-OR treshold', required=False, default=1, type=float)
args = vars(parser.parse_args())

#Main script
#hotspotfile = open(filename,'w+')
#detailsfile = open(details.tsv,'w+')
WTdict = {}
MTdict = {}
coordCounter = {}

ifile = open(args['input'], 'r')
for row in ifile:
    MatrixId,SequenceID,La,HitBegin,HitEnd,Sense = row.split(',')
    coord,info = SequenceID.split('|')
    mut,ranking = info.split('$')
    r,f,x = ranking.split(';')
    y,genotype = x.split('_')
    x,rank = r.split('=')
    x,freq = f.split('=')
    dictEntry = str(coord) + str("|") + str(MatrixId)
    coordCounter[coord] =+ int(freq)
    if genotype == "MT":
        if dictEntry in MTdict:
            if MTdict[dictEntry] > La:
                continue
            else:
                MTdict[dictEntry] = La
        else:
            MTdict[dictEntry] = La
    if genotype == "WT":
        if dictEntry in WTdict:
            if WTdict[dictEntry] > La:
                continue
            else:
                WTdict[dictEntry] = La
        else:
            WTdict[dictEntry] = La

    '''
    Fields now available:
    MatrixId,SequenceID,La,HitBegin,HitEnd,Sense,genotype,coord,mut,rank,freq(case)
    '''
ifile.close()

for entry,value in WTdict.items():
    coord,tfbs = entry.split('|')
    if entry in MTdict:
        if value == MTdict[entry]:
            continue
        else:
            if float(value) > (float(MTdict[entry]) + float(args['treshold'])):
                for _ in itertools.repeat(None, coordCounter[coord]):
                    print(str(coord) + str("\t") + str(tfbs) + str("\t") + str("LOSS"))
            elif float(value) < (float(MTdict[entry]) - float(args['treshold'])):
                for _ in itertools.repeat(None, coordCounter[coord]):
                    print(str(coord) + str("\t") + str(tfbs) + str("\t") + str("GAIN"))
    else:
        for _ in itertools.repeat(None, coordCounter[coord]):
            print(str(coord) + str("\t") + str(tfbs) + str("\t") + str("LOSS"))

for entry,value in MTdict.items():
    coord,tfbs = entry.split('|')
    if entry in WTdict:
        continue
    else:
        for _ in itertools.repeat(None, coordCounter[coord]):
            print(str(coord) + str("\t") + str(tfbs) + str("\t") + str("GAIN"))


'''
detailsfile:
coord   TFBS    loss/gain (tov wt)
'''
