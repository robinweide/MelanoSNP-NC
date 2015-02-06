#!/usr/bin/python

'''
to get all changed TFBS, use the following line:

python oPossumParser.py -i resFiles.lst -id y | sort -k2,2 -k4,4 -u -t";" | sed $'s/;/\t/g' | sed $'s/ /_/g' > all_non0-0_TFBS.lst

'''

import sys
import argparse
import os
from scipy import stats
import pprint

__author__ = 'Robin van der Weide'

# Command-line thingies
parser = argparse.ArgumentParser(
    description='Gets significant changes in motifs from oPossum-results.txtfile')
parser.add_argument('-i', '--input', help='List of results.txt path', required=True)
parser.add_argument('-fs', '--fisig', help='Set fisher significance-level (opossum)', default=float(100), required=False)
parser.add_argument('-ms', '--misig', help='Set fisher significance-level (TRIfle)', default=float(100), required=False)
parser.add_argument('-zs', '--Zsig', help='Set Z significance-level (opossum)', default=float(0.01), required=False)
parser.add_argument('-ks', '--Ksig', help='Set Ks significance-level (opossum)', default=float(0.01), required=False)
parser.add_argument('-Na', '--Ncase', help='Cases', default=int(1287), required=False)
parser.add_argument('-Nb', '--Ncontrol', help='controls', default=int(655), required=False)
parser.add_argument('-id', '--with_ID', help='display Ensembl-ids (y/n)', default=str("n"), required=False)

args = vars(parser.parse_args())

paths = open(args['input'], 'r')
Zsig = float(args['Zsig'])
Ksig = float(args['Ksig'])
sig = float(args['fisig'])
msig = float(args['misig'])
ensemblBiomartEnstEnsgName = open(os.path.dirname(os.path.abspath(__file__))+"/files/ensemblBiomartEnstEnsgName.tsv", 'r')

ensemblBiomartEnstEnsgNameDictEnsg = {} #enst | ensg
ensemblBiomartEnstEnsgNameDictName = {} # enst | name
for row in ensemblBiomartEnstEnsgName:
    t,g,n = row.split("\t")
    ensemblBiomartEnstEnsgNameDictEnsg[t] = g.rstrip()
    ensemblBiomartEnstEnsgNameDictName[t] = n.rstrip()

jasparFam = {}
jasparFam['MF0001.1'] = 'ETS_class'
jasparFam['MF0002.1 '] = 'bZIP_CREB/G-box-like_subclass'
jasparFam['MF0003.1'] = 'REL_class'
jasparFam['MF0004.1'] = 'Nuclear_Receptor_class'
jasparFam['MF0005.1'] = 'Forkhead_class'
jasparFam['MF0006.1'] = 'bZIP_cEBP-like_subclass'
jasparFam['MF0007.1'] = 'bHLH(zip)_class'
jasparFam['MF0008.1 '] = 'MADS_class'
jasparFam['MF0009.1'] = 'TRP(MYB)_class'
jasparFam['MF0010.1 '] = 'Homeobox_class'
jasparFam['MF0011.1'] = 'HMG_class'
for file in paths:
    fileName = os.path.dirname(file).split("/")[-1]
    filE = file.rstrip()
    entry = open(filE, 'r')
    for row in entry:
        if row.startswith('TF'):
            continue
        else:
            ID,TF,Class,Family,TaxGroup,IC,GCContent,Targetseqhits,Targetseqnonhits,Backgroundseqhits,Backgroundseqnonhits,TargetTFBShits,TargetTFBSnucleotiderate,BackgroundTFBShits,BackgroundTFBSnucleotiderate,Zscore,Fisherscore,KSscore = row.split("\t")
            print(fileName + "-" + TF)
            if float(Targetseqhits) + float(Backgroundseqhits) > float(0) and float(Targetseqnonhits) + float(Backgroundseqnonhits) > float(0):
                filee = open("./" + str(fileName) + "/" + str(TF) + ".txt.txt", 'r')
                for line in filee:
                    SeqID = ""
                    if not line.startswith("_"):
                        SeqID, Start, Stop, Strand, Score, pScore, Seq = line.split("\t")
                        Chr,coord = SeqID.split(":")
                        SeqStart, SeqStop = coord.split("-")
                        print(Chr + ":" + str(int(SeqStart)+int(Start)) + "-" + str(int(SeqStart) + int(Stop)))
                    else:
                        streep, Start, Stop, Strand, Score, pScore, Seq = line.split("\t")
                        Chr,coord = SeqID.split(":")
                        SeqStart, SeqStop = coord.split("-")
                        print(Chr + ":" + str(int(SeqStart)+int(Start)) + "-" + str(int(SeqStart) + int(Stop)))
                #print(ensemblBiomartEnstEnsgNameDictEnsg[fileName] + ";" + ensemblBiomartEnstEnsgNameDictName[fileName] + ";" + fileName + ";" + TF + ";" +Class + ";" +Family + ";" +str(Targetseqhits) + ";" + str(Targetseqnonhits) + ";" +str(Backgroundseqhits) + ";" + str(Backgroundseqnonhits) + ";" + str(Fisherscore) + ";" + str(Zscore)+ ";" + str(KSscore))
#            ratio = float(float(int(Targetseqhits)+int(Targetseqnonhits))/args['Ncase'])
#            thit = int(int(Targetseqhits)/ratio)
#            tbak = int(int(Targetseqnonhits)/ratio)
#            bhit = int(int(Backgroundseqhits)/ratio)
#            bbak = int( int(Backgroundseqnonhits)/ratio)
#            ratioCaseControl = float(float(args['Ncase'])/float(args['Ncontrol']))
#            GL = ""
#            if thit > (bhit*ratioCaseControl):
#                GL = "gain"
#            elif thit < (bhit*ratioCaseControl):
#                GL = "loss"
#            Moddsratio, Mpvalue = stats.fisher_exact([[thit, bhit], [tbak, bbak]])
#            pvalue = abs(float(Fisherscore))
#           #print(str(pvalue) + Fisherscore)
#            if Zscore == str("NA"):
#                Zscore = float(1)
#            #if abs(pvalue) <= float(sig) and abs(float(Zscore)) >= float(Zsig):
#            if KSscore.rstrip() == str("NA"):
#                continue
#            else:
#                if abs(float(KSscore)) <= float(Ksig) and float(0) <= abs(pvalue) <= float(sig) and abs(float(Zscore)) >= float(Zsig):
#                    if args['with_ID'] == "y":
#                        if TF in jasparFam:
#                        #print(str(thit) +str("\t") + str(bhit*ratioCaseControl) )
#                            print(ensemblBiomartEnstEnsgNameDictEnsg[fileName] + "\t" + ensemblBiomartEnstEnsgNameDictName[fileName] + "\t" + fileName + "\t" +jasparFam[TF] + "\t" +Class + "\t" +Family + "\t" +str(thit) + "\t" + str(tbak) + "\t" +str(bhit) + "\t" + str(bbak) + "\t" + str(GL) + "\t" + str(pvalue) + "\t" + str(Zscore)+ "\t" + str(KSscore))
#                        else:
#                            print(ensemblBiomartEnstEnsgNameDictEnsg[fileName] + "\t" + ensemblBiomartEnstEnsgNameDictName[fileName] + "\t" + fileName + "\t" + TF + "\t" +Class + "\t" +Family + "\t" +str(thit) + "\t" + str(tbak) + "\t" +str(bhit) + "\t" + str(bbak) + "\t" + str(GL) + "\t" + str(pvalue) + "\t" + str(Zscore)+ "\t" + str(KSscore))
#                    else:
#                        #print(str(thit) +str("\t") + str(bhit*ratioCaseControl) )
#                        if TF in jasparFam:
#                            print(ensemblBiomartEnstEnsgNameDictName[fileName] + "\t" + jasparFam[TF] + "\t" +Class + "\t" +Family + "\t" +str(thit) + "\t" + str(tbak) + "\t" +str(bhit) + "\t" + str(bbak) + "\t" + str(GL) + "\t" + str(pvalue)+ "\t" + str(Mpvalue))
#                        else:
#                            print(ensemblBiomartEnstEnsgNameDictName[fileName] + "\t" + TF + "\t" +Class + "\t" +Family + "\t" +str(thit) + "\t" + str(tbak) + "\t" +str(bhit) + "\t" + str(bbak) + "\t" + str(GL) + "\t" + str(pvalue)+ "\t" + str(Mpvalue))
