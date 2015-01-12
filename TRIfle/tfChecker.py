__author__ = 'robin'
import sys
import re
for line in sys.stdin:
    chr, start, stop, case, control, fisher, encode, foundGene, gainLoss, c,a,b,x,n,tfrow,q = re.split(r'\t+', line)
    tfs = tfrow.split(",")
    for tf in tfs:
        posstf = tf.split("(")
        if posstf[0].upper() == foundGene.upper():
            print(chr + str("\t") + start + str("\t") + stop + str("\t") + posstf[0] + str("\t") + foundGene + str("\t") + gainLoss)

'''
DIT WERKT, EN MOET GEBRUIKT WORDEN, MAAR IS AAN DE  KLEINE KANT.
GEBRUIK HIERVOOR WEL DE ENSEMBL IDS VOOR DE TF'S... ANDERS WORDT HET SNEL ERG ROMMELIG
PROBEER DE LIJST VAN HT.F.ASS OM GENEN NAAR FAMILIES OM TE ZETTEN (http://www.transcriptionfactor.org/index.cgi?Download)
DAN TESTEN FAM-FAM
PRINT DAN ZOWEL GENEN ALS FAMS
'''