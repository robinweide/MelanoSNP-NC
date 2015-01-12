#!/usr/bin/python
import argparse
import re
import os

__author__ = 'Robin van der Weide'
# Command-line thingies
parser = argparse.ArgumentParser(
    description='Expands BED-file to per-base lines')
parser.add_argument('-i', '--inputBED', help='BED-file with motifs', required=True)
parser.add_argument('-m', '--motifName', help='TSV with jaspar-code and gene names', required=False, default=os.path.dirname(os.path.abspath(__file__))+ "/files/JASPAR_CORE_2014.ids") #this file has: jaspar | genename
parser.add_argument('-g', '--go', help='File with GO-annotations of Homo sapiens (reference)', required=False, default=os.path.dirname(os.path.abspath(__file__))+ "/files/gene_association.goa_ref_human") #this file has: col3 (name) and col4(GO)
#parser.add_argument('-f', '--families', help='TSV with gene names and families', required=False, default=os.path.dirname(os.path.abspath(__file__))+ "/files/families.tsv") #this file has: genename | familyname
#parser.add_argument('-t', '--targets', help='TSV with gene names and targets from PAZAR', required=False, default=os.path.dirname(os.path.abspath(__file__))+ "/files/targets.tsv") #this file has: genename | targetname
args = vars(parser.parse_args())

inFile = open(args['inputBED'], 'r')
motifFile = open(args['motifName'], 'r')
goFile = open(args['go'], 'r')
#targetFile = open(args['targets'], 'r')

'''
Make dict of gofile:
col3 -> col4 (name->go)
if col3 exist, append Go-term to go (delim=,)
'''
goDict = {}
for row in goFile:
    row.rstrip()
    if row.startswith('!'):
            continue
    else:
        p = row.split("\t")
        if p[2] in goDict:
            goDict[p[2]].append(p[4])
        else:
            goDict[p[2]] = [p[4]]


'''
make dict of motifFile
col1 -> col2 (jaspar/name)
'''
MotifUniprotDict = {}
for row in motifFile:
    motif, uniprot = re.split(r'\t+', row)
    MotifUniprotDict[motif] = [uniprot.upper().rstrip()]

'''
make dict from input
window -> all uniprotnames (delim = ,)


windowMotifDict = {}
windowInfoDict = {}
for row in inFile:
    chr, start, stop, caseSNPs, controlSNPs,fisher, encodeTFregion, motif, gainloss = re.split(r'\t+', row)
    coord = str(str(chr)  + ":" + str(start)  + "-" + str(stop))
    info = str(str(caseSNPs)  + "/" + str(controlSNPs)  + ":" + str(fisher)+ "|" + str(encodeTFregion))
    windowInfoDict[coord] = info
    uni = MotifUniprotDict[motif]
    if coord in windowMotifDict:
        windowMotifDict[coord].append(uni)
    else:
        windowMotifDict[coord] = [uni]
#for k,v in goDict.items():
#    print(str(k) + "\t" + str(v))



Get Go terms of uniprotnames
If TF-names has ::*, split and try both.
'''
'''
golist = []
goList = []
longlist = []
nameList = []
for coord,motiflist in windowMotifDict.items():
    for motif in motiflist:
        if motif[0].find("::") == -1:
            if motif[0] in goDict:
                go = goDict[motif[0]]
                golist.append(go)
                longlist.append(str(coord) + str("\t") + str(windowInfoDict[coord]) + str("\t") + str(motif[0]))
                nameList.append(motif[0])
        else:
            p = motif[0].split("::")
            for entry in p:
                #print(entry)
                if entry in goDict:
                    go = goDict[entry]
                    golist.append(go)
                    nameList.append(entry)
                    longlist.append(str(coord) + str("\t") + str(windowInfoDict[coord]) + str("\t") + str(entry))
    for blok in golist:
        for entry in blok:
            goList.append(entry)

    #uncomment following line for  printing     9:21935019-21935035	0.0/4.0:0.0156165|0	[ GABPA    CEBPA    GFI1B ]

    #print(str(coord) + str("\t") + str(windowInfoDict[coord]) + str("\t") + str(nameList) )
    nameList = []
    golist = []
    goList = []'''
#'''
#Convert TF-names to families.
#If TF-names has ::*, split and try both.
#'''

#'''
#Get Go terms of families.
#'''

'''
Get target-genes of TFs.
'''

'''
Get go terms of targets
'''

#'''
#Get target-genes of families.
#'''

#'''
#Get go terms of targets
#'''

'''
for each line:
    print coord info tf gainloss
'''

for row in inFile:
    chr, start, stop, caseSNPs, controlSNPs,fisher, encodeTFregion, motif, gainloss = re.split(r'\t+', row)
    coord = str(str(chr)  + "\t" + str(start)  + "\t" + str(stop))
    blok = MotifUniprotDict[motif]
    gene = ""
    for i in blok:
        gene = i
    Info = str(str(caseSNPs)  + "\t" + str(controlSNPs)  + "\t" + str(fisher) + "\t" + str(encodeTFregion) + "\t" + str(gene)+ "\t" + str(gainloss))
    print(coord + "\t" + Info.rstrip())
#for k,v in goDict.items():
#    print(str(k) + "\t" + str(v))