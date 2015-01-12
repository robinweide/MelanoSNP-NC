#!/usr/bin/python
import argparse
import re
import sys
from scipy import stats
import math
import pprint

__author__ = 'Robin van der Weide'

#Command-line thingies
parser = argparse.ArgumentParser(
    description='Takes expanded BED-files of case and control (with frequencies as kmers of a location) and outputs enriched windows.')
parser.add_argument('-a', '--cases', help='A sorted (sort -k 1,1n -k2,2n) and expanded BED of cases', required=True)
parser.add_argument('-r', '--ref', help='A sorted (sort -k 1,1n -k2,2n) and expanded BED of all possible start-locations of the window (e.g. reference, targets)', required=True)
parser.add_argument('-b', '--controls', help='A sorted (sort -k 1,1n -k2,2n) and expanded BED of controls', required=False)
parser.add_argument('-t', '--treshold', help='Percentage-difference treshold', required=False, default=0, type=float)
parser.add_argument('-Na', '--CasesN', help='counts of case samples', required=False, default=1, type=float)
parser.add_argument('-Nb', '--ControlsN', help='counts of control samples', required=False, default=1, type=float)
parser.add_argument('-w', '--window', help='Sliding-window size', required=True)
args = vars(parser.parse_args())

reflist = []
cases = []
controls = []
rawcasesWindowsDict ={}
rawcontrolsWindowsDict ={}
casesWindowsDict ={}
controlsWindowsDict ={}
windowSize = args['window']

allname = str("./") + str(int(windowSize)) + str("_window_all") +  str(".tsv")
allfile = open(allname,'w+')
chiname = str("./") + str(int(windowSize)) + str("_window_CHI2Y") + str(".tsv")
chifile = open(chiname,'w+')
chinameS = str("./") + str(int(windowSize)) + str("_window_CHI2Ys") + str(".tsv")
finame = str("./") + str(int(windowSize)) + str("_window_FISHER") + str(".tsv")



sys.stderr.write("Loading reference" + str("\n"))
rfile = open(args['ref'], 'r')
for row in rfile:
    chr,start,stop = re.split(r'\t+', row)
    item = str(chr) + str(":") + str(start)
    reflist.append(item)


sys.stderr.write("Loading cases" + str("\n"))
afile = open(args['cases'], 'r')
for row in afile:
    achr,astart,astop = re.split(r'\t+', row)
    item = str(achr) + str(":") + str(astart)
    cases.append(item)
sys.stderr.write("Loading controls" + str("\n"))
if args['controls'] is not None:
    bfile = open(args['controls'], 'r')
    for row in bfile:
        bchr,bstart,bstop = re.split(r'\t+', row)
        item = str(bchr) + str(":") + str(bstart)
        controls.append(item)


sys.stderr.write("Running sliding window:\n")
perC = float(0)
locCounta = float(0)
locCountb = float(0)
refcount = int(0)
for coord in reflist:
    locCounta = float(0)
    locCountb = float(0)
    refcount += int(1)
    chr,start = coord.split(':')
    chR = chr.upper()
    if re.match('^CHR', chR):
        CHR = float(chr[3:])
    else:
        CHR = float(chr)
    MIN = float(start)
    MAX = float(start) + float(windowSize)
    entry = str(chr) + str(":") + str(int(MIN)) + str("-") + str(int(MAX))
    casesWindowsDict[entry] = float(0)
    controlsWindowsDict[entry] = float(0)
    leng = len(reflist)
    perc = 100 * float(refcount)/float(leng)
    if perc >= perC:
        #sys.stderr.write(str("\t") + str(perC) + str("%"))
        sys.stderr.write("\r%    d%%    " % perc)
        sys.stderr.flush()
        perC += float(0.5)
    if float(len(cases)) > float(0):
        testcase = cases[0]
    else:
        sys.stderr.write("\r%    d%%    " % int(100))
        sys.stderr.flush()
        break
    c,l = testcase.split(':')
    C = c.upper()
    if re.match('^CHR', C):
        snpCc = c[3:]
        cc = float(snpCc)
    else:
        cc = float(c)
    if CHR < float(cc):
        continue
    elif CHR == float(cc) and MAX < float(l):
        continue
    else:
        for snp in cases:
            sc,ss = snp.split(':')
            c = sc.upper()
            if re.match('^CHR', c):
                snpCc = c[3:]
                snpC = float(snpCc)
            else:
                snpC = float(c)
            snpS = float(ss)
            #print(str(snpC) + str(":") + str(int(snpS)))
            if CHR == snpC and MIN <= snpS <= MAX:
                locCounta += float(1)
                #casesWindowsDict[entry] += float(1)
                #print(str(casesWindowsDict[entry]) + str("\t") + snp)
            if CHR >= snpC and snpS < MIN:
                cases.remove(snp)
            elif CHR > snpC:
                cases.remove(snp)
            if snpC >= CHR and snpS > MAX:
                break
        if args['controls'] is not None:
            for snp in controls:
                sc,ss = snp.split(':')
                c = sc.upper()
                if re.match('^CHR', c):
                    snpC = float(c[3:])
                else:
                    snpC = float(c)
                snpS = float(ss)
                if CHR == snpC and MIN <= snpS <= MAX:
                    locCountb += float(1)
                if CHR >= snpC and snpS < MIN:
                    controls.remove(snp)
                elif CHR > snpC:
                    controls.remove(snp)
                if snpC >= CHR and snpS > MAX:
                    break

    if args['controls'] is not None and locCountb > float(0):
        controlsWindowsDict[entry] = (locCountb/float(args['ControlsN']))
        rawcontrolsWindowsDict[entry] = locCountb
    else:
        controlsWindowsDict[entry] = float(0)
        rawcontrolsWindowsDict[entry] = float(0)
    if locCounta > float(0):
        casesWindowsDict[entry] = (locCounta/args['CasesN'])
        rawcasesWindowsDict[entry] = locCounta
    else:
        rawcasesWindowsDict[entry] = float(0)
sys.stderr.write("\nFinding and saving relevant windows" + str("\n"))
fifile = open(finame,'w+')
for k,v in rawcasesWindowsDict.items():
    chr,loc = k.split(':')
    star,stop = loc.split('-')
    outputLine = str(str(chr)  + "\t" + str(star)  + "\t" + str(stop))
    OM = float(v+0.0000001)
    EM = float(rawcontrolsWindowsDict[k]+0.0000001)
    OWT = float(float(args['CasesN']*float(windowSize))-OM)
    EWT = float(float(args['ControlsN']*float(windowSize))-EM)
    EM1 = float((EM/float(args['ControlsN']))*float(args['CasesN']))
    EWT1 = float((EWT/float(args['ControlsN']))*float(args['CasesN']))
    CHI = float(((((OM - EM1)-0.5)**2)/EM1) + ((((OWT - EWT1)-0.5)**2)/EWT1))
    if CHI > float(3.841):
        chifile.write(outputLine + "\t" + str(OM-0.0000001) + "\t" + str(EM-0.0000001) + "\t" + str((OM-0.0000001)/(EM1-0.0000001)) + "\t" + str(CHI) + "\n")
    allfile.write(outputLine + "\t" + str(OM-0.0000001) + "\t" + str(EM-0.0000001) + "\n")
    oddsratio, pvalue = stats.fisher_exact([[OM, EM1], [OWT, EWT1]])
    if pvalue <= float(0.05):
        fifile.write(outputLine + "\t" + str(OM-0.0000001) + "\t" + str(EM-0.0000001) + "\t" + str((OM-0.0000001)/(EM1-0.0000001)) + "\t" + str(pvalue) + "\n")
sys.stderr.write("DONE!!!" + str("\n"))



chifile.close()

sys.stderr.write("Now do the following to steps to get rid of overlaps:" + str("\n"))
sys.stderr.write("||| \t cat " + chiname + "  | sed $'s/:/\t/' | sed $'s/-/\t/' | sort -k1,1n -k2,2n > tmp" + str("\t|||\n"))
sys.stderr.write("||| \t bedtools merge -i tmp -c 4,5,6,7 -o mean,mean,mean,mean -d 20 | sort -k6,6n > " + chiname + ".mergedNsorted.tsv ; rm tmp" + str("\t|||\n"))
sys.stderr.write("ENJOY :)" + str("\n"))


