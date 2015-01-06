#!/usr/bin/python
import argparse
import re
import sys
import math
import pprint

__author__ = 'Robin van der Weide'

#Command-line thingies
parser = argparse.ArgumentParser(
    description='Takes expanded BED-files of case and control (with frequencies as kmers of a location) and outputs enriched windows.')
parser.add_argument('-a', '--cases', help='A sorted (sort -k 1,1n -k2,2n) and expanded BED of cases', required=True)
parser.add_argument('-r', '--ref', help='A sorted (sort -k 1,1n -k2,2n) and expanded BED of all possible start-locations of the window (e.g. reference, targets)', required=True)
parser.add_argument('-b', '--controls', help='A sorted (sort -k 1,1n -k2,2n) and expanded BED of controls', required=False)
parser.add_argument('-t', '--treshold', help='Difference treshold', required=False, default=1, type=float)
parser.add_argument('-Na', '--CasesN', help='counts of case samples', required=False, default=1, type=float)
parser.add_argument('-Nb', '--ControlsN', help='counts of control samples', required=False, default=1, type=float)
parser.add_argument('-w', '--window', help='Sliding-window size', required=True)
args = vars(parser.parse_args())

reflist = []
cases = []
controls = []
casesWindowsDict ={}
controlsWindowsDict ={}
windowSize = args['window']

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



perC = int(0)
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
        sys.stderr.write("Progress: " + str(perC) + str("%\n"))
        perC += 10
    if float(len(cases)) > float(0):
        testcase = cases[0]
    else:
        sys.stderr.write("Progress: " + str("100") + str("%\n"))
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
                    #controlsWindowsDict[entry] += float(1)
    #print(locCountb)
    #print(locCounta)
    if args['controls'] is not None and locCountb > float(0):
        controlsWindowsDict[entry] = math.log10(locCountb/float(windowSize))/float(args['ControlsN'])
    else:
        controlsWindowsDict[entry] = float(0)
    if locCounta > float(0):
        casesWindowsDict[entry] = math.log10(locCounta/float(windowSize))/float(args['CasesN'])
    else:
        casesWindowsDict[entry] = float(0)
sys.stderr.write("Finding and saving relevant windows" + str("\n"))
for k,v in casesWindowsDict.items():
    if controlsWindowsDict[k] < v:
        print(str(k) + "\t" + str(v) + "\t" + str(controlsWindowsDict[k]))
#pprint.pprint(controlsWindowsDict)
sys.stderr.write("DONE!!!" + str("\n"))


#for window,freq in casesWindowsDict.items():
#    print(str(window) + str("\t") + str(freq) + str("\t") + str(controlsWindowsDict[window]))
    #if float(controlsWindowsDict[window] + args['treshold']) <= float(freq) <= float(controlsWindowsDict[window] - args['treshold']):
    #    print(str(window) + str("\t") + str(freq) + str("\t") + str(controlsWindowsDict[window]))