#!/usr/bin/python
import argparse
import re
from subprocess import PIPE,Popen, STDOUT
__author__ = 'Robin van der Weide'

#Command-line thingies
parser = argparse.ArgumentParser(
    description='This is TRIfle: integrated non-coding variant prioritisation. Please give ONE input-file.')
parser.add_argument('-v', '--vcf', help='TRIfle-filter output', required=True)
parser.add_argument('-r', '--ref', help='path to reference genome', required=True)
parser.add_argument('-o', '--out', help='output-name', required=True)
args = vars(parser.parse_args())
def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]
#Main script

filename = str("./") + str(args['out']) + str("_WtMt.fasta")
fastafile = open(filename,'w+')
filename = str("./") + str(args['out']) + str("_WtMt.ids")
idsfile = open(filename,'w+')

ifile = open(args['vcf'], 'r')
for row in ifile:
    fields = re.split(r'\t+', row)
    coords = str(fields[0]) + str(":") + str(int(fields[1])-29) + str("-") + str(int(fields[1])+29)
    mut = str(fields[4].lower())
    p = cmdline('samtools faidx ' + args['ref'] + " chr" + coords + " | tail -n1").strip().lower()
    seqW = list(p)
    seqM = list(p)
    seqM[29] = mut.upper()
    fastafile.write(">"+str(fields[0]) + str(":") + str(fields[1]) + str("-") + str(int(fields[1]) + 1) + str("|") + str(fields[3]) + "/" + str(fields[4]) + "$" + str(fields[6]) + "_WT\n")
    fastafile.write("".join(seqW).upper() + "\n")
    fastafile.write(">"+str(fields[0]) + str(":") + str(fields[1]) + str("-") + str(int(fields[1]) + 1) + str("|") + str(fields[3]) + "/" + str(fields[4]) + "$" + str(fields[6]) + "_MT\n")
    fastafile.write("".join(seqM).upper() +"\n")
    idsfile.write(str(fields[0]) + str(":") + str(fields[1]) + str("-") + str(int(fields[1]) + 1) + str("|") + str(fields[3]) + "/" + str(fields[4]) + "$" + str(fields[6]) +"_WT\t" +str(fields[0]) + str(":") + str(fields[1]) + str("-") + str(int(fields[1]) + 1) + str("|") + str(fields[3]) + "/" + str(fields[4]) + "$" + str(fields[6]) +"_MT\n")
ifile.close()

