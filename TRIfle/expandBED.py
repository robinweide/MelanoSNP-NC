#!/usr/bin/python
import argparse
import re

__author__ = 'Robin van der Weide'

# Command-line thingies
parser = argparse.ArgumentParser(
    description='Expands BED-file to per-base lines')
parser.add_argument('-i', '--inputBED', help='BED-file in need of expanding', required=True)
args = vars(parser.parse_args())

ifile = open(args['inputBED'], 'r')
for row in ifile:
    chr, start, stop = re.split(r'\t+', row)
    cunter = int(start)
    while (cunter <= int(stop)):
        print(str(chr) + str("\t") + str(cunter) + str("\t") + str(cunter + 1))
        cunter += 1