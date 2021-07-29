#!/usr/bin/env python
import sys
filename = sys.argv[1]
with open(filename, 'r') as infile:
    seqs = infile.read().split(">")[1:]
    for i in seqs:
        header = i.split()[0]
        seq = ''.join([x.strip() for x in i.split()[1:]])
        seqLen = len(seq)
        print("%s\t%s" % (header, seqLen))
