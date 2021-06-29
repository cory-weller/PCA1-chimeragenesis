#!/usr/bin/env py

import sys

infileName = sys.argv[1]

# 789 nucleotides of manually edited sequence
editedNucleotides =  """ATGAAACCTAAGAAGTTTGATTTCTCCCTTGGAACTTCAGACAATGATAGGAAGGCTGGCAACTCTGAAA\
ATATTCCTATTGACACGATGTATGGCTCCAATTGGCCTGTGGAGATGCACGCCTCCGATGACCAACGCCT\
ATCAAGTCGCAAACAAAAGTGCCTAAACAAATCCGAAGTTATGTGTAACGAAGGAGGTAACCGAAACAAA\
GCTGCCTCAGACTCAGACAGCTGTTGCGGTGATGCCCAGCGAGGCGAAAATTACAGTGACGAGTCTTGTG\
TAGACAAGTGCTGTGCTGAGAAGGAGAATGAAACGGAAGCCGCGTTTGATTCCGATAGTTGCTGTGGCGA\
TGCTCAACGCGGGGAAAACTATAGTGATGAATCATGTGTCAACGAGTGTTGTGCAAAGAAAGAGAACGAG\
ACTGAGGCGGCATCGGACTCGGACTCCTGTTGTGGTGACGCTCAACGTGGAGAGAATTACAGTGATGAAA\
GCTGTGTAGATAAGTGTTGCGCCGAAAAGGAGAACGAAACCGAAGCAGCTTTCGACTCTGACTCGTGCTG\
CGGGGATGCTCAGCGTGGCGAGAACTATTCGGACGAGTCATGTGTTAATGAATGCTGTGCGAAGAAAGAA\
AATGAGACCGAGGCCGCTTCAGACAGCGATTCGTGTTGCGGAGACGCGCAGCGCGGTGAAAACTACTCAG\
ATGAGTCATGCGTAAATGAGTGCTGCGCTAAGAAGGAAAACGAGACAGAAGCGGCTAGCGGTAGTGACTC\
ATGCTGTGGCGATGCACAG"""

with open(infileName, 'r') as infile:
    text = [x.strip() for x in infile.readlines()]
    header = text[0]
    seq = ''.join(text[1:]).replace("-","")

header = header.split("|")[0]
header = header + "|length=" + str(len(seq))


newSeq = editedNucleotides + seq[789:]

def foldSeq(seq, lineLength):
    return('\n'.join([seq[x:(x+lineLength)] for x in range(0,len(seq), lineLength)]))

assert len(seq) == len(newSeq), "sequence lengths do not match"

print(header + "\n" + foldSeq(newSeq, 60))