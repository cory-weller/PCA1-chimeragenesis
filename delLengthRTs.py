#!/usr/bin/env python3

cutSite = 5821

HAlength = 80

with open('seqs/PW5_min_construct.fasta', 'r') as infile:
    seq = ''.join([x.strip() for x in infile.readlines()[1:]])

breaks = [25] * 9 +  [100] * 15 + [250]*12

newBreaks = []
cumSum = 275


for i in breaks:
    cumSum += i
    newBreaks.append(cumSum)

leftArms = []
rightArms = []

for i in newBreaks:
    leftEnd = cutSite - i
    leftStart = leftEnd - HAlength
    rightStart = cutSite + i
    rightEnd = rightStart + HAlength
    leftArm = seq[leftStart : leftEnd].lower()
    rightArm = seq[rightStart : rightEnd].upper()
    leftArms.append(("-%s" % i, leftArm))
    rightArms.append(("+%s" % i, rightArm))

for L in leftArms:
    for R in rightArms:
        print(">%s:%s\n%s\n%s" % (L[0], R[0], L[1], R[1]))