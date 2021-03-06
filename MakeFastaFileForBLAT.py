#! /usr/bin/python

# Purpose:  To create the fasta file for the subsequent BLAT queries.
# Input:    1. A list of Illumina IDs for SNPs
#           2. A bpm.csv file
# Output:   test.out.fa file -> we can split this up into pieces and submit it to BLAT on the cluster
# 
# Author: et85, etsai@bwh.harvard.edu

from argparse import ArgumentParser
import re, sys

def checkID(thisID, thisDict):
    if(thisDict.has_key(thisID)):
        return 1
    else:
        return 0

def writeFasta(thisIdInd,thisSeqInd,thisList,thisOutFH):
#     print thisList[thisIdInd] + " " + thisList[thisSeqInd]
    seqList = re.split('\[|\/|\]', thisList[thisSeqInd])
    nameSeqA = thisList[thisIdInd] + '-' + seqList[1]
    nameSeqB = thisList[thisIdInd] + '-' + seqList[2]
    
    if(seqList[1] == '-'):
        seqA = seqList[0] + seqList[3]
        seqB = seqList[0] + seqList[2] + seqList[3]
    elif(seqList[2] == '-'):
        seqA = seqList[0] + seqList[1] + seqList[3]
        seqB = seqList[0] + seqList[3]
    else:
        seqA = seqList[0] + seqList[1] + seqList[3]
        seqB = seqList[0] + seqList[2] + seqList[3]
#     print nameSeqA + " " + seqA
#     print nameSeqB + " " + seqB
    
    thisOutFH.write('> ' + nameSeqA + '\n')
    thisOutFH.write(seqA + '\n')
    thisOutFH.write('> ' + nameSeqB + '\n')
    thisOutFH.write(seqB + '\n')

parser = ArgumentParser()

parser.add_argument("-B","--bpm",
                    type=str,
                    dest="bpm",
                    default='/Users/et85/testing/gtc/MEGA_Consortium_15063755_B1.csv',
                    help="csv-version of bpm file")

parser.add_argument("-I","--infile",
                    type=str,
                    dest="infile",
                    default='/Users/et85/testing/fwdstrand/allissues_x25145.txt',
                    help="input list of SNP names to format for the the output fasta file")

parser.add_argument("-O","--outfile",
                    type=str,
                    dest="outfile",
                    default='/Users/et85/testing/fwdstrand/test.out.fa',
                    help="output fasta file")

args = parser.parse_args()

if(args.bpm):
    thisBPM = args.bpm

if(args.infile):
    infile = args.infile

if(args.outfile):
    outfile = args.outfile

bpmFH = open(thisBPM,"r")
inFH = open(infile,"r")
outFH = open(outfile,"w")

# Save infile into hash

nameDict = {}
for thisLine in inFH.readlines():
    thisLine = thisLine.rstrip()
    nameDict[thisLine] = 1

# Header line right after '[Assay]' line
flagAssay = 0
idInd = nameInd = seqInd = ''

for thisLine in bpmFH.readlines():
    thisLine = thisLine.rstrip()
    if(flagAssay == 2):
        thisList = thisLine.split(',')
        if(checkID(thisList[idInd],nameDict)):
            writeFasta(idInd,seqInd,thisList,outFH)
    elif(flagAssay == 1):
        thisList = thisLine.split(',')
        idInd = thisList.index('IlmnID')
        nameInd = thisList.index('Name')
        seqInd = thisList.index('SourceSeq')
        flagAssay = 2
    elif(thisLine == '[Assay]'):
        flagAssay = 1
    elif(thisLine == '[Controls]'):  # end of for loop
        break











