#! /usr/bin/python
# Modified from zCall code by Jackie Goldstein by Jeremy Peirce
# Removed rare variant calling code so this now takes a GTC and BPM file (newer format eg CytoSNP 850K) and outputs PED

import sys
import re
from BPMjmod import *
from SAM import *  #you can import from the "SAM" or "BLAT" modules
#from BLAT import *
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-B","--bpm",type="string",
                  dest="bpm",action="store",
                  help="CSV format of BPM file")
parser.add_option("-A","--blat",type="string",dest="blat",action="store",help="BLAT alignments for SNPs")
parser.add_option("-O","--output",type="string",dest="output",action="store",
                  help="specify output format for further alignment work as text, fasta, or doublefasta, or as a strandlist of matched strands (name, strand) for gtc2ped.  Specifying diagnostics will generate summary statistics")

(options, args) = parser.parse_args()

if options.blat == None:
    print "specify blat alignment file with -A"

if options.bpm == None:
    print "specify BPM file path with -B"
    sys.exit()

if options.output == None:
    print "specity output format with -O"
    sys.exit

### Initialize BLAT and BPM file classes
bpm = BPM(options.bpm)
blat = BLAT(options.blat)
#print blat.names[0], blat.chr[0],blat.pos[0],blat.strand[0]

### Get Number of SNPs and blat alignments
numSNPs = len(bpm.names)
numBlatAlignments = len(blat.names)
if options.output == "diagnostics":
    print 'Diagnostics Mode.  Read BPM and Alignment File'
    print 'There are ', numSNPs, ' SNPs'
    print 'There are ', numBlatAlignments, ' BLAT alignments'

### Blat results Dictionary Setup
snpLookup={}

### Initialize dictionary with an empty list for each value to append later
for name in bpm.names:
    snpLookup.setdefault(name, [])

### Set up dictionary with key=SNP name, values will be blat entries in the format
### (chr, pos, strand) for each entry in a list    
for i in xrange(numBlatAlignments):
    snpLookup[blat.names[i]].append((blat.chr[i], blat.pos[i], blat.strand[i]))

if options.output == "diagnostics":
    print 'created dictionary'

"""For each SNP name in bpm find the snp name in the dictionary of blat results
Check the chromosome, then position in the blat results.
If they both match, position to within 50bp, then write the strand.
If no match check the next blat result.
If no more blat results write the strand as a u """

strandList = []
numMatchedPositions = 0
numUnmatchedPositions = 0
numUniqueUnmatchedPositions = 0
numZeroBpm = 0
numForward = 0
numReverse = 0
numUnknown = 0

for i in xrange(numSNPs): #numSNPs
  blatEntries = snpLookup[bpm.names[i]]  # gets the list of blat results for the ith SNP name
  numEntries = len(blatEntries) # how many entries for this snp
  for j in xrange(numEntries):
      noMatch = 0
      blatChr = blatEntries[j][0]
      blatPos = blatEntries[j][1]
      
      #handle mismatch in reporting of XY region
      if bpm.chr[i] == "XY":
          if blatChr == "X" or blatChr == "Y":
              blatChr = "XY"
              
      #check whether chromosome is the same and position are close
      if blatChr == bpm.chr[i] and (abs(int(blatPos) - int(bpm.pos[i]))<1000):
              thisBlatStrand = blatEntries[j][2] #forward/reverse in +/- notation
              strandList.append(thisBlatStrand) #makes a list of blat strands when one is found

              #gather forward/reverse stats
              if blatEntries[j][2] == "+":
                  numForward += 1
              elif blatEntries[j][2] == "-":
                  numReverse += 1
              elif blatEntries[j][2] == "u":
                  numUnknown += 1
              numMatchedPositions +=1
              
              break           
      else:
          noMatch = 1
          

          
  #if this is still true after the end of checking multiple alignments add 1 to unmatched positions        
  if noMatch == 1:
      numUnmatchedPositions +=1
      strandList.append('u') #adds strand of unknown if not found
      if numEntries <= 2:
          numUniqueUnmatchedPositions +=1
          
      #prepare output based on setting for unmatched only
      if options.output == 'text':
          print bpm.names[i] ,bpm.chr[i], bpm.pos[i], 'blat', blatEntries
      elif options.output == 'fasta':
          print '>'+bpm.names[i]+'\n'+bpm.sourceStrand[i]
      elif options.output == 'doublefasta':
          sourceSeq = re.split('\[|\/|\]',bpm.sourceStrand[i]) #splits to 4 element list
          fivePrimeSeq = sourceSeq[0]

          if sourceSeq[1] == '-':
              alleleASeq = ''
          else:
              alleleASeq = sourceSeq[1]

          if sourceSeq[2] == '-':
              alleleBSeq = ''
          else:
              alleleBSeq = sourceSeq[2]

          threePrimeSeq = sourceSeq[3]

          alleleAStrand = fivePrimeSeq + alleleASeq + threePrimeSeq
          alleleBStrand = fivePrimeSeq + alleleBSeq + threePrimeSeq

          print '>'+bpm.names[i]+'\n'+alleleAStrand
          print '>'+bpm.names[i]+'\n'+alleleBStrand

  if bpm.chr[i] == "0":
      numZeroBpm += 1

if options.output == 'strandlist': #prepare strand list output for use in gtc2ped.py
    for i in xrange(numSNPs):
          print bpm.names[i],strandList[i] #for each i prints a name and strand

if options.output == 'diagnostics':  
  print 'Forward   : ', numForward
  print 'Reverse   : ', numReverse
  print 'Unassigned: ', numUnmatchedPositions
  print 'Matched   : ', numMatchedPositions
  print 'Unique Aligned Unmatched:', numUniqueUnmatchedPositions
  print 'Chr in BPM is 0', numZeroBpm
