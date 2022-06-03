#!/usr/bin/env python3

# Maf block to bed regions by species name

import sys

import string
import re
from optparse import OptionParser

import bx.align.maf

parser = OptionParser(usage="usage: extract_bed_byspecies.py -m large.maf -s Sebastes_aleutianus")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")
parser.add_option("-s", "--speciesname", action="store", dest='speciesname', type="string")


(options, args) = parser.parse_args()
    
big_maffile1 = open(options.inputmaf, 'r')
maf_reader1 = bx.align.maf.Reader(big_maffile1)

filterspeciesname1 = options.speciesname

for block_maf1 in maf_reader1:

    for linenow1 in block_maf1.components:
        speciesname1, sequencename1 = linenow1.src.split(".",1)
        
        if(speciesname1 == filterspeciesname1):    
            seqname1 = linenow1.src
            seqstart1 = linenow1.forward_strand_start
            seqend1 = linenow1.forward_strand_end            

            print(seqname1,"\t",seqstart1,"\t",seqend1,"\t.\t.\t+", sep='')
            
big_maffile1.close()



