#!/usr/bin/env python3

# Maf block stats by species lists

import sys

import string
import re
from optparse import OptionParser

import bx.align.maf

parser = OptionParser(usage="usage: filter_maf_byspecieslist.py -i ABCD -s 3")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")
parser.add_option("-l", "--specieslist", action="store", dest='specieslist', type="string")


(options, args) = parser.parse_args()
    
big_maffile1 = open(options.inputmaf, 'r')
maf_reader1 = bx.align.maf.Reader(big_maffile1)

custom_species_list1 = open(options.specieslist, 'r')
custom_species_list1 = custom_species_list1.read().splitlines()

# Let's print the header
print("#blocknum\tsize\ttotalseq\tnumspecies")

block_number = 0
for block_maf1 in maf_reader1:
    block_number = block_number + 1
#    print("Processing Block - ", block_number)
    
    numseqs1 = len(block_maf1.components)

    custom_species_names1 = []
    
    curr_seqnum1 = 0
    for linenow1 in block_maf1.components:
        curr_seqnum1+=1
        
        if(curr_seqnum1==1):
            curr_blocksize1 = linenow1.size
            
        speciesname1, sequencename1 = linenow1.src.split(".",1)
        if(speciesname1 in custom_species_list1):    
            custom_species_names1.append(speciesname1)
    
    custom_species_total1 = sorted(set(custom_species_names1))
    
    print("block_#",block_number,"\t",curr_blocksize1,"\t",numseqs1,"\t",len(custom_species_total1), sep='')
    

big_maffile1.close()


