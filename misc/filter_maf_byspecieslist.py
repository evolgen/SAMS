#!/usr/bin/env python3

# Filter maf blocks by species lists and minimum species requirement

import sys

import string
import re
from optparse import OptionParser

import bx.align.maf

parser = OptionParser(usage="usage: filter_maf_byspecieslist.py -m file.maf -l names.list -n 5")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")
parser.add_option("-l", "--specieslist", action="store", dest='specieslist', type="string")
parser.add_option("-n", "--minspecies", action="store", dest='minspecies', type="int", default=10)


(options, args) = parser.parse_args()
    
big_maffile1 = open(options.inputmaf, 'r')
maf_reader1 = bx.align.maf.Reader(big_maffile1)

custom_species_list1 = open(options.specieslist, 'r')
custom_species_list1 = custom_species_list1.read().splitlines()

minspecies = options.minspecies

# Let's print the header
print("##maf version=1\n")

block_number = 0
for block_maf1 in maf_reader1:
    block_number = block_number + 1
#    print("Processing Block - ", block_number)
    
    numspecies1 = len(block_maf1.components)

    custom_species_names1 = []
    
    for linenow1 in block_maf1.components:
        speciesname1, sequencename1 = linenow1.src.split(".",1)
        if(speciesname1 in custom_species_list1):    
            custom_species_names1.append(speciesname1)
#        elif(re.match('^Anc[0-9]*$', speciesname1)) :
#            print(speciesname1)
    
    custom_species_total1 = sorted(set(custom_species_names1))
    if(len(custom_species_names1) >= minspecies):
        print(block_maf1)#.replace(" ","\t"))


big_maffile1.close()


