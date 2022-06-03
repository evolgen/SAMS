#!/usr/bin/env python3

import sys
import string
import re
import bx.align.maf
import pandas as pd
import collections
from optparse import OptionParser
import gzip
import pathlib

parser = OptionParser(usage="usage: filter_maf_singlecopy.py -m large.maf")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")

(options, args) = parser.parse_args()

#if '.gz' in pathlib.Path(options.inputmaf).suffixes:
#    big_maffile1 = gzip.open(options.inputmaf, 'r')
#else:
big_maffile1 = open(options.inputmaf, 'r')

maf_reader1 = bx.align.maf.Reader(big_maffile1)

print("##maf version=1\n")

block_counter = 0
for block_maf1 in maf_reader1: # We go block by block here - scope for parallelisation later
    block_counter = block_counter + 1
    block_now1 = pd.DataFrame(block_maf1.components) # creating df for the maf component
    
    block_now2 = block_now1[0].apply(lambda x: pd.Series(str(x).split(" ")))[[1,6]] # extracting seq header and sequence
    block_now2.columns = ["header", "sequence"]
    block_now2['species'] = block_now2['header'].apply(lambda x: pd.Series(str(x).split(".")))[[0]] # creating species names as column
    
    sp_all1 = collections.Counter(block_now2['species'])    
    single_species1 = [i1 for i1,j1 in sp_all1.items() if j1==1] # Getting single copy species
    
    start_single_species1 = ["s "+str1 for str1 in single_species1] # Modifying start for matching later

    if len(single_species1) <= 1: # Filtering blocks with <2 single copy species in blocks 
        next
    
    print("a score=0")
    for lines1 in block_maf1.components:
        test_truer = list(filter(str(lines1).startswith, start_single_species1)) != []
        if test_truer:
            print(str(lines1))#.replace(" ","\t"))

maf_reader1.close()            
        
