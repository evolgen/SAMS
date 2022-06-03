#!/usr/bin/env python3

import sys

import string
import re
import bx.align.maf
import pandas as pd
import numpy as np
import re
import collections
import subprocess

from optparse import OptionParser

import sized_slidingwindow
import species_info_fromtree
import data_extract

parser = OptionParser(usage="usage: calculate_nucl_transitions.py -m file.maf -l names.list -n 5")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")
parser.add_option("-t", "--inputtree", action="store", dest='inputtree', type="string")
parser.add_option("-l", "--specieslist", action="store", dest='specieslist', type="string")
parser.add_option("-o", "--outtsv", action="store", dest='outtsv', type="string", default=sys.stdout)

(options, args) = parser.parse_args()

#
accepted_nucleotides = "ATGCatgc"


#big_maffile1 = open('/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/cactus2/test/test.maf', 'r')
big_maffile1 = open(options.inputmaf, 'r')
maf_reader1 = bx.align.maf.Reader(big_maffile1)

#custom_species_list1 = open('/global/scratch2/rohitkolora/Rockfish/Genomes/traits/species_with_ages.list', 'r')
custom_species_listfile = open(options.specieslist, 'r') # Selected species list
custom_species_list1 = custom_species_listfile.read().splitlines()
custom_species_listfile.close()

# Extract ancestral info for species from tree
#main_tree1_path = "/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/cactus2/test/tree_seb77.nwk"
main_tree1_path = str(options.inputtree) # Selected species list
species_to_ancestor_sortdist_df1 = species_info_fromtree.species_to_ancestor_sortdist_df(main_tree1_path)

# Create df of species and ancestors
anc_species_anc_dict1 = species_info_fromtree.conv_sp2ancdist_df2dict(species_to_ancestor_sortdist_df1)

### Fun starts here
# Create dict for counting transitions
list_triplet_combinations = sized_slidingwindow.string_to_sized_combinations('ATGC',3)
species_count_dict1 = {}
for pattern_elem2 in list_triplet_combinations:
    species_count_dict1[pattern_elem2] = {}
    for pattern_elem3 in list_triplet_combinations:
        if(pattern_elem2[0]==pattern_elem3[0] and pattern_elem2[2]==pattern_elem3[2]):
            species_count_dict1[pattern_elem2][pattern_elem3] = 0

initial_species_counter1 = {}
for key_elem1 in custom_species_list1: # anc_species_anc_dict1.keys(): # Creating nested dictionary of string patterns for each species
    initial_species_counter1[key_elem1] = dict(species_count_dict1)

# Play with MAF blocks
block_counter = 0
for block_maf1 in maf_reader1: # We go block by block here - scope for parallelisation later
    block_counter = block_counter + 1
    block_now1 = pd.DataFrame(block_maf1.components) # creating df for the maf component
    
    block_now2 = block_now1[0].apply(lambda x: pd.Series(str(x).split(" ")))[[1,6]] # extracting seq header and sequence
    block_now2.columns = ["header", "sequence"]
    block_now2['species'] = block_now2['header'].apply(lambda x: pd.Series(str(x).split(".")))[[0]] # creating species names as column
    
    sp_all1 = collections.Counter(block_now2['species'])    
    all_species_block_now2 = [i1 for i1,j1 in sp_all1.items() if j1==1] # Non redundant species alone
    
    anc_species_block_now2 = [speciesname for speciesname in all_species_block_now2 if re.search('Anc[0-9]*$', speciesname)] # filtering for non ancestral alone
    
    filt_species_block_now = [speciesname for speciesname in all_species_block_now2 if not re.search('Anc[0-9]*$', speciesname)] # filtering for non ancestral alone
    filt_species_block_now2 = list(set(custom_species_list1) & set(filt_species_block_now)) # retaining only selected species
    
    for curr_species1 in filt_species_block_now2: # Extracting data from tree
        sequence_curr_species1 = data_extract.extract_sequence_4df(block_now2, curr_species1) # Extracting sequence for species
        
        blocks_sequence_curr_species1 = sized_slidingwindow.create_list_slidingwindow(sequence_curr_species1,3)
        
        curr_species_anc_list1 = anc_species_anc_dict1[curr_species1] # extract ancestral list for species
        for anc_species1 in curr_species_anc_list1: # Iterating through anc. list for those in maf block
            if anc_species1 not in anc_species_block_now2 or anc_species1 not in all_species_block_now2:
                next
            sequence_curr_ancestral1 = data_extract.extract_sequence_4df(block_now2, anc_species1) # Extracting sequence for species
            break
            
        blocks_sequence_curr_ancestral1 = sized_slidingwindow.create_list_slidingwindow(sequence_curr_ancestral1,3)
        
        for count1 in range(0, len(blocks_sequence_curr_species1), 1):
            curr_anc_nucl1 = blocks_sequence_curr_ancestral1[count1] # Storing ancestral pattern
            curr_species_nucl1 = blocks_sequence_curr_species1[count1] # Storing species pattern
            
#            if(len(blocks_sequence_curr_species1) != len(blocks_sequence_curr_ancestral1)): 
#                next
                
            if all(nucl1 in accepted_nucleotides for nucl1 in curr_species_nucl1): # Checking for ATGCs alone
                if all(nucl1 in accepted_nucleotides for nucl1 in curr_anc_nucl1):
                    #initial_species_counter1[curr_species1][curr_anc_nucl1][curr_species_nucl1] = min(initial_species_counter1[curr_species1][curr_anc_nucl1].get(curr_species_nucl1,0)) + 1 #Counting instances of nucleotides
                    if(curr_anc_nucl1[0] == curr_mrca1_nucl1[0] and
                            curr_anc_nucl1[2] == curr_anc_nucl1[2]):
                        initial_species_counter1[curr_species1][curr_anc_nucl1][curr_species_nucl1] = initial_species_counter1[curr_species1][curr_anc_nucl1][curr_species_nucl1] + 1 #Counting instances of nucleotides
                    
    if(block_counter%100 == 0):
        print(block_counter, end=" ", flush=True)

print("\n\nProcessed all blocks - ", block_counter, "in total\n")

final_species_counter1 = data_extract.counter_dict_conv_values(initial_species_counter1)    
final_species_counter_df = data_extract.counter_dict_2_df(final_species_counter1)

big_maffile1.close()

#final_species_counter_df.to_csv('/global/scratch2/rohitkolora/Rockfish/Genomes/alignments/cactus2/test/outtest.tsv', sep='\t', index=False) # Tab separated file
final_species_counter_df.to_csv(options.outtsv, sep='\t', index=False)

print("Saved to output file too!!!")

