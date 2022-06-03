#!/usr/bin/env python3.6

import sys
import Bio
from Bio import Phylo

from optparse import OptionParser
import itertools
import string
import re
import bx.align.maf
import pandas as pd
import numpy as np
import re
import collections

import sized_slidingwindow
import species_info_fromtree
import data_extract

parser = OptionParser(usage="usage: calculate_nucl_transitions_speciespairs.py -m file.maf -l names.list -t tree.nwk")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")
parser.add_option("-t", "--inputtree", action="store", dest='inputtree', type="string")
parser.add_option("-l", "--specieslist", action="store", dest='specieslist', type="string")
parser.add_option("-o", "--outtsv", action="store", dest='outtsv', type="string", default=sys.stdout)

(options, args) = parser.parse_args()

#
accepted_nucleotides = "ATGCatgc"

list_triplet_combinations = sized_slidingwindow.string_to_sized_combinations('ATGC',3)

big_maffile1 = open(options.inputmaf, 'r')
maf_reader1 = bx.align.maf.Reader(big_maffile1)

main_tree1_path = str(options.inputtree) # Selected species list
species_pairs_ancestors_list1 = species_info_fromtree.get_ancestor_4speciespairs(main_tree1 = Phylo.read(main_tree1_path, "newick"))

custom_species_list1 = open(options.specieslist, 'r') # Selected species list
custom_species_list1 = custom_species_list1.read().splitlines()

print("Read input maf, tree, list")

species_pairs_combination_list1 = list(itertools.combinations(custom_species_list1, 2))
species_pairs_combination_df1 = pd.DataFrame(species_pairs_combination_list1, columns =['species1', 'species2']) 

species_pairs_ancestors_df1 = pd.merge(species_pairs_combination_df1, species_pairs_ancestors_list1,
         how='inner', on=['species1', 'species2'])

#species_pairs_ancestors_df1_nh = species_pairs_ancestors_df1.to_string(header=False)
species_pairs_ancestors_df1_list1 = species_pairs_ancestors_df1.values.tolist()

triplets_thrice = list(itertools.product(list_triplet_combinations,list_triplet_combinations,list_triplet_combinations))
triplets_thrice_list1 = []
for curr_elem_lst1 in range(0, len(triplets_thrice)):
    now_value1 = str('__'.join(map(str, triplets_thrice[curr_elem_lst1])))
    triplets_thrice_list1.append(now_value1)
    
for key_elem1 in species_pairs_ancestors_df1_list1:
    key_elem1 = str('__'.join(map(str, key_elem1)))
    for pattern_elem_thrice in triplets_thrice_list1: # Creating string pattern dictionary
        pattern_elem_thrice = str(pattern_elem_thrice)
        value_init1 = 0
        exec(f'{key_elem1}__{pattern_elem_thrice} = value_init1')
print("Species pairs and ancestor list created")
                
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

    df_sp1 = pd.DataFrame(filt_species_block_now2, columns =['species1'])
    df_sp2 = pd.DataFrame(filt_species_block_now2, columns =['species2'])
    df_mrca = pd.DataFrame(anc_species_block_now2, columns =['mrca'])
    
    df_merge_sp_1 = pd.merge(species_pairs_ancestors_df1, df_sp1, how='inner', on=['species1'])
    df_merge_sp_1_2 = pd.merge(df_merge_sp_1, df_sp2, how='inner', on=['species2'])
    df_merge_sp_1_2_mrca = pd.merge(df_merge_sp_1_2, df_mrca, how='inner', on=['mrca'])
    
    for curr_row1 in df_merge_sp_1_2_mrca.index:
        curr_species1 = str(df_merge_sp_1_2_mrca['species1'][curr_row1])
        curr_species2 = str(df_merge_sp_1_2_mrca['species2'][curr_row1])
        curr_mrca1 = str(df_merge_sp_1_2_mrca['mrca'][curr_row1])

#        print(curr_species1, curr_species2, curr_mrca1)
        sequence_curr_species1 = data_extract.extract_sequence_4df(block_now2, curr_species1) # Extracting sequence for species1
        blocks_sequence_curr_species1 = sized_slidingwindow.create_list_slidingwindow(sequence_curr_species1,3)

        sequence_curr_species2 = data_extract.extract_sequence_4df(block_now2, curr_species2) # Extracting sequence for species2
        blocks_sequence_curr_species2 = sized_slidingwindow.create_list_slidingwindow(sequence_curr_species2,3)

        sequence_curr_mrca1 = data_extract.extract_sequence_4df(block_now2, curr_mrca1) # Extracting sequence for mrca
        blocks_sequence_curr_mrca1 = sized_slidingwindow.create_list_slidingwindow(sequence_curr_mrca1,3)
        
        for count1 in range(0, len(blocks_sequence_curr_species1), 1):
            curr_species1_nucl1 = blocks_sequence_curr_species1[count1] # Storing species pattern
            curr_species2_nucl1 = blocks_sequence_curr_species2[count1] # Storing species pattern
            curr_mrca1_nucl1 = blocks_sequence_curr_mrca1[count1] # Storing species pattern

            if all(nucl1 in accepted_nucleotides for nucl1 in curr_species1_nucl1): # Checking for ATGCs alone
                if all(nucl1 in accepted_nucleotides for nucl1 in curr_species2_nucl1):
#                    if(pattern_elem2 == pattern_elem3):
#                        next
                    if all(nucl1 in accepted_nucleotides for nucl1 in curr_mrca1_nucl1):
                        #print(curr_species1, curr_species2, curr_mrca1, curr_species1_nucl1, curr_species2_nucl1, value_now1)
                        value_now1 = int(globals()[f'{curr_species1}__{curr_species2}__{curr_mrca1}__{curr_species1_nucl1}__{curr_species2_nucl1}__{curr_mrca1_nucl1}']) + 1
                        exec(f'{curr_species1}__{curr_species2}__{curr_mrca1}__{curr_species1_nucl1}__{curr_species2_nucl1}__{curr_mrca1_nucl1} = value_now1')
    
    if(block_counter%100 == 0):
        print(block_counter, end=",", flush=True)

    
print("Processed blocks :", block_counter)
initial_species_counter1 = {}
species_count_dict1 = {}
for key_elem1 in species_pairs_ancestors_df1_list1:
    key_elem1 = str('__'.join(map(str, key_elem1)))
    initial_species_counter1[key_elem1] = {}
    for pattern_elem_thrice in triplets_thrice_list1: # Creating string pattern dictionary
        pattern_elem_thrice = str(pattern_elem_thrice)
        curr_value1 = int(globals()[f'{key_elem1}__{pattern_elem_thrice}'])
        initial_species_counter1[key_elem1][pattern_elem_thrice] = curr_value1

final_species_counter1 = data_extract.pairspecies_counter_dict_conv_values(initial_species_counter1)
final_species_counter_df = data_extract.pairwise_counter_dict_2_df(final_species_counter1)

final_species_counter_df.to_csv(options.outtsv, sep='\t', index=False) # Tab separated file

print("Saved to output file too!!!")


