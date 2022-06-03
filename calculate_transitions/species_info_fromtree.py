#!/usr/bin/env python3

import sys
import itertools
from Bio import Phylo
from tabulate import tabulate
import pandas as pd
import re

import sized_slidingwindow

def generate_pairs(self):
    pairs = itertools.tee(self)
    pairs[1].__next__()
    return zip(pairs[0], pairs[1])

def terminal_neighbor_dists(self):
    """Return a list of distances between adjacent terminals."""
    print("#species1\tspecies2\tdistance")
    for iname1 in generate_pairs(main_tree1.find_clades()):
        if re.search('^Anc[0-9]+', iname1[0].name):
            continue
        print(iname1[0].name, iname1[1].name, main_tree1.distance(*iname1), sep="\t")

def species_to_ancestor_sortdist_df(curr_tree_nwk_path):
    main_tree1 = Phylo.read(curr_tree_nwk_path, "newick")
    
    terms = [term for term in main_tree1.get_terminals()]
    list_terms1 = []
    for curr1 in range(0,len(terms)):
        list_terms1.append(terms[curr1].name)
    
    nonterms = [nonterm for nonterm in main_tree1.get_nonterminals()]
    list_nonterms1 = []
    for curr2 in range(0,len(nonterms)):
        list_nonterms1.append(nonterms[curr2].name)
        
    dist_term_nonterm_all = []    
    for curr_term1 in list_terms1:        
        for curr_nonterm1 in list_nonterms1:
            arr_curr_pair1 = [curr_term1, curr_nonterm1, main_tree1.distance(curr_term1,curr_nonterm1)]
            dist_term_nonterm_all.append(arr_curr_pair1)
    
    dist_term_nonterm_all = sorted(dist_term_nonterm_all, key=lambda x: (x[0], x[2]), reverse=False)
    
    return dist_term_nonterm_all
    
def conv_sp2ancdist_df2dict(df_dist_term_nonterm_df_pd):    
    dist_term_nonterm_all_df = pd.DataFrame(df_dist_term_nonterm_df_pd, columns=['species1','species2','distance'], dtype=float)
    
    #anc_species_anc_dict1 = dist_term_nonterm_all_df.groupby(['species1'])['species2'].agg({'species1': ','.join})['species1'].to_dict() # Creates problem with python versions
    anc_species_anc_dict1 = dist_term_nonterm_all_df.groupby(['species1'])['species2'].apply(','.join).to_dict()


    for anc_key1, anc_value1 in anc_species_anc_dict1.items():
        anc_species_anc_dict1[anc_key1] = anc_value1.split(',')
    
    return anc_species_anc_dict1

def initial_species_counter(anc_species_anc_dict1):
    species_count_dict1 = {}
    list_triplet_combinations = sized_slidingwindow.string_to_sized_combinations('ATGC',3)
    
    for pattern_elem2 in list_triplet_combinations: # Creating string pattern dictionary
        species_count_dict1[pattern_elem2] = [0,0]

    species_count_dict2 = {}
    for key_elem1 in anc_species_anc_dict1.keys(): # Creating nested dictionary of string patterns for each species
        species_count_dict2[key_elem1] = species_count_dict1
    
    return species_count_dict2

def conv_sp2ancdist_df2dict_local(df_dist_term_nonterm_df_pd):    
    dist_term_nonterm_all_df = pd.DataFrame(df_dist_term_nonterm_df_pd, columns=['species1','species2','distance'], dtype=float)
    
    anc_species_anc_dict1_local = dist_term_nonterm_all_df.groupby(['species1']).agg(species2 = ('species2', lambda x: ','.join(x))).reset_index()
    anc_species_anc_dict1 = dict(zip(anc_species_anc_dict1_local.species1, anc_species_anc_dict1_local.species2))
    ###anc_species_anc_dict1 = pd.Series(anc_species_anc_dict1_local.species2.values, index=anc_species_anc_dict1_local.species1).to_dict()
    
    for anc_key1, anc_value1 in anc_species_anc_dict1.items():
        anc_species_anc_dict1[anc_key1] = anc_value1.split(',')
    
    return anc_species_anc_dict1
def get_ancestor_4speciespairs(main_tree1):

    species_pair_anc_df = []
    for present_name1 in main_tree1.get_terminals():
        if re.search('^Anc[0-9]+', present_name1.name) or re.search('^OROOT', present_name1.name):
            continue
        for present_name2 in main_tree1.get_terminals():
            if re.search('^Anc[0-9]+', present_name2.name) or re.search('^OROOT', present_name2.name):
                continue
            if present_name1==present_name2:
                continue

            current_bl1 = main_tree1.common_ancestor(present_name1.name, present_name2.name)
            #print(present_name1.name, present_name2.name, current_bl1.name, sep="\t")
            species_pair_anc_df.append([present_name1.name, present_name2.name, current_bl1.name])
    output_df1 = pd.DataFrame(species_pair_anc_df, columns=['species1', 'species2', 'mrca'])
    return(output_df1)
