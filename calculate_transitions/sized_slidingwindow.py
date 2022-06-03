#!/usr/bin/env python3

# Creating sliding windows of a sequence string

import sys
import itertools
import string
from optparse import OptionParser
import itertools

def string_to_sized_combinations(string_dna, size_permutation):

    if(size_permutation > len(string_dna)):
        print("The window size is",size_permutation,"\n This cannot be larger than input-seq length :",len(string_dna))
        exit()
    if type(string_dna) != str:
        print("Please use string in the function such as : string_to_sized_combinations('XYZ',2) \n")
        exit()
    if type(size_permutation) != int:
        print("Please use integer in the function such as : string_to_sized_combinations('XYZ',2) \n")
        exit()
    
    #trinuc_dnastring = list(itertools.permutations(string_dna_new, size_permutation))
    trinuc_dnastring = list([p1 for p1 in itertools.product(string_dna, repeat=size_permutation)])


    permuate_strings1 = []
    for elem1 in trinuc_dnastring:
        permuate_strings1.append(''.join(elem1))
        
    return permuate_strings1


if __name__ == "__main__":
    __main__()    

def prepend(list, str):
    str += '{0}'
    list = [str.format(i) for i in list]
    return(list)

def string_slidingwindow(stringseq, window_size=3):
    for elem2 in range(len(stringseq) - window_size + 1):
        yield stringseq[elem2:elem2+window_size]

def create_dict_slidingwindow(seqstring1, window_size=3):      
    triplet_blocks = []
    seqstring1 = seqstring1.strip('\n')
    for seq_curr1 in string_slidingwindow(seqstring1, window_size):
        triplet_blocks.append(seq_curr1)
    
    index_list1 = range(1,len(triplet_blocks)+1,1)
    index_list1 = [*index_list1]
    prepend_str2 = "B"
    index_list1 = prepend(index_list1, prepend_str2)
    index_list1 = [*index_list1]

    sequence_dict1 = dict(zip(index_list1, triplet_blocks))
    return sequence_dict1

def create_list_slidingwindow(seqstring1, window_size=3):      
    triplet_blocks = []
    seqstring1 = seqstring1.replace('\n','')
    seqstring1 = seqstring1.replace(' ','')
    for seq_curr1 in string_slidingwindow(seqstring1, window_size):
        triplet_blocks.append(seq_curr1)

    return triplet_blocks    

#nucl_seq1 = "ATACACAGCAGCGAGGAGGAGT"    
#create_dict_slidingwindow(nucl_seq1, 4)

if __name__ == "__main__":
    string_to_sized_combinations()
    prepend()
    string_slidingwindow()
    create_dict_slidingwindow()
    create_list_slidingwindow()
    