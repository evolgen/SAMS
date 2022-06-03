#!/usr/bin/env python3

# This is to create permutations of strings, here we are attempting permutate trinucleotides

import string
import sys
from optparse import OptionParser
import itertools

#string_dna = "ATGC"

required="inputstring".split()

def __main__():

    parser = OptionParser(usage="usage: string_permutations.py -i ABCD -s 3")
    parser.add_option("-i", "--inputstring", action="store", dest='inputstring', type="string")
    parser.add_option("-s", "--permutesize", action="store", type="int", default=2)
#    parser.add_option("-o", "--outputfile", action="store", default="")

    (options, args) = parser.parse_args()
    
    for opt1a in required:
        if options.__dict__[opt1a] is None:
            parser.error("Parameter required : %s"%opt1a)
    
    string_dna = options.inputstring
    size_permutation = options.permutesize
    
    if(size_permutation > len(string_dna)):
        print("The window size is",size_permutation,"\n Oopsies this cannot be larger than input sequence length :",len(string_dna))
        exit()

    trinuc_dnastring = list(itertools.permutations(string_dna, size_permutation))

    for elem1 in trinuc_dnastring:
        print(''.join(elem1))


if __name__ == "__main__":
    __main__()    