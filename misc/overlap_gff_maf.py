#!/usr/bin/env python3

# Maf block to bed regions by species name

import sys

import string
import re
from optparse import OptionParser
import io

import bx.align.maf
import pandas
import numpy

parser = OptionParser(usage="usage: extract_bed_byspecies.py -m large.maf -g large.gff")
parser.add_option("-m", "--inputmaf", action="store", dest='inputmaf', type="string")
parser.add_option("-g", "--inputgff", action="store", dest='inputgff', type="string")

(options, args) = parser.parse_args()

big_gfffile1 = open(options.inputgff, 'r')
sys.stderr.write("Processing the GFF file \n")

input_gfffile1 = pandas.read_csv(big_gfffile1, delimiter="\t", sep='\t', chunksize=500000, iterator=True, skip_blank_lines=True, header=None)#.dropna()
input_gfffile2 = pandas.concat(input_gfffile1)
input_gfffile1 = []
big_gfffile1.close()

input_bedfile1 = pandas.DataFrame(input_gfffile2.iloc[:, [0,3,4]])
input_bedfile1.columns = ['seqname','seqstart_bed','seqend_bed']
sys.stderr.write("BED file created from GFF\n")


big_maffile1 = open(options.inputmaf, 'r')
maf_reader1 = bx.align.maf.Reader(big_maffile1)

print("##maf version=1\n")

for block_maf1 in maf_reader1:

    current_maf_bed = pandas.DataFrame([], columns = ["seqname", "seqstart_maf","seqend_maf"])
    for linenow1 in block_maf1.components:
#        speciesname1, sequencename1 = linenow1.src.split(".",1)
        
        seqname1 = linenow1.src
        seqstart1 = linenow1.forward_strand_start
        seqend1 = linenow1.forward_strand_end            

        line_now_bed1 = [seqname1, seqstart1, seqend1]
        linenow1_series = pandas.Series(line_now_bed1, current_maf_bed.columns)
        current_maf_bed = current_maf_bed.append(linenow1_series, ignore_index=True)
    merged_maf_coords = current_maf_bed.merge(input_bedfile1, on='seqname', how='inner')#.drop_duplicates(keep=False,inplace=True)

    filt_merged_maf_coords = merged_maf_coords.query("(seqstart_maf >= seqstart_bed & seqend_maf <= seqend_bed) | (seqstart_maf <= seqstart_bed & seqend_maf >= seqend_bed) | (seqstart_maf <= seqstart_bed & seqend_maf <= seqend_bed & seqend_maf >= seqstart_bed) | (seqstart_maf >= seqstart_bed & seqend_maf >= seqend_bed & seqstart_maf <= seqend_bed)") 
    if(len(filt_merged_maf_coords.index) > 0):
        print(block_maf1)
        ###outputnowstd=io.StringIO(); filt_merged_maf_coords.to_csv(outputnowstd, sep="\t", encoding='utf-8', index=False); outputnowstd.seek(0); print(outputnowstd.read());

big_maffile1.close()

