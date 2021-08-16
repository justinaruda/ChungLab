# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 12:26:35 2021

@author: justi
"""

import HTSeq
import pandas as pd
import getopt, sys

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]

short_options = "b:o:"
long_options = ["bed=", "output="]

try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    # Output error, and return with an error code
    print (str(err))
    sys.exit(2)

bed_file = None
bed_output = None
out_bed = False

for current_argument, current_value in arguments:
    if current_argument in ("-b", "--bed"):
        bed_file = current_value
    if current_argument in ("-o", "--output"):
        bed_output = current_value
        out_bed = True
    
difference_chrom = []
difference_start = []
difference_end = []
difference_name = []
difference_score = []
difference_strand = []
    
def log_tag(tag):
    global difference_chrom
    global difference_start
    global difference_end
    global difference_score
    global diference_strand
    difference_chrom.append(tag.iv.chrom)
    difference_start.append(tag.iv.start)
    difference_end.append(tag.iv.end)
    difference_name.append(tag.name)
    difference_score.append(tag.score)
    difference_strand.append(tag.iv.strand)

    
bed_reader = HTSeq.BED_Reader(bed_file)


k = 0
used = set()

for tag in bed_reader:
    if tag.name in used:
        k = k+1
        print("DUPLICATE FOUND!")
        print(tag)
        continue
    else:
        log_tag(tag)

print("Removed " + str(k) + " duplicates.")
bed = pd.DataFrame()
bed['chrom'] = difference_chrom
bed['start'] = difference_start
bed['end'] = difference_end
bed['name'] = difference_name
bed['score'] = difference_score
bed['strand'] = difference_strand

f = open(bed_output, 'w')
bed.to_csv(f, sep="\t", line_terminator="\n", index=False, header=False)
f.close()

# python remove_duplicates.py -b A1KO_Mock_Down.quantified_tags.bed -o A1KO_Mock_Down.quantified_tags.dedup.bed