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

short_options = "a:b:o:kx"
long_options = ["output=","keep_duplicates", "retain_ratio"]

try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    # Output error, and return with an error code
    print (str(err))
    sys.exit(2)

my_file = None
bed_file = None
bed_output = None
keep_duplicates = False
retain_ratio = False

for current_argument, current_value in arguments:
    if current_argument in ("-a"):
        my_file = current_value
    if current_argument in ("-b"):
        bed_file = current_value
    if current_argument in ("-o", "--output"):
        bed_output = current_value
    if current_argument in ("-k", "--keep_duplicates"):
        keep_duplicates = True
    if current_argument in ("-x", "--retain_ratio"):
        retain_ratio = True

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

      
my_df = pd.read_table(my_file, header=None, names=["chrom","start","end","name","score","strand"])

bed_reader = HTSeq.BED_Reader(bed_file)

name_count = {}

if retain_ratio:
    all_names = my_df["name"]
    for name in all_names:
        if name in name_count.keys():
            name_count[name] += 1
        else:
            name_count[name] = 1
uniq_names = set(my_df["name"].unique())    
i = 0
k = 0
used = set()
for tag in bed_reader:
    name = tag.name
    if name in uniq_names:
        i = i+1
        if retain_ratio:
            name_count[name] -= 1
            print(name_count[name])
            if name_count[name] < 0:
                log_tag(tag)
            continue
        if name in used:
            k = k+1
            print("DUPLICATE FOUND!")
            print(tag)
            if keep_duplicates:
                log_tag(tag)
        used.add(name)
    else:
        log_tag(tag)
action = "Removed "
if keep_duplicates:
    action = "Retained "

print("Processed " + str(i) + " entries.")
print(action + str(k) + " duplicates.")
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

# python find_differences.py -a A1KO_Mock_Down.quantified_tags.bed -b A1KO_Mock_Down.repeats.quantified.bed -o difference.bed -k