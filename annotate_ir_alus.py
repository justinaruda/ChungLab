# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 23:02:35 2021

@author: justi
"""

import pandas as pd
import HTSeq
import numpy as np
import getopt, sys
import re

def remove_empty_sets(tuples):
    new_tuples = []
    for tup in tuples:
        if len(tup[1])!=0:
             new_tuples.append(tup)
    return new_tuples

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]

short_options = "t:a:o:"
long_options = ["table=", "annotation_file=", "output="]

try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    # Output error, and return with an error code
    print (str(err))
    sys.exit(2)

table_file = None
annotation_file = None
out_file = None

for current_argument, current_value in arguments:
    if current_argument in ("-t", "--table"):
        table_file = current_value
    if current_argument in ("-a", "--annotation_file"):
        annotation_file = current_value
    if current_argument in ("-o", "--output"):
        out_file = current_value

print('Reading input files...')

clip_df = pd.read_table(table_file)
repeat_annots = pd.read_table(annotation_file)

print('Loading repeat annotation file...')
family_gas = HTSeq.GenomicArrayOfSets("auto",stranded=False)
for row in repeat_annots.iterrows():
    row = row[1]
    # if row['repFamily'] != 'Alu':
    #     continue
    iv = HTSeq.GenomicInterval(row['#genoName'], row['genoStart'], row['genoEnd'], strand='.')
    family_gas[iv] += (row['repName'], row['repFamily'],iv,row['strand'])

print('Repeat annotation file loaded...')
print('Starting annotation...')

iv_regex = re.compile(r'>(chr\d*):(\d*)-(\d*)')

clip_df['Repeat'] = 'None'
clip_df['Repeat_Context'] = 'None'
for idx,row in clip_df.iterrows():
    iv_match = iv_regex.search(row['Interval'])
    if iv_match is not None:
        region = iv_match.group(1)
        start = int(iv_match.group(2))
        end = int(iv_match.group(3))
    else:
        continue
    clip_iv = HTSeq.GenomicInterval(region, start, end)
    families  = []
    sets = [entry for entry in family_gas[clip_iv].steps()]
    sets = remove_empty_sets(sets)
    repeats_tuples = [list(val)[0] for iv,val in sets]
    families.extend([family for name,family,iv,strand in repeats_tuples])
    families_set = set(families)
    if len(families_set) != 0:
        repeat_families_str = ' '.join(sorted(families_set))
        clip_df.at[idx,'Repeat'] = repeat_families_str
        if 'Alu' in repeat_families_str:
            repeat_context = 'individual'
            ir_counter = 0
            alu_counter = 0
            last_strand = None
            current_strand = None
            for name,family,iv,strand in repeats_tuples:
                if family == 'Alu':
                    if current_strand == None:
                        last_strand = strand
                        current_strand = strand
                    else:
                        current_strand = strand
                    if strand != '.' and current_strand != last_strand:
                        ir_counter += 1
                        last_strand = current_strand
                    alu_counter += 1
                if alu_counter > 1:
                    repeat_context = 'cluster'
                if ir_counter > 0:
                    repeat_context = 'ir_cluster'
            clip_df.at[idx,'Repeat_Context'] = repeat_context

f = open(out_file, 'w')
clip_df.to_csv(f, sep="\t", line_terminator="\n", index=False, header=True)
f.close()
print('Done.')