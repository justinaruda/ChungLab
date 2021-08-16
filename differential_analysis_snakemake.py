# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 16:36:00 2021

@author: justi
"""

import pandas as pd
from collections import OrderedDict
import re
import sys
import os
from pathlib import Path

def compare(first,second):
    global master_dict
    first_name,first_ext = first
    second_name,second_ext = second
    print("Comparing {0} and {1}...".format(first_name,second_name))
    first_dict = master_dict[first]
    first_total = sum(list(first_dict.values()))
    second_dict = master_dict[second]
    second_total = sum(list(second_dict.values()))
    difference_df = pd.DataFrame()
    first_entries = set(first_dict.keys())
    second_entries = set(second_dict.keys())
    combined_entries = first_entries | second_entries
    difference_df['Annotation'] = sorted(combined_entries)
    difference_df[first_name+"_Raw_Tags"] = 0
    difference_df[second_name+"_Raw_Tags"] = 0
    difference_df[first_name+"_Normalized_Tags"] = 0
    difference_df[second_name+"_Normalized_Tags"] = 0
    difference_df['Raw_Difference'] = 0
    difference_df['Normalized_Difference'] = 0
    difference_df['Raw_FC'] = 0
    difference_df['Normalized_FC'] = 0
    for entry in sorted(combined_entries):
        if entry in first_entries:
            first_val = first_dict[entry]
            if entry not in second_entries:
                second_val = 0
            else:
                second_val = second_dict[entry]
        elif entry in second_entries:
            first_val = 0
            second_val = second_dict[entry]
        first_normalized_val = first_val/first_total
        second_normalized_val = (second_val/second_total)
        match = (difference_df['Annotation']==entry)
        difference_df.loc[match,first_name+"_Raw_Tags"] = first_val
        difference_df.loc[match,second_name+"_Raw_Tags"] = second_val
        difference_df.loc[match,first_name+"_Normalized_Tags"] = first_normalized_val
        difference_df.loc[match,second_name+"_Normalized_Tags"] = second_normalized_val
        difference_df.loc[match,"Raw_Difference"]=second_val-first_val
        difference_df.loc[match,"Normalized_Difference"]=second_normalized_val-first_normalized_val
        if first_val != 0:
            difference_df.loc[match,"Raw_FC"]=second_val/first_val
        if first_val != 0:
            difference_df.loc[match,"Normalized_FC"]=second_normalized_val/first_normalized_val
            
    return difference_df

name_regex = re.compile(r'(?:.*\/|^)(.*?)\.')
ext_regex = re.compile(r'(?:.*\/|^).*?\.(.*)')
path_regex = re.compile(r'(.*\/|^).*?\.')

master_dict = OrderedDict()
comparisons = set()

for output_file in snakemake.output:
    name_match = name_regex.search(output_file)
    ext_match = ext_regex.search(output_file)
    if name_match != None and ext_match != None:
        name = name_match.group(1)
        ext = ext_match.group(1)
    name_match = name_regex.search(output_file)
    if name_match != None:
        name = name_match.group(1)
        comparisons.add(name)
    else:
        print("Output file name {0} in an unexpected format, please try again!".format(output_file))
        
for input_file in snakemake.input:
    name_match = name_regex.search(input_file)
    path_match = path_regex.search(input_file)
    if name_match != None and path_match != None:
        name = name_match.group(1)
        ext = '.'.join([ext for ext in path_match.group(1).split('/') if ext != ''])
    else:
        print("Input file name {0} in an unexpected format, please try again!".format(input_file))
        sys.exit()
    df = pd.read_table(input_file)
    sample_dict = dict()
    for row in zip(df["Annotation"],df["Total_Tags"]):
        sample_dict[row[0]] = row[1]
    master_dict[(name,ext)] = sample_dict

compared = dict()

for first_sample in master_dict:
    ext = first_sample[1]
    for second_sample in master_dict:
        if second_sample[1] != ext:
            continue
        comparison = "{0}_vs_{1}".format(first_sample[0],second_sample[0])
        if comparison not in comparisons:
            continue
        else:
            comparison_string = comparison+ext
            compared[comparison_string] = compare(first_sample,second_sample)

for output_file in snakemake.output:
    Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
    name_match = name_regex.search(output_file)
    ext_match = ext_regex.search(output_file)
    if name_match != None and ext_match != None:
        name = name_match.group(1)
        ext = ext_match.group(1).strip('.txt')
        df = compared[name+ext]
        df.sort_values(by='Normalized_Difference', axis=0, ascending=False, inplace=True)
        f = open(output_file,'w')
        df.to_csv(f, sep="\t", line_terminator="\n", index=False)
        f.close()