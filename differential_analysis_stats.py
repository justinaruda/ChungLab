# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 11:42:11 2021

@author: justi
"""

import pandas as pd
from collections import OrderedDict
import re
import sys
import os
from pathlib import Path
import multiprocessing

def parse_stats(stats_file):
    global total_tags_regex
    global unique_assigned_tags_regex
    global multiple_assigned_tags_regex
    global total_assigned_tags_regex
    with open(stats_file) as f:
        stats = f.read()
        return {"Total_Tags":int(total_tags_regex.search(stats).group(1)),"Total_Assessed_Tags":int(total_assessed_tags_regex.search(stats).group(1)),"Unique_Assigned_Tags":int(unique_assigned_tags_regex.search(stats).group(1)),"Multiple_Assigned_Tags":int(multiple_assigned_tags_regex.search(stats).group(1)),"Total_Assigned_Tags":int(total_assigned_tags_regex.search(stats).group(1))}

def align_dicts(dicts):
    master_dict = dict()
    for single_dict in dicts:
        master_dict = {**master_dict,**single_dict}
    master_df = pd.DataFrame(master_dict).fillna(int(0))
    master_dict = master_df.to_dict(orient='index')
    return master_dict

def compare(compare_string,first,second):
    global master_dict
    global stats_dict
    first_name,first_ext = first
    second_name,second_ext = second
    print("Comparing {0} and {1}...".format(first_name,second_name))
    first_dict = master_dict[first]
    first_stats = stats_dict[first]
    first_assigned_total = first_stats['Total_Assigned_Tags']
    first_total = first_stats['Total_Tags']
    second_dict = master_dict[second]
    second_stats = stats_dict[second]
    second_assigned_total = second_stats['Total_Assigned_Tags']
    second_total = second_stats['Total_Tags']
    difference_df = pd.DataFrame()
    first_entries = set(first_dict.keys())
    second_entries = set(second_dict.keys())
    combined_entries = first_entries | second_entries
    difference_df['Annotation'] = sorted(combined_entries)
    difference_df[first_name+"_Raw_Tags"] = 0
    difference_df[second_name+"_Raw_Tags"] = 0
    difference_df[first_name+"_Normalized_Tags_Assigned"] = 0
    difference_df[second_name+"_Normalized_Tags_Assigned"] = 0
    difference_df[first_name+"_Normalized_Tags_Total"] = 0
    difference_df[second_name+"_Normalized_Tags_Total"] = 0
    difference_df['Raw_Difference'] = 0
    difference_df['Normalized_Difference_Assigned'] = 0
    difference_df['Normalized_Difference_Total'] = 0
    difference_df['Raw_FC'] = 0
    difference_df['Normalized_FC_Assigned'] = 0
    difference_df['Normalized_FC_Total'] = 0
    
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
        
        first_normalized_assigned = first_val/first_assigned_total
        first_normalized_total = first_val/first_total
        second_normalized_assigned = second_val/second_assigned_total
        second_normalized_total = second_val/second_total
        match = (difference_df['Annotation']==entry)
        difference_df.loc[match,first_name+"_Raw_Tags"] = first_val
        difference_df.loc[match,second_name+"_Raw_Tags"] = second_val
        difference_df.loc[match,first_name+"_Normalized_Tags_Assigned"] = first_normalized_assigned
        difference_df.loc[match,second_name+"_Normalized_Tags_Assigned"] = second_normalized_assigned
        difference_df.loc[match,first_name+"_Normalized_Tags_Total"] = first_normalized_total
        difference_df.loc[match,second_name+"_Normalized_Tags_Total"] = second_normalized_total
        difference_df.loc[match,"Raw_Difference"]=second_val-first_val
        difference_df.loc[match,"Normalized_Difference_Assigned"]=second_normalized_assigned-first_normalized_assigned
        difference_df.loc[match,"Normalized_Difference_Total"]=second_normalized_total-first_normalized_total
        if first_val != 0:
            difference_df.loc[match,"Raw_FC"]=second_val/first_val
            difference_df.loc[match,"Normalized_FC_Assigned"]=second_normalized_assigned/first_normalized_assigned
            difference_df.loc[match,"Normalized_FC_Total"]=second_normalized_total/first_normalized_total
            
    return (compare_string, difference_df)

name_regex = re.compile(r'(?:.*\/|^)(.*?)\.')
ext_regex = re.compile(r'(?:.*\/|^).*?\.(.*)')
path_regex = re.compile(r'(.*\/|^).*?\.')
total_tags_regex = re.compile(r'(?:Total_Library_Tags:) (\d*)')
total_assessed_tags_regex = re.compile(r'(?:Total_Assessed_Tags:) (\d*)')
unique_assigned_tags_regex = re.compile(r'(?:Unique_Assigned_Tags:) (\d*)')
multiple_assigned_tags_regex = re.compile(r'(?:Multiple_Assigned_Tags:) (\d*)')
total_assigned_tags_regex = re.compile(r'(?:Total_Assigned_Tags:) (\d*)')

master_dict = OrderedDict()
stats_dict = OrderedDict()
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
        
    stats_file = '{0}.stats'.format(input_file)
    
    df = pd.read_table(input_file)
    sample_dict = dict()
    for row in zip(df["Annotation"],df["Total_Tags"]):
        sample_dict[row[0]] = row[1]
    master_dict[(name,ext)] = sample_dict
    stats_dict[(name,ext)] = parse_stats(stats_file)
    
compared = dict()
to_compare = set() 

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
            to_compare.add((comparison_string,first_sample,second_sample))

with multiprocessing.Pool(processes=snakemake.threads) as pool:
        compared_list = pool.starmap(compare, to_compare)
        
dict_list = list()
for comparison_string, comparison_df in compared_list:
    compared[comparison_string] = comparison_df
    cols = comparison_df.columns
    compared_samples_df = pd.DataFrame()
    for i in range(6):
        compared_samples_df[cols[i]] = comparison_df.iloc[:,i]
    compared_samples_dict = compared_samples_df.to_dict(orient='index')
    dict_list.append(compared_samples_dict)

summary_df = pd.DataFrame().from_dict(align_dicts(dict_list),orient='index')

summary_file = snakemake.params.summary
f = open(summary_file,'w')
summary_df.to_csv(f, sep="\t", line_terminator="\n", index=True)
f.close()

for output_file in snakemake.output:
    Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
    name_match = name_regex.search(output_file)
    ext_match = ext_regex.search(output_file)
    if name_match != None and ext_match != None:
        name = name_match.group(1)
        ext = ext_match.group(1).strip('.txt')
        df = compared[name+ext]
        df.sort_values(by='Normalized_Difference_Total', axis=0, ascending=False, inplace=True)
        f = open(output_file,'w')
        df.to_csv(f, sep="\t", line_terminator="\n", index=False)
        f.close()