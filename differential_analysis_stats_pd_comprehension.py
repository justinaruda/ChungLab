# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 13:36:21 2021

@author: justi
"""

import pandas as pd
import numpy as np
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
    master_df = pd.DataFrame().from_dict(master_dict).fillna(0)
    return master_df

def compare(compare_string,first,second):
    global master_dict
    global stats_dict
    first_name,first_ext = first
    second_name,second_ext = second
    print("Comparing {0} and {1}...".format(first_name,second_name))
    first_dict = {first_name:master_dict[first]}
    first_stats = stats_dict[first]
    first_assigned_total = first_stats['Total_Assigned_Tags']
    first_assessed_total = first_stats['Total_Assessed_Tags']
    first_total = first_stats['Total_Tags']
    second_dict = {second_name:master_dict[second]}
    second_stats = stats_dict[second]
    second_assigned_total = second_stats['Total_Assigned_Tags']
    second_assessed_total = second_stats['Total_Assessed_Tags']
    second_total = second_stats['Total_Tags']

    header = ['{0}_Raw_Tags'.format(first_name),'{0}_Raw_Tags'.format(second_name),'Annotation']
    
    aligned_df = align_dicts([first_dict,second_dict])
    aligned_df['Annotation'] = aligned_df.index.values
    aligned_df.columns = header
    aligned_df = aligned_df[['Annotation','{0}_Raw_Tags'.format(first_name),'{0}_Raw_Tags'.format(second_name)]]
    
    aligned_df['{0}_Normalized_Tags_Assigned'.format(first_name)] = aligned_df['{0}_Raw_Tags'.format(first_name)]/first_assigned_total
    aligned_df['{0}_Normalized_Tags_Assigned'.format(second_name)] = aligned_df['{0}_Raw_Tags'.format(second_name)]/second_assigned_total
    aligned_df['{0}_Normalized_Tags_Total'.format(first_name)] = aligned_df['{0}_Raw_Tags'.format(first_name)]/first_total
    aligned_df['{0}_Normalized_Tags_Total'.format(second_name)] = aligned_df['{0}_Raw_Tags'.format(second_name)]/second_total
    aligned_df['{0}_Normalized_Tags_Assessed'.format(first_name)] = aligned_df['{0}_Raw_Tags'.format(first_name)]/first_assessed_total
    aligned_df['{0}_Normalized_Tags_Assessed'.format(second_name)] = aligned_df['{0}_Raw_Tags'.format(second_name)]/second_assessed_total
    aligned_df['Raw_Difference'] = aligned_df['{0}_Raw_Tags'.format(second_name)] - aligned_df['{0}_Raw_Tags'.format(first_name)]
    aligned_df['Normalized_Difference_Assigned'] = aligned_df['{0}_Normalized_Tags_Assigned'.format(second_name)] - aligned_df['{0}_Normalized_Tags_Assigned'.format(first_name)]
    aligned_df['Normalized_Difference_Total'] = aligned_df['{0}_Normalized_Tags_Total'.format(second_name)] - aligned_df['{0}_Normalized_Tags_Total'.format(first_name)]
    aligned_df['Normalized_Difference_Assessed'] = aligned_df['{0}_Normalized_Tags_Assessed'.format(second_name)] - aligned_df['{0}_Normalized_Tags_Assessed'.format(first_name)]
    aligned_df['Raw_FC'] = aligned_df['{0}_Raw_Tags'.format(second_name)]/aligned_df['{0}_Raw_Tags'.format(first_name)]
    aligned_df['Normalized_FC_Assigned'] = aligned_df['{0}_Normalized_Tags_Assigned'.format(second_name)]/aligned_df['{0}_Normalized_Tags_Assigned'.format(first_name)]
    aligned_df['Normalized_FC_Total'] = aligned_df['{0}_Normalized_Tags_Total'.format(second_name)]/aligned_df['{0}_Normalized_Tags_Total'.format(first_name)]
    aligned_df['Normalized_FC_Assessed'] = aligned_df['{0}_Normalized_Tags_Assessed'.format(second_name)]/aligned_df['{0}_Normalized_Tags_Assessed'.format(first_name)]
    
    aligned_df.replace(np.inf,0,inplace=True)
        
    return (compare_string, aligned_df)

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
    for i in range(9):
        compared_samples_df[cols[i]] = comparison_df.iloc[:,i]
    compared_samples_dict = compared_samples_df.to_dict()
    dict_list.append(compared_samples_dict)
summary_df = align_dicts(dict_list)
summary_df['Annotation'] = summary_df.index.values
summary_file = snakemake.params.summary
f = open(summary_file,'w')
summary_df.to_csv(f, sep="\t", line_terminator="\n", index=False)
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