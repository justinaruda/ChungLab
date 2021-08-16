# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 16:36:00 2021

@author: justi
"""

import pandas as pd
import numpy as np
from collections import OrderedDict

def parse_name(name):
    
    treatments = ["Mock","IFN"]
    groups = ["Up","Down"]
    sample_dict = OrderedDict()
    for treatment in treatments:
        if treatment in name:
            sample_dict["treatment"] = treatment
    for group in groups:
        if group in name:
            sample_dict["group"] = group
    
    return sample_dict

def compare(first,second):
    global master_dict
    print("Comparing {0} and {1}".format(first,second))
    first_dict = master_dict[first]
    first_total = sum(list(first_dict.values()))
    second_dict = master_dict[second]
    second_total = sum(list(second_dict.values()))
    difference_df = pd.DataFrame()
    first_entries = set(first_dict.keys())
    second_entries = set(second_dict.keys())
    combined_entries = first_entries | second_entries
    difference_df['Annotation'] = sorted(combined_entries)
    difference_df[first+"_Raw_Tags"] = 0
    difference_df[second+"_Raw_Tags"] = 0
    difference_df[first+"_Normalized_Tags"] = 0
    difference_df[second+"_Normalized_Tags"] = 0
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
        
        difference_df.loc[difference_df['Annotation']==entry,first+"_Raw_Tags"] = first_val
        difference_df.loc[difference_df['Annotation']==entry,second+"_Raw_Tags"] = second_val
        difference_df.loc[difference_df['Annotation']==entry,first+"_Normalized_Tags"] = first_normalized_val
        difference_df.loc[difference_df['Annotation']==entry,second+"_Normalized_Tags"] = second_normalized_val
        difference_df.loc[difference_df['Annotation']==entry,"Raw_Difference"]=second_val-first_val
        difference_df.loc[difference_df['Annotation']==entry,"Normalized_Difference"]=first_normalized_val-second_normalized_val
        if first_val != 0:
            difference_df.loc[difference_df['Annotation']==entry,"Raw_FC"]=second_val/first_val
        if first_val != 0:
            difference_df.loc[difference_df['Annotation']==entry,"Normalized_FC"]=second_normalized_val/first_normalized_val
            
    return difference_df

comparisons = set()
file_names = []
master_dict = OrderedDict()
for f in ["A1KO_Mock_Down","A1KO_Mock_Up","A1KO_IFN_Up","A1KO_IFN_Down"]:
    file_name=(r"C:\Users\justi\Desktop\PKR HITS-CLIP\CTK\quantified\refSeq\\"+f+r".tags.quantified.txt")
    df = pd.read_table(file_name)
    sample_dict = dict()
    for row in zip(df["Annotation"],df["Total_Tags"]):
        sample_dict[row[0]] = row[1]
    master_dict[f] = sample_dict

compared = []

for first_sample in master_dict:
    first_sample_info = parse_name(first_sample)
    for second_sample in master_dict:
        second_sample_info = parse_name(second_sample)
        if first_sample_info["treatment"] != second_sample_info["treatment"] and first_sample_info["group"] == second_sample_info["group"]:
            comparison = str(sorted([first_sample,second_sample]))
            if comparison in comparisons:
                continue
            else:
                comparisons.add(comparison)
                comparison_string = "{0}_vs_{1}".format(first_sample,second_sample)
                compared.append((comparison_string,compare(first_sample,second_sample)))

for comparison,df in compared:
    df.sort_values(by='Normalized_Difference', axis=0, ascending=False, inplace=True)
    f = open(r"C:\Users\justi\Desktop\PKR HITS-CLIP\CTK\quantified\refSeq\\"+comparison+r".txt",'w')
    df.to_csv(f, sep="\t", line_terminator="\n", index=False)
    f.close()