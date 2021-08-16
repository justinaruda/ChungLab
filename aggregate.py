# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:36:06 2021

@author: justi
"""

import pandas as pd
import re
import os
from pathlib import Path
import matplotlib.pyplot as plt
        
def compare(first,second):
    global comparison_dict
    first_name = first
    second_name = second
    print("Comparing {0} and {1}...".format(first_name,second_name))
    first_dict = comparison_dict[first]
    second_dict = comparison_dict[second]
    difference_df = pd.DataFrame()
    first_entries = set(first_dict.keys())
    second_entries = set(second_dict.keys())
    combined_entries = first_entries | second_entries
    difference_df['Interval'] = sorted(combined_entries)
    difference_df[first_name+"_mfe"] = 0
    difference_df[second_name+"_mfe"] = 0
    difference_df['Raw_Difference'] = 0
    difference_df['Raw_FC'] = 0
    for entry in sorted(combined_entries):
        if entry in first_entries:
            first_val = float(first_dict[entry][0])
            if entry not in second_entries:
                second_val = 0
            else:
                second_val = float(second_dict[entry][0])
        elif entry in second_entries:
            first_val = 0
            second_val = second_dict[entry]
        match = (difference_df['Interval']==entry)
        difference_df.loc[match,first_name+"_mfe"] = first_val
        difference_df.loc[match,second_name+"_mfe"] = second_val
        difference_df.loc[match,"Raw_Difference"]=second_val-first_val
        if first_val != 0:
            difference_df.loc[match,"Raw_FC"]=second_val/first_val
        difference_df.sort_values(by='Raw_FC',inplace=True)
    return ("{0}_vs_{1}".format(first_name,second_name),difference_df)


def sort_idx(index):
    def sort(name):
        treatment_match = re.search(r'^[^_]*_([^_]*)_[^_]*',name)
        pos_match = re.search(r'^[^_]*_[^_]*_([^_]*?)\.',name)
        if treatment_match == None:
            return 1000
        if pos_match == None:
            return 1000
        else:
            treatment = treatment_match.group(1)
            pos = pos_match.group(1)
            string = "{a}_{b}".format(a=treatment,b=pos)
            order = {"Mock_Up":0,"IFN_Up":1,"Mock_Down":2,"IFN_Down":3}
            return order[string]
    list_sorted = sorted(index.tolist(),key=sort)
    new_idx = pd.Index(list_sorted)
    return new_idx

plt.figure(figsize=(10,12),dpi=600)

mfe_regex = re.compile(r'.* \(.*(-\d*\.\d*).*\).*')
name_regex= re.compile(r'\/(.*).fold')

Path(os.path.dirname(snakemake.output[0])).mkdir(parents=True, exist_ok=True)

output_dict = dict()
overall_df = pd.DataFrame()
comparison_dict = dict()

for input_file in snakemake.input:
    name = name_regex.search(input_file).group(1)
    rna_fold_output = pd.read_table(input_file,header=None)
    rna_fold_output.replace({mfe_regex : r'\1'},regex=True,inplace=True)
    rna_fold_regions = rna_fold_output[rna_fold_output[0].str.startswith(">")]
    rna_fold_mfes = rna_fold_output[rna_fold_output[0].str.startswith("-")]
    mean = rna_fold_mfes[0].astype(float).mean()
    overall_df[name] = rna_fold_mfes[0].astype(float)
    sample_dict = dict()
    for region,mfe in zip(rna_fold_regions[0],rna_fold_mfes[0]):
        sample_dict[region] = [mfe]
    f = open(os.path.join(os.path.dirname(snakemake.output[0]),name+'.out'),'w')
    sample_df = pd.DataFrame().from_dict(sample_dict,orient='index').to_csv(f,sep="\t",line_terminator="\n",index=True,header=False)
    f.close()
    comparison_dict[name] = sample_dict
    output_dict[input_file] = mean

for first_name in comparison_dict.keys():
    for second_name in comparison_dict.keys():
        if first_name.replace("unedited","edited") == second_name and second_name.replace("edited","unedited") == first_name:
            name,compared_df = compare(first_name,second_name)
            f = open(os.path.join(os.path.dirname(snakemake.output[0]),name+'_compared.out'),'w')
            compared_df.to_csv(f,sep="\t",line_terminator="\n",index=False,header=True)
            f.close()

df = pd.DataFrame().from_dict(output_dict,orient="index")

boxplot = overall_df.sort_index(axis=1,key=sort_idx).boxplot(notch=True)

plt.xticks(rotation = 45, ha='right')
plt.tight_layout()
plt.savefig('MFE_boxplots.svg',format='svg')

f = open(snakemake.output[0], 'w')
df.to_csv(f, sep="\t", line_terminator="\n", header=False)
f.close()