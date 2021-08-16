# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:45:26 2021

@author: justi
"""

import pandas as pd
import re
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
        
def align_dicts(dicts):
    master_dict = dict()
    for single_dict in dicts:
        master_dict = {**master_dict,**single_dict}
    master_df = pd.DataFrame().from_dict(master_dict).fillna(0)
    return master_df

def compare(first,second):
    global comparison_dict
    first_name = first
    second_name = second
    print("Comparing {0} and {1}...".format(first_name,second_name))
    first_dict = comparison_dict[first]
    second_dict = comparison_dict[second]
    difference_df = align_dicts([first_dict,second_dict])
    difference_df.columns = ['Interval','{0}_mfe'.format(first_name),'{0}_mfe'.format(second_name)]
    difference_df['Interval'] = difference_df.index.values
    difference_df['Raw_Difference'] = difference_df['{0}_mfe'.format(second_name)]-difference_df['{0}_mfe'.format(first_name)]
    difference_df['Raw_FC'] = difference_df['{0}_mfe'.format(second_name)]/difference_df['{0}_mfe'.format(first_name)]
    difference_df.replace(np.inf,99999)
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