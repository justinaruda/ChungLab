# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 08:25:09 2021

@author: justi
"""

import pandas as pd
import re
import csv

expressed_df = pd.read_table(r'/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/expressed_genes.csv',header=None)
expressed_series = pd.Series(expressed_df[0])
name_set = set()
for element in expressed_series:
    for name in element.split(" "):
        name_set.add(name)

gtf_df = pd.read_table(r'/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/gencode.v38.primary_assembly.annotation.gene_name.gtf',header=None)

name_regex = re.compile(r'gene_name "(.*?)"')

total = len(gtf_df.index)

names = []    
for values in gtf_df[8]:
    names.append(name_regex.search(values).group(1))

gtf_df['Gene_Name'] = names
gtf_df = gtf_df.loc[gtf_df['Gene_Name'].isin(name_set)]
gtf_df.drop('Gene_Name',axis=1,inplace=True)

with open(r'/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/gencode.v38.primary_assembly.annotation.gene_name_expressed.gtf','w') as f:
    gtf_df.to_csv(f, sep="\t", line_terminator="\n", index=False, header=False, quoting=csv.QUOTE_NONE)