# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:36:06 2021

@author: justi
"""

import pandas as pd
import re
import os
from pathlib import Path

mfe_regex = re.compile(r'.* \(.*(-\d*\.\d*).*\).*')
name_regex= re.compile(r'\/(.*).fold')

Path(os.path.dirname(snakemake.output[0])).mkdir(parents=True, exist_ok=True)

output_dict = dict()

for input_file in snakemake.input:
    rna_fold_output = pd.read_table(input_file,header=None)
    rna_fold_output.replace({mfe_regex : r'\1'},regex=True,inplace=True)
    rna_fold_output = rna_fold_output[rna_fold_output[0].str.startswith("-")]
    mean = rna_fold_output[0].astype(float).mean()
    f = open(os.path.join(os.path.dirname(snakemake.output[0]),name_regex.search(input_file).group(1)+'.out'),'w')
    rna_fold_output[0].to_csv(f,index=False,header=False)
    f.close()
    output_dict[input_file] = mean

df = pd.DataFrame().from_dict(output_dict,orient="index")

f = open(snakemake.output[0], 'w')
df.to_csv(f, sep="\t", line_terminator="\n", header=False)
f.close()