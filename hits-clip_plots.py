# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 20:39:42 2021

@author: justi
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import mplcursors
import re
import matplotlib.cm as cm

table_df = pd.read_table(r'C:/Users/justi/Desktop/PKR HITS-CLIP/CTK/differential_analysis_gencode/snakemake/new/aggregated_sample_summary.txt')
#table_df = pd.read_table(r'C:/Users/justi/Desktop/PKR HITS-CLIP/CTK/CITS/snakemake_gencode/aggregated_sample_summary.txt')

name_regex = re.compile('(.*)_[^_]*$')

names = []
for annot in table_df['Annotation']:
    name_match = name_regex.search(annot)
    if name_match is not None:
        names.append(name_match.group(1))
    else:
        print(annot)
table_df['Name'] = names

genes = table_df[table_df['Annotation'].str.contains('_gene')]
exons = table_df[table_df['Annotation'].str.contains('_exon')]
utrs = table_df[table_df['Annotation'].str.contains('_UTR')]
introns = table_df[table_df['Annotation'].str.contains('_transcript')]
sno_rnas = table_df[table_df['Annotation'].str.contains('SNO')]

def stacked_bar(input_file,libraries,groups):
    input_df = pd.read_table(input_file)
    for library in libraries:
        genes = input_df[input_df['Annotation'].str.contains('_gene')]
        exons = input_df[input_df['Annotation'].str.contains('_exon')]
        utrs = input_df[input_df['Annotation'].str.contains('_UTR')]
        introns = input_df[input_df['Annotation'].str.contains('_transcript')]
        
        
        
    

add_groups = ['SNHG11', 'SNHG12', 'SNHG15', 'SNHG16', 'SNHG17', 'SNHG19', 'SNHG25', 'SNHG26', 'SNHG29', 'SNHG3', 'SNHG30', 'SNHG4', 'SNHG5', 'SNHG6', 'SNHG7', 'SNHG8', 'SNHG1', 'GAS5', 'RPS12', 'RPS2', 'RPS8', 'RPSA', 'RPL17', 'RPL23A', 'RPL4', 'RPL7A', 'GNL3', 'RCC1', 'RPL13A', 'RPL27A']

def scatter_plot(group1,group2,key,groups):
    global table_df
    
    add_dict = dict()
    parsed_groups = []
    for group in groups:
        if type(group) == tuple:
            add_dict[group[0]] = group[1]
            group = group[0]
        parsed_groups.append(group)
    
    colors = cm.viridis(np.linspace(0, 1, len(groups)+1))[1:]
    
    table_df = table_df[table_df['Annotation'].str.contains(key)]
    plt.figure()
    ax = plt.axes()
    
    leftover_df = table_df
   
    scatters = dict()
    for series_key,color in zip(parsed_groups,colors):
        if series_key in add_dict.keys():
            cond = table_df['Annotation'].str.contains(series_key) | table_df['Name'].isin(add_dict[series_key])
            leftover_cond = leftover_df['Annotation'].str.contains(series_key) | leftover_df['Name'].isin(add_dict[series_key])
        else:
            cond = table_df['Annotation'].str.contains(series_key)
            leftover_cond = leftover_df['Annotation'].str.contains(series_key)
        
        series_df = table_df[cond]
        leftover_df = leftover_df[~leftover_cond]
        x = series_df[group1]
        y = series_df[group2]
        labels = list(series_df['Name'])
        scatters[series_key] = (ax.scatter(x,y,color=color,zorder=10,label=series_key),labels)
    
    mplcursors.cursor(scatters['SNORD'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['SNORD'][1][sel.target.index]))
    mplcursors.cursor(scatters['SNORA'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['SNORA'][1][sel.target.index]))
    mplcursors.cursor(scatters['MT-'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['MT-'][1][sel.target.index]))
    mplcursors.cursor(scatters['Embedded snoRNA'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['Embedded snoRNA'][1][sel.target.index]))
    
    x = leftover_df[group1]
    y = leftover_df[group2]
    labels = list(leftover_df['Name'])
    scatters['leftover'] = ax.scatter(x,y,color='#1c1c1c',label='Other'),labels
    ax.legend()
    
    ax.set_xlabel(group1)
    ax.set_ylabel(group2)
    
    mplcursors.cursor(scatters['SNORD'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['SNORD'][1][sel.target.index]))
    mplcursors.cursor(scatters['SNORA'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['SNORA'][1][sel.target.index]))
    mplcursors.cursor(scatters['MT-'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['MT-'][1][sel.target.index]))
    mplcursors.cursor(scatters['leftover'][0],multiple=True).connect("add", lambda sel: sel.annotation.set_text(scatters['leftover'][1][sel.target.index]))

    
    lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),
    np.max([ax.get_xlim(), ax.get_ylim()]),
    ]
    
    lims = np.array(lims)

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.plot(lims, lims*2, color='firebrick', linestyle = 'dotted', alpha=0.5, zorder=10)
    ax.plot(lims, lims*0.5, color='firebrick', linestyle = 'dotted', alpha=0.5, zorder=10)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    plt.title('{0} vs. {1}'.format(group1,group2))

#scatter_plot('A1KO_Mock_Down_Normalized_Tags_Total','A1KO_IFN_Down_Normalized_Tags_Total','_gene',[("Embedded snoRNA",add_groups),'SNORD','SNORA','MT-'])
#scatter_plot('A1KO_Mock_Down_Normalized_Tags_Assigned','A1KO_IFN_Down_Normalized_Tags_Assigned','_gene',[("Embedded snoRNA",add_groups),'SNORD','SNORA','MT-'])
#scatter_plot('A1KO_Mock_Down_Normalized_Tags_Assessed','A1KO_IFN_Down_Normalized_Tags_Assessed','_gene',[("Embedded snoRNA",add_groups),'SNORD','SNORA','MT-'])