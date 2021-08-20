# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 23:02:35 2021

@author: justi
"""

import pandas as pd
import HTSeq
import re
import multiprocessing
import HTSeq

def remove_empty_sets(tuples):
    new_tuples = []
    for tup in tuples:
        if len(tup[1])!=0:
             new_tuples.append(tup)
    return new_tuples

def annotate(input_file,clip_df):
    global family_gas
    global gene_gas
    print("Annotating {0}...".format(input_file))
    clip_df['Annotation'] = 'None'
    clip_df['Repeat'] = 'None'
    clip_df['Repeat_Context'] = 'None'
    for idx,row in clip_df.iterrows():
        iv_match = iv_regex.search(row['Interval'])
        if iv_match is not None:
            region = iv_match.group(1)
            start = int(iv_match.group(2))
            end = int(iv_match.group(3))
            strand = iv_match.group(4)
            print(strand)
        else:
            continue
        clip_iv = HTSeq.GenomicInterval(region, start, end)
        crosslink_pos = HTSeq.GenomicPosition(region,start+round((end-start)/2),strand=strand)
        annot_set = gene_gas[crosslink_pos]
        if len(annot_set) != 0:
            annot_str = str(annot_set)
            clip_df.at[idx,'Annotation'] = annot_str
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
    return (input_file,clip_df)

annotation_file = snakemake.params.annotation
genome_annotation_file = snakemake.params.genome_annotation

clip_dfs = list()
print('Reading input files...')
for input_file in snakemake.input.files:
    clip_df = pd.read_table(input_file)
    clip_dfs.append((input_file,clip_df))

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


print('Loading genomic annotation file')

gene_gas = HTSeq.GenomicArrayOfSets("auto",stranded=True)
genomic_annotations = HTSeq.GFF_Reader(genome_annotation_file,end_included=True)
for annot in genomic_annotations:
    gene_gas[annot.iv] = annot.name

print('Starting annotation...')

iv_regex = re.compile(r'(chr\d*):(\d*)-(\d*).*\(([+-])\)')

with multiprocessing.Pool(processes=snakemake.threads) as pool:
        annotated_dfs = pool.starmap(annotate, clip_dfs)

for in_file,clip_df in annotated_dfs:
    out_file = '{0}.ir.out'.format(in_file.strip('.out'))
    print('Writing {0} to {1}...'.format(in_file,out_file))
    f = open(out_file, 'w')
    clip_df.to_csv(f, sep="\t", line_terminator="\n", index=False, header=True)
    f.close()