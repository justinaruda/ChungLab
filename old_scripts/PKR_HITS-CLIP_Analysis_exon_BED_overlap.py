# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 11:43:41 2021

@author: justi
"""
import HTSeq
import pysam
import pandas as pd
import numpy as np
import getopt, sys
import io
import csv
import re

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]

short_options = "b:a:o:p:q:v"
long_options = ["bam=", "bed=", "annotation_file=", "output=", "peak_file=", "mapq=", "verbose"]

try:
    arguments, values = getopt.getopt(argument_list, short_options, long_options)
except getopt.error as err:
    # Output error, and return with an error code
    print (str(err))
    sys.exit(2)

table_file = None
annotation_file = None
out_file = None
mapq = 60
peak_file = None
no_peaks = False
verbose = False

for current_argument, current_value in arguments:
    if current_argument in ("-b", "--bam", "--bed"):
        clip_file = current_value
        
    if current_argument in ("-a", "--annotation_file"):
        annotation_file = current_value
    if current_argument in ("-o", "--output"):
        out_file = current_value
    if current_argument in ("-p", "--peak_fie"):
        peak_file = current_value
    if current_argument in ("-q", "--mapq"):
        mapq = int(current_value)
    if current_argument in ("-v", "--verbose"):
        verbose = True
        
if peak_file == None:
    no_peaks = True
    print('Running in no peaks mode')

if clip_file.endswith('.bam'):
    clip_type = 'bam'
if clip_file.endswith('.gz'):
    clip_type = 'bed'

if clip_type == 'bam':
    print('MAPQ minimum set to '+str(mapq)+'.')
elif clip_type == 'bed':
    print('CLIP tags detected in .BED format. Ignoring MAPQ...')
print('Reading input files...')

# peak_file = r'Peaks_Uniq_13T.gtf'
# #input_bam = r'/home/hchunglab/data/justin/PKR_fCLIP/SRR6456378Aligned.sortedByCoord.out.bam'
# clip_file = r'/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/Align/13T/Peaks_uniq_dedup_13T.sorted.out.bam'
# annotation_file = r'/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/hg38.ncbiRefSeq.measles.gtf'
# #annotation_file = r'repeat_masker_grch38.fix.bed' #r'Alu_elements.bed' #r'alu_peak_intersect.bed' # 
# out_file = 'test_output.txt'

def get_intersection(iv1,iv2):
    start = iv1.start if iv1.start > iv2.start else iv2.start
    end = iv1.end if iv1.end < iv2.end else iv2.end
    return HTSeq.GenomicInterval(iv1.chrom,start,end,strand=iv1.strand)

def drop_peaks(label,iv):
    
    global gene_counts   
    
    gene_name = re.match(r'(.*)_[^_]*',label).group(1)
    feature_type = re.match(r'.*_([^_]*)',label).group(1)
    
    if feature_type == 'gene':
        new_entry_list = []
        gene_label = gene_name+"_gene"
        i = 0
        for entry in gene_counts[gene_label]:
            if entry[0] == iv:
                i+=1
            if i > 1:
                continue
            new_entry_list.append(entry)
        gene_counts[gene_label] = new_entry_list
    
    order = ['UTR','CDS','exon','transcript']
    position = [i for i, feature in enumerate(order) if feature in feature_type]
    if len(position)==0:
        return
    position = position[0]    
    to_del = set()
    dropped_from = set()
    for entry_label in gene_counts.keys():
        entry_name = re.match(r'(.*)_[^_]*',entry_label).group(1)
        entry_type = re.match(r'.*_([^_]*)',entry_label).group(1)
        if entry_label in dropped_from:
            continue
        if gene_name == entry_name:
            i = 0
            drop = False
            while i < len(order):
                if order[i] in entry_type:
                    if not drop and i < position:
                        i = position
                        drop = True
                    if not drop and i > position:
                        drop = True
                    if drop:
                        new_entry_list = []
                        current_label = gene_name+'_'+order[i]
                        for entry in gene_counts[current_label]:
                            old_iv = entry[0]
                            if iv.overlaps(old_iv):
                                new_iv = get_intersection(iv,old_iv)
                                count_tags(new_iv)
                                continue
                            new_entry_list.append(entry)
                        if len(new_entry_list) == 0:
                            to_del.add(current_label)
                        gene_counts[current_label] = new_entry_list
                        dropped_from.add(current_label)
                i+=1
    for key in to_del:
        print(to_del)
        del gene_counts[key]                       

def remove_empty_sets(tuples):
    new_tuples = []
    for tup in tuples:
        if len(tup[1])!=0:
             new_tuples.append(tup)
    return new_tuples

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0])) != 0

def process_cigar(clip_tag, start, end):
    for cigar_operation in clip_tag.cigar:
        if cigar_operation.type == 'M':  
            if getOverlap([cigar_operation.ref_iv.start,cigar_operation.ref_iv.end], [start,end]):
                return 1
    return 0

def get_header(file_name):
    if file_name.endswith('.gtf'):
        header = ['Region','Source','Type','Start','End','Score','Strand','Frame','Attributes']
    elif file_name.endswith('.bed'):
        header = ['Region','Start','End','Attributes','Count','Strand']
    return header

# def divide_feature:
#     global annots:
#         for feature in annots.fetch()

def get_peaks(feature):
    global no_peaks
    if no_peaks:
        return [(feature.iv,feature.iv.chrom + ':' + str(feature.iv.start) + '-' + str(feature.iv.end))]
    else:
        return [iv_id_tuple for iv_id_tuple in gas[feature.iv].steps()]

if no_peaks:
    peaks_df = pd.DataFrame()
else:
    peaks_df = pd.read_table(peak_file, header=None, names=get_header(peak_file))

gas = HTSeq.GenomicArrayOfSets('auto',stranded=True)

for idx, peak in peaks_df.iterrows():
    if '_' in peak['Region']:
        peak['Region'] = peak['Region'].replace("_","")
    iv = HTSeq.GenomicInterval(peak['Region'],peak['Start'],peak['End'],strand=peak['Strand'])
    gas[iv] += peak['Region'] + ':' + str(peak['Start']) + '-' + str(peak['End'])

gtf_df = pd.read_table(annotation_file, header=None, names=get_header(annotation_file))

gtf_df['Region'].replace(regex="_",value="",inplace=True)

gtf_df = gtf_df[~gtf_df['Region'].str.startswith('#')]

# for idx, feature in gtf_df:
#     gtf_df['Region'].iloc[idx] = feature['Region'].replace("_","")

gtf_io = io.StringIO()

gtf_df.to_csv(gtf_io, sep="\t", line_terminator="\n", index=False, header=False, quoting=csv.QUOTE_NONE)

gtf_io.seek(0)

if annotation_file.endswith('.gtf'):
    annots = HTSeq.GFF_Reader(gtf_io,end_included=True)
    annotation_type = 'gtf'
if annotation_file.endswith('.bed'):
    annots = HTSeq.BED_Reader(gtf_io)
    annotation_type = 'bed'

annots_df = pd.read_table(annotation_file, header=None, names=get_header(annotation_file))

#input_bam_reader = HTSeq.BAM_Reader(input_bam)

if clip_type == 'bam':
    clip_reader = HTSeq.BAM_Reader(clip_file)
if clip_type == 'bed':
    clip_reader = pysam.TabixFile(clip_file)

gene_counts = {}

i = 0
l = len(pd.read_table(annotation_file,header=None).index)
checked_features = set()
if no_peaks:
    chroms = set(annots_df['Region'].values)
else:
    chroms = set(peaks_df['Region'].values)
checked_region_ivs = {}
checked_region_reads = {}


def count_tags(iv):
    global clip_type
    global mapq
    global label
    global gene_counts
    global clip_read_names
    global checked_region_ivs
    global checked_region_reads
    global annotation_type
    global clip_reader
    global clip_tags
    
    chrom = iv.chrom
    start = iv.start
    end = iv.end

    if clip_type == 'bam':
        clip_tags = clip_reader.fetch(region=chrom+':'+str(start)+'-'+str(end))
    if clip_type == 'bed' and chrom in clip_reader.contigs:
        clip_tags = clip_reader.fetch(chrom,start,end,parser=pysam.asBed())
#        input_count = 0
    clip_count = 0
#        input_read_names = set()
    # for clip_tag in input_reads:
    #     if clip_tag.read.name in input_read_names:
    #         continue
    #     if clip_tag.aQual < 60:
    #         continue
    #     if clip_tag.aligned:
    #         input_count += process_cigar(clip_tag,start,end)
    for clip_tag in clip_tags:
        if clip_type == 'bam':
            tag_name = clip_tag.read.name
        elif clip_type == 'bed':
            tag_name = clip_tag.name
        if tag_name in clip_read_names:
            continue
        if clip_type == 'bam' and clip_tag.aQual < mapq:
            continue
        if clip_type == 'bam' and clip_tag.aligned:
            clip_count += process_cigar(clip_tag,start,end)
            clip_read_names.add(tag_name)
        elif clip_type == 'bed':
            clip_count += 1
            clip_read_names.add(tag_name)
    if clip_count == 0:
        return True
    if label not in gene_counts.keys():
        gene_counts[label] = [(HTSeq.GenomicInterval(chrom,start,end),clip_count)] #[clip_count-input_count] #[clip_count/input_count] for FC
        checked_region_ivs[label] = checked_ivs
        checked_region_reads[label] = clip_read_names
    else:
        gene_counts[label].append((HTSeq.GenomicInterval(chrom,start,end),clip_count)) #.append(clip_count-input_count) #.append(clip_count/input_count)
        checked_region_ivs[label] = set.union(checked_region_ivs[label],checked_ivs)
        checked_region_reads[label] = set.union(checked_region_reads[label],clip_read_names)
        
    if annotation_type == 'gtf':
        
        if gene_label not in gene_counts.keys():
            gene_counts[gene_label] = [(HTSeq.GenomicInterval(chrom,start,end),clip_count)]
        else:
            gene_counts[gene_label].append((HTSeq.GenomicInterval(chrom,start,end),clip_count))
    return False

for feature in annots:
    i+=1
    if feature.name+'_'+feature.iv.chrom+':'+str(feature.iv.start)+'-'+str(feature.iv.end) not in checked_features:
        if verbose:
            print('Processing {0}, {1} of {2} ({3:.2f}%)'.format(feature.name,i,l,(i/l)*100))
        elif i%10000 == 0:
            print('Processing {0}, {1} of {2} ({3:.2f}%)'.format(feature.name,i,l,(i/l)*100))
        checked_ivs = set()
        clip_read_names = set()
        checked_features.add(feature.name+'_'+feature.iv.chrom+':'+str(feature.iv.start)+'-'+str(feature.iv.end))
    else:
        continue
    
    if annotation_type == 'bed':
        feature.iv.strand = '+'
        peak_ivs = get_peaks(feature)
        feature.iv.strand = '-'
        peak_ivs = peak_ivs + get_peaks(feature)
    if annotation_type == 'gtf':
        peak_ivs = get_peaks(feature)
        
    peak_ivs = remove_empty_sets(peak_ivs)
    
    if annotation_type == 'gtf':
        label = feature.name+'_'+feature.type
    if annotation_type == 'bed':
        label = feature.name
    if label in gene_counts.keys():
        checked_ivs = checked_region_ivs[label]
        clip_read_names = checked_region_reads[label]
        
    gene_label = feature.name+'_'+'gene'
    
    for iv,val in peak_ivs:
        
        chrom = iv.chrom
        if chrom not in chroms:
            continue
        start = feature.iv.start if feature.iv.start > iv.start else iv.start
        end = feature.iv.end if feature.iv.end < iv.end else iv.end
        
        new_iv = HTSeq.GenomicInterval(iv.chrom,start,end,strand=iv.strand)
        
        if new_iv in checked_ivs:
            continue
        
        checked_ivs.add(new_iv)
        
#        input_reads = input_bam_reader.fetch(region=chrom+':'+str(start)+'-'+str(end))
        if count_tags(new_iv):
            continue
        
        drop_peaks(label,(HTSeq.GenomicInterval(chrom,start,end)))
        drop_peaks(gene_label,(HTSeq.GenomicInterval(chrom,start,end)))

final_gene_count = {}

gene_list = list(gene_counts.keys())

final_gene_count_df = pd.DataFrame()
final_gene_count_df['Gene'] = gene_list
final_gene_count_df['Peak_Count'] = 0
final_gene_count_df['Average_Tags'] = 0
final_gene_count_df['Total_Tags'] = 0
final_gene_count_df['Peak_Tags'] = "" #FC

for gene in gene_list:
    
    gene_count_array = [tup[1] for tup in gene_counts[gene]]
    if len(gene_count_array) > 0:
        final_gene_count_df.loc[(final_gene_count_df['Gene']==gene),'Peak_Count'] = len(gene_count_array)
        final_gene_count_df.loc[(final_gene_count_df['Gene']==gene),'Average_Tags'] = np.average(gene_count_array) #FC
        final_gene_count_df.loc[(final_gene_count_df['Gene']==gene),'Total_Tags'] = np.sum(gene_count_array)
        final_gene_count_df.loc[(final_gene_count_df['Gene']==gene),'Peak_Tags'] = '['+', '.join(map(str,[(str(tup[0]),tup[1]) for tup in gene_counts[gene]]))+']'

final_gene_count_df.sort_values(by='Total_Tags', axis=0, ascending=False, inplace=True)
f = open(out_file, 'w')
final_gene_count_df.to_csv(f, sep="\t", line_terminator="\n", index=False, header=True)
f.close()