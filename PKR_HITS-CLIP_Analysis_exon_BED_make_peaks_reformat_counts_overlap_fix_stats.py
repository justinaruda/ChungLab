# -*- coding: utf-8 -*-
"""
Created on Sun Aug  8 13:05:39 2021

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
import gzip

full_cmd_arguments = sys.argv
argument_list = full_cmd_arguments[1:]

short_options = "b:a:o:p:q:x:svzl"
long_options = ["bam=", "bed=", "annotation_file=", "output=", "peak_file=", "mapq=", "bed_output", "strict", "verbose", "keep_zeros", "stats"]

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
strict = False
bed_output = None
out_bed = False
keep_zeros = False
stats = False

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
    if current_argument in ("-x", "--bed_output"):
        bed_output = current_value
        out_bed = True
    if current_argument in ("-s", "--strict"):
        strict = True
    if current_argument in ("-v", "--verbose"):
        verbose = True
    if current_argument in ("-z", "--keep_zeros"):
        keep_zeros = True
    if current_argument in ("-l", "--stats"):
        stats = True
        
if peak_file == None:
    no_peaks = True
    print('Running in no peaks mode.')
if strict:
    print('Running in strict mode.')
if keep_zeros:
    print('Writing intervals with zero tags.')
if clip_file.endswith('.bam'):
    clip_type = 'bam'
if clip_file.endswith('.gz'):
    clip_type = 'bed'

if clip_type == 'bam':
    print('MAPQ minimum set to '+str(mapq)+'.')
elif clip_type == 'bed':
    print('CLIP tags detected in .BED format. Ignoring MAPQ...')
print('Reading input files...')

if out_bed:
    bed_list = []

# peak_file = r'Peaks_Uniq_13T.gtf'
# #input_bam = r'/home/hchunglab/data/justin/PKR_fCLIP/SRR6456378Aligned.sortedByCoord.out.bam'
# clip_file = r'/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/Align/13T/Peaks_uniq_dedup_13T.sorted.out.bam'
# annotation_file = r'/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/hg38.ncbiRefSeq.measles.gtf'
# #annotation_file = r'repeat_masker_grch38.fix.bed' #r'Alu_elements.bed' #r'alu_peak_intersect.bed' # 
# out_file = 'test_output.txt'

def check_label(label):
    global gene_counts
    global name_regex
    global feature_regex
    gene_name = label[0]
    feature_type = label[1]
    if gene_name not in gene_counts.keys():
        return False
    if feature_type in gene_counts[gene_name].keys():
        return True
    else:
        return False

def get_label(label):
    global gene_counts
    global name_regex
    global feature_regex
    gene_name = label[0]
    feature_type = label[1]
    if feature_type in gene_counts[gene_name].keys():
        return gene_counts[gene_name][feature_type]
    
def set_label(label,value):
    global gene_counts
    global name_regex
    global feature_regex
    gene_name = label[0]
    feature_type = label[1]
    if gene_name not in gene_counts.keys():
        gene_counts[gene_name] = {feature_type:value}
    else:
        gene_counts[gene_name][feature_type] = value
    
def delete_label(label):
    global gene_counts
    global name_regex
    global feature_regex
    gene_name = label[0]
    feature_type = label[1]
    del gene_counts[gene_name][feature_type]

def append_label(label,value):
    global gene_counts
    global name_regex
    global feature_regex
    gene_name = label[0]
    feature_type = label[1]
    gene_counts[gene_name][feature_type].append(value)
    
def flatten_labels():
    global gene_counts
    full_labels = set()
    for gene in gene_counts.keys():
        for feature in gene_counts[gene].keys():
            full_labels.add((gene,feature))
    return full_labels
    
def drop_peaks(label,iv):
    
    global gene_counts
    
    gene_name = label[0]
    feature_type = label[1]
    
    if feature_type == 'gene':
        new_entry_list = []
        gene_label = (gene_name,"gene")
        i = 0
        for entry in get_label(gene_label):
            if entry[0] == iv:
                i+=1
            if i > 1:
                continue
            new_entry_list.append(entry)
        set_label(gene_label,new_entry_list)
    
    order = ['force','UTR','CDS','exon','transcript']
    position = [i for i, feature in enumerate(order) if feature in feature_type]
    if len(position)==0:
        return
    position = position[0]    
    to_del = set()
    dropped_from = set()
    for entry_type in gene_counts[gene_name].keys():
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
                    current_label = (gene_name,order[i])
                    for entry in get_label(current_label):
                        if entry[0] != iv:
                            new_entry_list.append(entry)
                    if len(new_entry_list) == 0:
                        to_del.add(current_label)
                    set_label(current_label,new_entry_list)
                    dropped_from.add(current_label)
            i+=1
    for key in to_del:
        delete_label(key)                   

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

def get_peaks(feature):
    return [iv_id_tuple for iv_id_tuple in gas[feature.iv].steps()]

def file_len(fname):
    i = 0
    if fname.endswith('.gz'):
        with gzip.open(fname) as f:
            for i, l in enumerate(f,1):
                pass
    else:
        with open(fname) as f:
            for i, l in enumerate(f,1):
                pass
    return i

def log_stats(tag_name):
    global checked_tags
    global unique_assigned_tags
    global multi_assigned
    global multiple_assigned_tags
    global checked
    global new_iv
    global clip_tag
    global counted
    
    if tag_name not in checked:
        print("{0}     {1}:   {2}:{3}-{4} not in checked!!".format(tag_name,new_iv,clip_tag.contig,clip_tag.start,clip_tag.end))
        import sys
        sys.exit()
    
    if tag_name not in checked_tags:
        unique_assigned_tags += 1
        checked_tags.add(tag_name)
        counted.add(tag_name)
        counted_dict[tag_name] = (new_iv,(clip_tag.contig,clip_tag.start,clip_tag.end))
    else:
        if tag_name not in multi_assigned:
            multiple_assigned_tags += 1
            multi_assigned.add(tag_name)

name_regex = re.compile(r'(.*)_[^_]*') 
feature_regex = re.compile(r'.*_([^_]*)')  

if no_peaks:
    peak_file = annotation_file

peaks_df = pd.read_table(peak_file, header=None, names=get_header(peak_file))
gas = HTSeq.GenomicArrayOfSets('auto',stranded=True)

peaks_df = peaks_df[~peaks_df['Region'].str.startswith('#')]
peaks_df['Region'].replace(regex="_",value="",inplace=True)

if annotation_file.endswith('gtf'):
    annotation_type = 'gtf'
  
#if no_peaks and annotation_type == 'gtf':
#   peaks_df['End'] = [x+1 for x in peaks_df['End']]

peak_io = io.StringIO()
peaks_df.to_csv(peak_io, sep="\t", line_terminator="\n", index=False, header=False, quoting=csv.QUOTE_NONE)
peak_io.seek(0)
    
end_included = True

if peak_file.endswith('.bed'):
    peaks = HTSeq.BED_Reader(peak_io)
else:
    peaks = HTSeq.GFF_Reader(peak_io,end_included=end_included)



for peak in peaks:
    # if peak.iv.chrom == 'chr10' and abs(peak.iv.start-97456891) < 10:
    #     print(peak.iv)
    #     import sys
    #     sys.exit()
    gas[peak.iv] += peak.iv.chrom + ':' + str(peak.iv.start) + '-' + str(peak.iv.end)

if no_peaks:
    annot_df = peaks_df
    peak_io.seek(0)
    annot_io = peak_io
else:
    annot_df = pd.read_table(annotation_file, header=None, names=get_header(annotation_file))
    annot_df['Region'].replace(regex="_",value="",inplace=True)
    annot_df = annot_df[~annot_df['Region'].str.startswith('#')]
    annot_io = io.StringIO()
    annot_df.to_csv(annot_io, sep="\t", line_terminator="\n", index=False, header=False, quoting=csv.QUOTE_NONE)
    annot_io.seek(0)
    
if annotation_type == 'gtf':
    annots = HTSeq.GFF_Reader(annot_io,end_included=end_included)
    
if annotation_file.endswith('.bed'):
    annots = HTSeq.BED_Reader(annot_io)
    annotation_type = 'bed'

#input_bam_reader = HTSeq.BAM_Reader(input_bam)

if stats:
    total_lines = file_len(clip_file)

if clip_type == 'bam':
    clip_reader = HTSeq.BAM_Reader(clip_file)
    bam_chroms = set()
    for alignment in HTSeq.BAM_Reader(clip_file):
        chrom = alignment.iv.chrom.replace("_","")
        if chrom not in bam_chroms:
            bam_chroms.add(chrom)
elif clip_type == 'bed':
    clip_reader = pysam.TabixFile(clip_file)
    
gene_counts = {}

l = len(pd.read_table(annotation_file,header=None).index)
checked_features = set()
if no_peaks:
    chroms = set(annot_df['Region'].values)
else:
    chroms = set(peaks_df['Region'].values)
checked_region_ivs = {}
checked_region_reads = {}

if strict or stats:
    checked_tags = set()
    if stats:
        unique_assigned_tags = 0
        multiple_assigned_tags = 0
        multi_assigned = set()

if stats:
    total_tags = int(total_lines if clip_type == 'bed' else total_lines/4 if clip_type == 'bam' else 0)
    k = 0
    checked = set()
    checked_dict = dict()
    
    if no_peaks and annotation_type == 'gtf':
        if end_included:
            end_increment = 0
        else:
            end_increment = 1
    
    for idx,peak in peaks_df.iterrows():
        if clip_type == 'bed' and peak['Region'] not in clip_reader.contigs:
            continue
        if clip_type == 'bam' and peak['Region'] not in bam_chroms:
            continue
        start = int(peak['Start']) - 1
        end = int(peak['End']) + end_increment
        if clip_type == 'bed':
            tags = clip_reader.fetch(peak['Region'],start,end,parser=pysam.asBed())
        if clip_type == 'bam':
            tags = clip_reader.fetch(peak['Region'],start,end)
        for tag in tags:
        #     if tag.name == 'NS500289:958:HGFMTBGXJ:1:22302:4489:15069#3#AACCGCCCTTGTATAC':
        #         print(peak['Region'],int(peak['Start'])-1,int(peak['End']))
        #         print(HTSeq.GenomicInterval(peak['Region'],int(peak['Start'])-1,int(peak['End'])))
        #         print((tag.contig,tag.start,tag.end),HTSeq.GenomicInterval(tag.contig,tag.start,tag.end))
        #         print(HTSeq.GenomicInterval(peak['Region'],int(peak['Start'])-1,int(peak['End'])).overlaps(HTSeq.GenomicInterval(tag.contig,tag.start,tag.end)))
        #         print(HTSeq.GenomicInterval(peak['Region'],int(peak['Start'])-1,int(peak['End'])),HTSeq.GenomicInterval(tag.contig,tag.start,tag.end))
        #         print(HTSeq.GenomicInterval(peak['Region'],int(peak['Start'])-1,int(peak['End'])-1).overlaps(HTSeq.GenomicInterval(tag.contig,tag.start,tag.end)))
        #         print(HTSeq.GenomicInterval(peak['Region'],int(peak['Start'])-1,int(peak['End'])-1),HTSeq.GenomicInterval(tag.contig,tag.start,tag.end))
        #         import sys
        #         sys.exit()
        #     if peak['Region']=='chr11':
        #         import sys
        #         sys.exit()
            if tag.name not in checked:
                k+=1
                checked.add(tag.name)
                checked_dict[tag.name] = (HTSeq.GenomicInterval(peak['Region'],int(peak['Start']-1),int(peak['End']),peak['Strand']),(tag.contig,tag.start,tag.end))
    total_peak_tags = k

counted = set()
counted_dict = dict()

new_chrom = None
current_chrom = None
i = 0
for feature in annots:
    i+=1
    new_chrom = feature.iv.chrom
    if new_chrom != current_chrom:
        chrom_labels = set()
        checked_region_ivs[new_chrom] = dict()
        checked_region_reads[new_chrom] = dict()
        current_chrom = new_chrom
        
    if feature.name+'_'+current_chrom+':'+str(feature.iv.start)+'-'+str(feature.iv.end) not in checked_features:
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
        label = (feature.name,feature.type)
        gene_label = (feature.name,'gene')
    if annotation_type == 'bed':
        label = (feature.name,None)
    if label in chrom_labels:
        checked_ivs = checked_region_ivs[current_chrom][label]
        clip_read_names = checked_region_reads[current_chrom][label]
    
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
        if clip_type == 'bam':
            clip_tags = clip_reader.fetch(chrom,start,end)
        if clip_type == 'bed' and chrom.replace('_','') in clip_reader.contigs:
            clip_tags = clip_reader.fetch(chrom,start,end,parser=pysam.asBed())
        else:
            if no_peaks:
                continue
            print('Error! Tag file missing contig! Please verify file integrity and try again.')
            continue
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
                tag_strand = clip_tag.iv.strand
            elif clip_type == 'bed':
                tag_name = clip_tag.name
                tag_strand = clip_tag.strand
                
            if strict and tag_name in checked_tags:
                    continue
            
            if tag_name in clip_read_names:
                continue
            if clip_type == 'bam' and clip_tag.aQual < mapq:
                continue
            if clip_type == 'bam' and clip_tag.aligned:
                if tag_strand == new_iv.strand:
                    increment = process_cigar(clip_tag,start,end)
                    clip_count += increment
                    clip_read_names.add(tag_name)
                    if stats and increment > 0:
                        log_stats(tag_name)
                else:
                    continue
            elif clip_type == 'bed':
                if True:#tag_strand == new_iv.strand:
                    clip_count += 1
                    clip_read_names.add(tag_name)
                    if stats:
                        log_stats(tag_name)
                else:
                    continue
            
            if out_bed:
                if clip_type == 'bam':
                    row = (clip_tag.iv.chrom,clip_tag.iv.start,clip_tag.iv.end,clip_tag.name,'.',clip_tag.iv.strand)
                if clip_type == 'bed':
                    row = (clip_tag.contig,clip_tag.start,clip_tag.end,clip_tag.name,clip_tag.score,clip_tag.strand)
                bed_list.append(row)
        
        if clip_count == 0 and not keep_zeros:
            continue
        
        if label not in chrom_labels:
            chrom_labels.add(label)
            checked_region_ivs[current_chrom][label] = checked_ivs
            checked_region_reads[current_chrom][label] = clip_read_names
        else:
            checked_region_ivs[current_chrom][label] = set.union(checked_region_ivs[current_chrom][label],checked_ivs)
            checked_region_reads[current_chrom][label] = set.union(checked_region_reads[current_chrom][label],clip_read_names)
            
        if not check_label(label):
            set_label(label,[(HTSeq.GenomicInterval(chrom,start,end),clip_count)]) #[clip_count-input_count] #[clip_count/input_count] for FC
        else:
            append_label(label,(HTSeq.GenomicInterval(chrom,start,end),clip_count)) #.append(clip_count-input_count) #.append(clip_count/input_count)
        
        if annotation_type == 'gtf':
            if not check_label(gene_label):
                set_label(gene_label,[(HTSeq.GenomicInterval(chrom,start,end),clip_count)])
            else:
                append_label(gene_label,(HTSeq.GenomicInterval(chrom,start,end),clip_count))
            drop_peaks(label,(HTSeq.GenomicInterval(chrom,start,end)))
            drop_peaks(gene_label,(HTSeq.GenomicInterval(chrom,start,end)))
            

diff = checked - counted
for tag in diff:
    print(checked_dict[tag])
# import sys
# sys.exit()


final_gene_count = {}

gene_list = flatten_labels()

if annotation_type == 'gtf':
    gene_string_list = list(map(lambda x:"{0}_{1}".format(x[0],x[1]),gene_list))
if annotation_type == 'bed':
    gene_string_list = list(map(lambda x:"{0}".format(x[0]),gene_list))

if no_peaks:
    peak_label = "Feature_Count"
    peak_tag_label = "Feature_Tags"
else:
    peak_label = "Peak_Count"
    peak_tag_label = "Peak_Tags"
    
final_gene_count_df = pd.DataFrame()
final_gene_count_df['Annotation'] = gene_string_list
final_gene_count_df[peak_label] = 0
final_gene_count_df['Average_Tags'] = 0
final_gene_count_df['Total_Tags'] = 0
final_gene_count_df[peak_tag_label] = ""

for gene_label in gene_list:
    if annotation_type == 'gtf':
        label_string = "{0}_{1}".format(gene_label[0],gene_label[1])
    if annotation_type == 'bed':
        label_string = "{0}".format(gene_label[0])
    gene_count_array = [tup[1] for tup in get_label(gene_label)]
    if len(gene_count_array) > 0:
        matches = (final_gene_count_df['Annotation']==label_string)
        final_gene_count_df.loc[matches,peak_label] = len(gene_count_array)
        final_gene_count_df.loc[matches,'Average_Tags'] = np.average(gene_count_array) #FC
        final_gene_count_df.loc[matches,'Total_Tags'] = np.sum(gene_count_array)
        final_gene_count_df.loc[matches,peak_tag_label] = '['+', '.join(map(str,[(str(tup[0]),tup[1]) for tup in get_label(gene_label)]))+']'

if out_bed:
    bed = pd.DataFrame(bed_list)
    f = open(bed_output, 'w')
    bed.to_csv(f, sep="\t", line_terminator="\n", index=False, header=False)
    f.close()

if stats:
    stats_series = pd.Series({"Total_Library_Tags:":total_tags,"Total_Assessed_Tags:":total_peak_tags,"Unique_Assigned_Tags:":unique_assigned_tags,"Multiple_Assigned_Tags:":multiple_assigned_tags,"Total_Assigned_Tags:":final_gene_count_df['Total_Tags'].sum()})
    f = open(out_file+'.stats', 'w')
    stats_series.to_csv(f, sep=" ", line_terminator="\n", index=True, header=False)
    f.close()
    
    stats_series.to_csv
    
final_gene_count_df.sort_values(by='Total_Tags', axis=0, ascending=False, inplace=True)
f = open(out_file, 'w')
final_gene_count_df.to_csv(f, sep="\t", line_terminator="\n", index=False, header=True)
f.close()