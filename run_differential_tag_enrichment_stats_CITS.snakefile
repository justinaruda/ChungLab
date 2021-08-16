samples = ["A1KO_Mock_Down","A1KO_Mock_Up","A1KO_IFN_Down","A1KO_IFN_Up"]
comparisons = ["A1KO_Mock_Down_vs_A1KO_IFN_Down","A1KO_Mock_Up_vs_A1KO_IFN_Up","A1KO_Mock_Down_vs_A1KO_Mock_Up","A1KO_IFN_Down_vs_A1KO_IFN_Up"]
tags = ["all_tags","peak_tags"]
mutations = ["del."]
annotations = ["gencode"] 

def get_annotations(wildcards):
	if wildcards.annotation == "repeatMasker":
		return "/home/hchunglab/data/justin/RNA*/genome/pri*/repeat_masker_grch38.fix.bed"
	if wildcards.annotation == "ncbiRefSeq":
		return "/home/hchunglab/data/justin/RNA*/genome/pri*/hg38.ncbiRefSeq.measles.gtf"
	if wildcards.annotation == "gencode":
		return "/home/hchunglab/data/justin/RNA*/genome/pri*/gencode.v38.primary_assembly.annotation.gene_name.gtf"

def get_peaks(wildcards):
	if wildcards.tags == "all_tags":
		return ""
	if wildcards.tags == "peak_tags":
		return "-p /home/hchunglab/data/heeg*/B*/N*/20*/cluster/{0}.pool.tag.uniq.peak.sig.halfPH.bed".format(wildcards.sample)

def get_memory(wildcards):
	if wildcards.tags == "all_tags":
		return 14
	if wildcards.tags == "peak_tags":
		return 8

rule dedup:
	input:
		"../{sample}.pool.tag.uniq.{mutation}CITS.s30.singleton.bed"
	output:
		"files/{sample}.pool.tag.uniq.{mutation}CITS.s30.singleton.dedup.bed"
	shell:
		"python /home/hchunglab/tmp/remove_duplicates.py -b {input} -o {output}"

rule sort:
	input:
		rules.dedup.output
	output:
		temp("files/{sample}.pool.tag.uniq.{mutation}CITS.s30.singleton.bed")
	shell:
		"sort -V -k1,1 -k2,2 {input} > {output}"

rule bgzip:
	input:
		rules.sort.output
	output:
		"files/{sample}.pool.tag.uniq.{mutation}CITS.s30.singleton.bed.gz"
	shell:
		"bgzip {input}"

rule tabix:
	input:
		rules.bgzip.output
	output:
		"files/{sample}.pool.tag.uniq.{mutation}CITS.s30.singleton.bed.gz.tbi"
	shell:
		"tabix -p bed {input}"

rule quantify:
	input:
		bgzip = rules.bgzip.output,
		tabix = rules.tabix.output
	output:
		quant = "{annotation}/{tags}/{sample}.{mutation}CITS.quantified.txt",
		bed = "{annotation}/{tags}/{sample}.{mutation}CITS.quantified.bed"
	params:
		annotation = get_annotations,
		peaks = get_peaks
	resources:
		memory = get_memory
	shell:
		"python /home/hchunglab/tmp/quantify_hitsclip_stats.py -b {input.bgzip} -a {params.annotation} {params.peaks} -o {output.quant} -x {output.bed} -l"

rule aggregate:
	input:
		expand("{annotation}/{tags}/{sample}.{mutation}CITS.quantified.txt",annotation=annotations,tags=tags,sample=samples,mutation=mutations)
	output:
		expand("{annotation}/{tags}/comparisons/{comparison}.{annotation}.{tags}.txt",annotation=annotations,tags=tags,comparison=comparisons,mutation=mutations)
	params:
		summary="aggregated_sample_summary.txt"
	threads: workflow.cores
	script:
		"differential_analysis_stats.py"