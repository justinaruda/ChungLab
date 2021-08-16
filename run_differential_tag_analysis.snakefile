samples = ["A1KO_Mock_Down","A1KO_Mock_Up","A1KO_IFN_Down","A1KO_IFN_Up"]
comparisons = ["A1KO_Mock_Down_vs_A1KO_IFN_Down","A1KO_Mock_Up_vs_A1KO_IFN_Up","A1KO_Mock_Down_vs_A1KO_Mock_Up","A1KO_IFN_Down_vs_A1KO_IFN_Up"]
tags = ["all_tags","peak_tags"]
annotations = ['repeatMasker','ncbiRefSeq'] 

def get_annotations(wildcards):
	if wildcards.annotation == "repeatMasker":
		return "/home/hchunglab/data/justin/RNA*/genome/pri*/repeat_masker_grch38.fix.bed"
	if wildcards.annotation == "ncbiRefSeq":
		return "/home/hchunglab/data/justin/RNA*/genome/pri*/hg38.ncbiRefSeq.measles.gtf"

def get_peaks(wildcards):
	if wildcards.tags == "all_tags":
		return ""
	if wildcards.tags == "peak_tags":
		return "-p /home/hchunglab/data/heeg*/B*/N*/20*/cluster/{0}.pool.tag.uniq.peak.sig.halfPH.bed".format(wildcards.sample)

def get_memory(wildcards):
	if wildcards.tags == "all_tags":
		return 27
	if wildcards.tags == "peak_tags":
		return 15

rule dedup:
	input:
		"/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/parsing/merge2/{sample}.pool.tag.uniq.rgb.bed"
	output:
		"/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/parsing/merge2/{sample}.pool.tag.uniq.rgb.dedup.bed"
	shell:
		"python /home/hchunglab/tmp/remove_duplicates.py -b {input} -o {output}"

rule sort:
	input:
		rules.dedup.output
	output:
		temp("/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/parsing/merge2/snakemake/{sample}.pool.tag.uniq.rgb.dedup.sorted.bed")
	shell:
		"sort -V -k1,1 -k2,2 {input} > {output}"

rule bgzip:
	input:
		rules.sort.output
	output:
		"/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/parsing/merge2/snakemake/{sample}.pool.tag.uniq.rgb.dedup.sorted.bed.gz"
	shell:
		"bgzip {input}"

rule tabix:
	input:
		rules.bgzip.output
	output:
		"/home/hchunglab/data/heegwon_shin/Basespace/Nextseq/20210510/parsing/merge2/snakemake/{sample}.pool.tag.uniq.rgb.dedup.sorted.bed.gz.tbi"
	shell:
		"tabix -p bed {input}"

rule quantify:
	input:
		bgzip = rules.bgzip.output,
		tabix = rules.tabix.output
	output:
		quant = "{annotation}/{tags}/{sample}.tags.quantified.txt",
		bed = "{annotation}/{tags}/{sample}.tags.quantified.bed"
	params:
		annotation = get_annotations,
		peaks = get_peaks
	resources:
		memory = get_memory
	shell:
		"python /home/hchunglab/tmp/quantify_hitsclip_bed.py -b {input.bgzip} -a {params.annotation} {params.peaks} -o {output.quant} -x {output.bed}"

rule aggregate:
	input:
		expand("{annotation}/{tags}/{sample}.tags.quantified.txt",annotation=annotations,tags=tags,sample=samples)
	output:
		expand("{annotation}/{tags}/comparisons/{comparison}.{annotation}.{tags}.txt",annotation=annotations,tags=tags,comparison=comparisons)
	script:
		"differential_analysis.py"