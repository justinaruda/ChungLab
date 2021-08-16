samples = ["A1KO_Mock_Down","A1KO_Mock_Up","A1KO_IFN_Down","A1KO_IFN_Up"]
mutations = ["sub"]
ivs = [1001]
genomes = ["edited","unedited"]

def get_genome(wildcards):
	if wildcards.genome == "unedited":
		return "grch38_primary_assembly_measles"
	if wildcards.genome == "edited":
		return "grch38_primary_assembly_measles_REDIportal_fixed_contigs"

rule get_intervals:
	input:
		"../{sample}.pool.tag.uniq.{mutation}.CIMS.s30.bed"
	output:
		"bed_intervals/{sample}.pool.tag.uniq.{mutation}.CIMS.s30.{iv}nt.slop.bed"
	params:
		half_window=lambda wildcards:round((int(wildcards.iv)-1)/2),
		genome=get_genome
	shell:
		"bedtools slop -b {params.half_window} -i {input} -g /home/hchunglab/data/justin/RNA*/genome/prim*/{params.genome}.genome > {output}"""

rule get_random_interval:
	output:
		"bed_intervals/random_intervals.{iv}nt.bed"
	shell:
		"bedtools random -l {wildcards.iv} -n 1500 -g /home/hchunglab/data/justin/RNA*/genome/prim*/grch38_primary_assembly_measles.genome > {output}"

rule get_fastas:
	input:
		"bed_intervals/{sample}.pool.tag.uniq.{mutation}.CIMS.s30.{iv}nt.slop.bed"
	output:
		"fastas/{sample}.pool.tag.uniq.{mutation}.CIMS.s30.{iv}nt.{genome}.slop.fasta"
	params:
		genome=get_genome
	shell:
		"bedtools getfasta -fi /home/hchunglab/data/justin/RNA*/genome/pri*/{params.genome}.fna -bed {input} > {output}"

rule get_random_fastas:
	input:
		"bed_intervals/random_intervals.{iv}nt.bed"
	output:
		"fastas/random_intervals.{iv}nt.{genome}.fasta"
	params:
		genome=get_genome
	shell:
		"bedtools getfasta -fi /home/hchunglab/data/justin/RNA*/genome/pri*/{params.genome}.fna -bed {input} > {output}"

rule fold:
	input:
		"fastas/{sample}.pool.tag.uniq.{mutation}.CIMS.s30.{iv}nt.{genome}.slop.fasta"
	output:
		"RNAfold/{sample}.pool.tag.uniq.{mutation}.CIMS.s30.{iv}nt.{genome}.noGU.fold"
	threads: 8
	shell:
		"RNAfold --jobs=8 -noGU < {input} > {output} && rm *.ps"

rule random_fold:
	input:
		"fastas/random_intervals.{iv}nt.{genome}.fasta"
	output:
		"RNAfold/random_intervals.{iv}nt.{genome}.noGU.fold"
	threads:8
	shell:
		"RNAfold --jobs=8 -noGU < {input} > {output} && rm *.ps"

def combine_inputs(wildcards):
	a = expand("RNAfold/{sample}.pool.tag.uniq.{mutation}.CIMS.s30.{iv}nt.{genome}.fold",sample=samples,mutation=mutations,iv=ivs,genome=genomes)
	b = expand("RNAfold/random_intervals.{iv}nt.{genome}.fold",iv=ivs,genome=genomes)
	return a + b

rule aggregate_ivs:
	input:
		combine_inputs
	output:
		"aggregated/mfe.summary.txt"
	script:
		"aggregate.py"