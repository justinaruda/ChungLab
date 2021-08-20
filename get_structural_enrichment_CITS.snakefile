samples = ["13T","14T","A1KO_Mock_Down","A1KO_Mock_Up","A1KO_IFN_Down","A1KO_IFN_Up"]
ivs = [1001]
mutations = ["del."]
genomes = ["edited","unedited"]

def get_genome(wildcards):
	if wildcards.genome == "unedited":
		return "grch38_primary_assembly_measles"
	if wildcards.genome == "edited":
		return "grch38_primary_assembly_measles_REDIportal_fixed_contigs"

def get_random_genome(wildcards):
	if wildcards.genome == "unedited":
		return "gencode.v38.expressed"
	if wildcards.genome == "edited":
		return "gencode.v38.expressed.edited"

def get_inputs(wildcards):
	comparisons, = glob_wildcards("aggregated/{comparison}_compared.out")
  	return expand("aggregated/{comparison}_compared.out",comparison=comparisons)

rule get_singleton:
	input:
		"../{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.bed"
	output:
		"../{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.bed"
	shell:
		"""awk '{{if($3-$2==1) {{print $0}}}}' {input} > {output}"""

rule get_intervals:
	input:
		"../{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.bed"
	output:
		"bed_intervals/{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.{iv}nt.{genome}.slop.bed"
	params:
		half_window=lambda wildcards:round((int(wildcards.iv)-1)/2),
		genome=get_genome
	shell:
		"bedtools slop -b {params.half_window} -i {input} -g /home/hchunglab/data/justin/RNA*/genome/prim*/{params.genome}.genome > {output}"""

rule get_random_interval:
	output:
		"bed_intervals/random_intervals.{iv}nt.bed"
	shell:
		"bedtools random -l {wildcards.iv} -n 3000 -g /home/hchunglab/data/justin/RNA*/genome/prim*/gencode.v38.expressed.genome > {output}"

rule get_fastas:
	input:
		"bed_intervals/{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.{iv}nt.{genome}.slop.bed"
	output:
		"fastas/{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.{iv}nt.{genome}.slop.fasta"
	params:
		genome=get_genome
	shell:
		"bedtools getfasta -fi /home/hchunglab/data/justin/RNA*/genome/pri*/{params.genome}.fna -s -bed {input} > {output}"

rule get_random_fastas:
	input:
		"bed_intervals/random_intervals.{iv}nt.bed"
	output:
		"fastas/random_intervals.{iv}nt.{genome}.fasta"
	params:
		genome=get_random_genome
	shell:
		"bedtools getfasta -fi /home/hchunglab/data/justin/RNA*/genome/pri*/{params.genome}.fna -s -bed {input} > {output}"

rule fold:
	input:
		"fastas/{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.{iv}nt.{genome}.slop.fasta"
	output:
		"RNAfold/{sample}.pool.tag.uniq.{mutation,.*}CITS.s30.singleton.{iv}nt.{genome}.fold"
	threads: workflow.cores
	shell:
		"RNAfold --jobs={threads} < {input} > {output} && rm *.ps"

rule random_fold:
	input:
		"fastas/random_intervals.{iv}nt.{genome}.fasta"
	output:
		"RNAfold/random_intervals.{iv}nt.{genome}.fold"
	threads: workflow.cores
	shell:
		"RNAfold --jobs={threads} < {input} > {output} && rm *.ps"

def combine_inputs(wildcards):
	a = expand("RNAfold/{sample}.pool.tag.uniq.{mutation}CITS.s30.singleton.{iv}nt.{genome}.fold",sample=samples,iv=ivs,genome=genomes,mutation=mutations)
	b = expand("RNAfold/random_intervals.{iv}nt.{genome}.fold",iv=ivs,genome=genomes)
	return a + b

rule aggregate_ivs:
	input:
		combine_inputs
	output:
		"aggregated/mfe.summary.txt"
	script:
		"aggregate.py"

rule annotate:
        input:
                aggregate = rules.aggregate_ivs.output,
                files = get_inputs
        resources:
                memory=13
        params:
                annotation="/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/repeat_masker_grch38.txt",
                genome_annotation="/home/hchunglab/data/justin/RNAseq/genome/primary_assembly/gencode.v38.primary_assembly.annotation.gene_name.gtf"
        script:
                "annotate_ir.py"