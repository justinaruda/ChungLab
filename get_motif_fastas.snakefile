samples = ["A1KO_IFN_Down","A1KO_IFN_Up","A1KO_Mock_Down","A1KO_Mock_Up"]
mutations = ["del"]
types = ["CITS.s30.singleton"] #"CIMS.s30"
ivs = [101]

rule all:
	input:
		expand("fastas/{sample}.pool.tag.uniq.{mutation}.{type}.{iv}nt.slop.fasta",sample=samples,mutation=mutations,type=types,iv=ivs)

rule get_intervals:
	input:
		"../{sample}.pool.tag.uniq.{mutation}.{type}.bed"
	output:
		"bed_intervals/{sample}.pool.tag.uniq.{mutation}.{type}.{iv}nt.slop.bed"
	params:
		half_window=lambda wildcards:round((int(wildcards.iv)-1)/2),
	shell:
		"bedtools slop -b {params.half_window} -i {input} -g /home/hchunglab/data/justin/RNA*/genome/prim*/grch38_primary_assembly_measles.genome > {output}"""

rule getfasta:
	input:
		rules.get_intervals.output
	output:
		"fastas/{sample}.pool.tag.uniq.{mutation}.{type}.{iv}nt.slop.fasta"
	shell:
		"bedtools getfasta -fi /home/hchunglab/data/justin/RNA*/genome/pri*/grch38_primary_assembly_measles.fna -s -bed {input} -fo {output}"