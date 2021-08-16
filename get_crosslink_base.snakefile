samples = ["A1KO_IFN_Down","A1KO_IFN_Up","A1KO_Mock_Down","A1KO_Mock_Up"]
mutations = ["ins.","sub.","del."]
types = ["CIMS.s30"] #"CIMS.s30"

rule all:
	input:
		expand("bases/{sample}.pool.tag.uniq.{mutation}{type}.quant",sample=samples,mutation=mutations,type=types)

rule getfasta:
	input:
		"../{sample}.pool.tag.uniq.{mutation}{type}.bed"
	output:
		"bases/{sample}.pool.tag.uniq.{mutation}{type}.fasta"
	shell:
		"bedtools getfasta -fi /home/hchunglab/data/justin/RNA*/genome/pri*/grch38_primary_assembly_measles.fna -s -rna -bed {input} -fo {output}"

rule quantify:
	input:
		"bases/{sample}.pool.tag.uniq.{mutation}{type}.fasta"
	output:
		"bases/{sample}.pool.tag.uniq.{mutation}{type}.quant"
	shell:
		"""awk 'BEGIN {{a=0;u=0;c=0;g=0;other=0}} $1~/^>/ {{next}} $1~/A/ {{a=a+1;next}} $1~/(U|T)/ {{u=u+1;next}} $1~/G/ {{g=g+1;next}} $1~/C/ {{c=c+1;next}} {{other=other+1}} END {{total=a+c+u+g+other;print "A: "a" ("a/total*100"%)","\\nU: "u" ("u/total*100"%)","\\nC: "c" ("c/total*100"%)","\\nG: "g" ("g/total*100"%)","\\nOther: "other" ("other/total*100"%)"}}' {input} > {output}"""