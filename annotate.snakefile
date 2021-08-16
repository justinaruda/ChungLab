#samples = ["A1KO_IFN_Down.pool.tag.uniq.peak.sig.1001nt.unedited_vs_A1KO_IFN_Down.pool.tag.uniq.peak.sig.1001nt.edited","A1KO_IFN_Up.pool.tag.uniq.peak.sig.1001nt.unedited_vs_A1KO_IFN_Up.pool.tag.uniq.peak.sig.1001nt.edited","A1KO_Mock_Down.pool.tag.uniq.peak.sig.1001nt.unedited_vs_A1KO_Mock_Down.pool.tag.uniq.peak.sig.1001nt.edited","A1KO_Mock_Up.pool.tag.uniq.peak.sig.1001nt.unedited_vs_A1KO_Mock_Up.pool.tag.uniq.peak.sig.1001nt.edited","random_intervals.1001nt.unedited_vs_random_intervals.1001nt.edited"]
#samples = ["13T.pool.tag.uniq.CITS.s30.singleton.1001nt.unedited_vs_13T.pool.tag.uniq.CITS.s30.singleton.1001nt.edited","14T.pool.tag.uniq.CITS.s30.singleton.1001nt.unedited_vs_14T.pool.tag.uniq.CITS.s30.singleton.1001nt.edited","A1KO_IFN_Down.pool.tag.uniq.CITS.s30.singleton.1001nt.unedited_vs_A1KO_IFN_Down.pool.tag.uniq.CITS.s30.singleton.1001nt.edited","A1KO_IFN_Up.pool.tag.uniq.CITS.s30.singleton.1001nt.unedited_vs_A1KO_IFN_Up.pool.tag.uniq.CITS.s30.singleton.1001nt.edited","A1KO_Mock_Down.pool.tag.uniq.CITS.s30.singleton.1001nt.unedited_vs_A1KO_Mock_Down.pool.tag.uniq.CITS.s30.singleton.1001nt.edited","A1KO_Mock_Up.pool.tag.uniq.CITS.s30.singleton.1001nt.unedited_vs_A1KO_Mock_Up.pool.tag.uniq.CITS.s30.singleton.1001nt.edited"]
samples = glob_wildcards("{sample}_compared.out")

rule all:
	input:
		expand("{sample}_compared.ir.out",sample = samples)

rule annotate:
	input:
		"{sample}_compared.out"
	output:
		"{sample}_compared.ir.out"
	resources:
		memory=15
	shell:
		"python annotate_ir.py -t {input} -a /home/hchunglab/data/justin/RNA*/genome/pri*/repeat_masker_grch38.txt -o {output}"