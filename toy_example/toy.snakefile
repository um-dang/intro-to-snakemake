
rule all:
	input:
		expand("output/{bname}_fsorted.txt", bname=config['basenames']),
		expand("output/{bname}_randsorted.txt", bname=config['basenames'])

rule rsort:
	input:
		"{base}.txt"
	output:
		"output/{base}_rsorted.txt"
	shell:
		"sort -r {input} > {output}"

rule append_value:
	input:
		"output/{base}_rsorted.txt"
	output:
		"output/{base}_appended.txt"
	params:
		append_val = config['append_val']
	shell:
		"cat {input} > {output} ; "
		"echo {params.append_val} >> {output}"

rule randsort:
	input:
		"{base}.txt"
	output:
		"output/{base}_randsorted.txt"
	shell:
		"sort -R {input} > {output}"

rule fsort:
	input:
		"output/{base}_appended.txt"
	output:
		"output/{base}_fsorted.txt"
	shell:
		"sleep 2 ; sort {input} > {output}"
