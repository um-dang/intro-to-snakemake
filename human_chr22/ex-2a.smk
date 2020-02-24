rule all:
    input:
        expand("aligned/{sample}.sam", sample=config['samples'])

# rule fasta_gunzip:
#     input:
#
#     output:
#
#     shell:
#

rule bwa_index:
    input:
        config['ref_prefix'] + ".fa"
    output:
        config['ref_prefix'] + ".amb",
        config['ref_prefix'] + ".ann",
        config['ref_prefix'] + ".bwt",
        config['ref_prefix'] + ".pac",
        config['ref_prefix'] + ".sa"
    params:
    	ref_prefix = config['ref_prefix']
    shell:
        "bwa index -p {params.ref_prefix} {input}"

rule bwa_mem:
    input:
        sa = "reference/chr22.sa",
        r1 = "read_data/{sample}/{sample}_R1.fastq.gz",
        r2 = "read_data/{sample}/{sample}_R2.fastq.gz"
    output:
        "aligned/{sample}.sam"
    threads: 2
    params:
        ref_prefix = config['ref_prefix']
    shell:
        "bwa mem -t {threads} {params.ref_prefix} {input.r1} {input.r2} > {output}"
