rule all:
    input:
        "output/sampleA_rsorted.txt",
        "output/sampleB_rsorted.txt",
        "output/sampleC_rsorted.txt"
        
rule rsort:
    input:
        "{basename}.txt"
    output:
        "output/{basename}_rsorted.txt"
    shell:
        "sort -r {input} > {output}"
