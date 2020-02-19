rule all:
    input:
        "output/123_rsorted.txt",
        "output/345_rsorted.txt",
        "output/567_rsorted.txt"

rule rsort:
    input:
        "{basename}.txt"
    output:
        "output/{basename}_rsorted.txt"
    shell:
        "sort -r {input} > {output}"
