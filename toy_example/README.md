# Snakemake toy example

This is a simple toy example that can be used to start learning the basics of snakemake. This folder contains some extremely simple input files. For example, sampleA.txt contains three lines of data, with values 1, 2, and 3, respectively. We'll use these simple data along with core command-line tools, included by default in most unix-like environments, to illustrate the basics of snakemake.

Note: If you haven't installed snakemake yet, I'd suggest using [conda](https://docs.conda.io/en/latest/miniconda.html) to do so.

Topics covered:
1. Setting up workflow targets and dependencies
2. The basics of writing rules
3. The workflow & DAG

Get started by expanding Example 1a

<details><summary>Expand - Ex. 1a</summary>

ex-1a.smk has the following contents:

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


Review the contents of the file, and consider the following topics:
* Targets & dependencies
* Writing rules
* Intro to wildcards


To perform a dry-run:

    snakemake --snakefile ex-1a.smk --dry-run

If all looks well, run example 1a (remove `--dry-run` from the previous command)

    snakemake --snakefile ex-1a.smk

Review the ouputs. Were they generated successfully? What happens if we try running it again?

    snakemake --snakefile ex-1a.smk

Try deleting the output directory. What happens if we run it then?

    rm -r output/
    snakemake --snakefile ex-1a.smk


## You have reached the end of example 1a ✅

</details>


Now we'll make the workflow a bit more interesting. We'll add more rules, use a configuration file, and more!

<details><summary>Expand - Ex. 1b</summary>


ex-1b.smk contents:

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
        "output/{base}_appended.txt"
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
        "sleep 2 ; sort -n {input} > {output}"

config-ex1.yml contents:

    basenames:
      - 'sampleA'
      - 'sampleB'
      - 'sampleC'
    append_val: 42

<details><summary>Ex. 1b rulegraph</summary>

![Ex. 1b rulegraph](https://github.com/um-dang/intro-to-snakemake/blob/master/img/rg-ex-1b.pdf)

</details>

Review these files and consider the following topics:
* The DAG/rulegraph
* The configuration file
* The expand statement



To dry-run the new snakefile/configfile

    snakemake --snakefile ex-1b.smk --configfile config-ex1.yml --dry-run

If all looks well, run the workflow (remove the `--dry-run` flag)

    snakemake --snakefile ex-1b.smk --configfile config-ex1.yml

Review the directory of results. Did snakemake run the workflow, and successfully create the desired targets?

Try deleting an intermediate file, then running the pipeline again. What happens?

    rm output/sampleA_appended.txt
    snakemake --snakefile ex-1b.smk --configfile config-ex1.yml

What about modifying an intermediate file, and running the pipeline again?

    echo 101 >> output/sampleB_appended.txt
    snakemake --snakefile ex-1b.smk --configfile config-ex1.yml


Try creating the DAG or rulegraph

    #DAG (file-level granularity)
    snakemake --snakefile ex-1b.smk --configfile config-ex1.yml --dag | dot -T pdf > dag-ex-1b.pdf
    #Rulegraph (rule-level granularity)
    snakemake --snakefile ex-1b.smk --configfile config-ex1.yml --rulegraph | dot -T pdf > rg-ex-1b.pdf

## You have reached the end of example 1b ✅

</details>
