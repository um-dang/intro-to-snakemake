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

<details><summary>Ex. 1a rulegraph</summary>

![Ex. 1a rulegraph](https://github.com/um-dang/intro-to-snakemake/blob/master/img/rg-ex-1a.pdf)

</details>

Topics covered:
* Targets & dependencies
* Writing rules

Targets - Have to think about the pipeline backwards - What do we want to end up with?

Dependencies - What rules must be in place (inputs/outputs) for these targets to be generated?
 - This will be more obvious with the addition of more rules

Writing rules - Generally will have 'input', 'output', and 'shell' blocks (more to the story)

The rule 'all' is placed at the top of the file (the first rule, anyway), and this is always executed by default. It's being used to define the targets for the workflow.

To perform a dry-run:

    snakemake --snakefile ex-1a.smk --dry-run

Notice that snakemake keeps track of the wildcards during the evaluation of each rule
* experiment by changing the targets so they don't match the input files

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

config-ex-1b.yml contents:

    basenames:
      - 'sampleA'
      - 'sampleB'
      - 'sampleC'
    append_val: 42

<details><summary>Ex. 1b rulegraph</summary>
![Ex. 1b rulegraph](../img/rg-ex-1b.pdf)
</details>

The config
* How is the config used with this snakefile?

The expand statement
* The various uses of curly braces can be confusing at first (at least for me)
* `expand` is distinct from `wildcards`
* can be thought of as "expand this string (arg 1) into an array of strings, filling in all combinations of values (args 2+ as key-value pairs)

Running the new snakefile/configfile

    snakemake --snakefile ex-1b.smk --configfile config-ex-1b.yml

Did snakemake run the workflow, and successfully create the desired targets?
* View the directory of results

More about the core tenets of snakemake (also gnu make, make-like things)
* Try running the workflow to completion, then running it again. What happens?
* Delete `output`, then try again. Isn't this cool?
* Try modifying an intermediate file, then running the pipeline again. How is this beneficial? How can it be problematic?

The workflow & DAG
* Directed Acyclic Graph - how snakemake 'knows' how to produce the desired targets
* It can be useful to see the workflow DAG, and imagine how snakemake 'thinks' about executing it

Viewing the DAG (or rulegraph)

    #DAG (file-level granularity)
    snakemake --snakefile ex-1b.smk --configfile config-ex-1b.yml --dag | dot -T pdf > dag-ex-1b.pdf
    #Rulegraph (rule-level granularity)
    snakemake --snakefile ex-1b.smk --configfile config-ex-1b.yml --rulegraph | dot -T pdf > rg-ex-1b.pdf

## You have reached the end of example 1b ✅

</details>
