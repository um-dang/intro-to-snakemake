# Snakemake example 01

This is a simple toy example that can be used to start learning the basics of snakemake. This folder contains some data files with simple contents. For example, file 123.txt contains three lines of data, with values 1, 2, and 3. We'll use these simple data along with core command-line tools, included by default in most unix environments, to illustrate the basics of snakemake.

Note: If you haven't installed snakemake yet, I'd suggest using [conda](https://docs.conda.io/en/latest/miniconda.html) to do so.

Topics covered:
1. Setting up workflow targets and dependencies
2. Getting started writing rules
3. The workflow & DAG

<details><summary>Expand - Ex. 1a</summary>

Create a file named toy.snakefile with the following contents:

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

Topics covered:
* Targets & dependencies
* Writing rules

Targets - Have to think about the pipeline backwards - What do we want to end up with?

Dependencies - What rules must be in place (inputs/outputs) for these targets to be generated?
 - This will be more obvious with the addition of more rules

Writing rules - Generally will have 'input', 'output', and 'shell' blocks (more the whole story)

The rule 'all' is placed at the top of the file (the first rule, anyway), and is always executed. It's being used to define the targets for the workflow. I.e. by default, this workflow will generate these targets.

Now perform a dry-run:

    snakemake --snakefile toy.snakefile --dry-run

Notice that snakemake keeps track of the wildcards during the evaluation of each rule
* experiment by changing the targets so they don't match the input files
            
</details>

## ✅

Now we'll make the workflow a bit more interesting. We'll add more rules, use a configuration file, and more!

<details><summary>Expand - Ex. 1b</summary>


toy.snakefile contents:

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
        "sleep 2 ; sort -n {input} > {output}"
        
toy_config.yml contents:

    basenames:
      - '123'
      - '345'
      - '567'
    append_val: 42
        
The expand statement
* The various uses of curly braces can be confusing at first (at least for me)
* `expand` is distinct from `wildcards`
* can be thought of as "expand this string (arg 1) into an array of strings, filling in all combinations of values (using arg 2)

Running the new snakefile/configfile

    snakemake --snakefile toy.snakefile --configfile toy_config.yml
        
Did snakemake run the workflow, and successfully create the desired targets?
* View the directory of results

More about the core tenets of snakemake (also gnu make, make-like things)
* Try running the workflow to completion, then running it again. What happens? 
* Delete `output`, then try again. Isn't this cool?
* Try deleting an intermediate file, then running the pipeline again. How is this beneficial? How can it be problematic?

The workflow & DAG
* Directed Acyclic Graph - how snakemake 'knows' how to produce the desired targets
* It can be useful to see the workflow DAG, and imagine how snakemake 'thinks' about executing it

Viewing the DAG (or rulegraph)

    #DAG (file-level granularity)
    snakemake --snakefile toy.snakefile --configfile toy_config.yml --dag | dot -T pdf > toy_dag.pdf
    #Rulegraph (rule-level granularity)
    snakemake --snakefile toy.snakefile --configfile toy_config.yml --rulegraph | dot -T pdf > toy_rulegraph.pdf

</details>

## ✅
