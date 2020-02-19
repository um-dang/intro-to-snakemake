# Snakemake RNAseq example 

This is a simple bioinformatics workflow example. This example goes a bit more in-depth with snakemake features such as environment management and cluster submission that are very useful when e.g. running them on a high-performance compute cluster with SLURM. This folder contains some simulated RNAseq read data which maps to hg38 (human) chromosome 22. The chr22 reference from ENSEMBL is also included in this repository.

##### Note: The basics of snakemake are covered in the [toy example](../toy_example/).

Main topics covered in this example:
1. Using Singularity for environment management
2. Submitting workflow jobs to the SLURM cluster


Get started by expanding Example 3a:

<details><summary>Expand - Ex. 3a</summary>

Here we'll make a simple workflow for aligning the reads. There are two rules, one for creating a bwa index from the given fasta, and one to use bwa mem to align the reads to the index. 

##### Note: You'll have to unzip the Chr22 reference before running the pipeline - Otherwise snakemake will fail to find the required input file (doesn't expect .gz)

Create a file named rnaseq.snakefile with the following contents:

    rule all:
        input:
            expand("aligned/{sample}.sam", sample=config['samples'])

    rule bwa_index:
        input:
            config['ref_fasta']
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
            
            
Also create a config file named rnaseq_config.yml with the following contents:

    ref_fasta: Homo_sapiens.GRCh38.dna_sm.chr22.fa
    ref_prefix: reference/chr22
    samples:
      - sample_01
      - sample_02
      - sample_03
      - sample_04
      - sample_05
      - sample_06
      - sample_07
      - sample_08


##### Things to notice/think about:
* Outputs for bwa_index rule are pretty verbose. Can we use an expand statement to define these more succinctly?
* Named inputs are used in bwa_mem step. How is this useful?
* `reference/chr22.sa` is listed as an input for bwa_mem. Why?
  * bonus: It's listed explicitly here to illustrate the above point. Maybe it should be expressed in a different way. How else could this be written?
* We're setting number of threads on a per-rule basis (for bwa_mem)

Perform a dry-run

    snakemake --snakefile rnaseq.snakefile --configfile rnaseq_config.yml --dry-run
    
Now try running it (run the previous command without the --dry-run flag).
What happens? Why?

## You have reached the end of example 3a ✅
            
</details>



Now we'll make a few changes to enable us to run this on the compute cluster

<details><summary>Expand - Ex. 3b</summary>

Modify rnaseq.snakefile so that it has the following contents:

    rule all:
        input:
            expand("aligned/{sample}.sam", sample=config['samples'])

    rule bwa_index:
        input:
            config['ref_fasta']
        output:
            expand(config['ref_prefix'] + ".{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
        params:
            ref_prefix = config['ref_prefix']
        singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
        shell:
            "bwa index -p {params.ref_prefix} {input}"

    rule bwa_mem:
        input:
            sa = config['ref_prefix'] + ".sa",
            r1 = "read_data/{sample}/{sample}_R1.fastq.gz",
            r2 = "read_data/{sample}/{sample}_R2.fastq.gz"
        output:
            "aligned/{sample}.sam"
        threads: 10
        params:
            ref_prefix = config['ref_prefix']
        singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
        shell:
            "bwa mem -t {threads} {params.ref_prefix} {input.r1} {input.r2} > {output}"

Create a cluster configuration file called greatlakes_config.yml with the following contents:

    __default__:
        name: '{rule}_{wildcards}'
        account: cgates1
        partition: standard
        nodes: '1'
        ntask: '1'
        memory: '1024'
        time:   '00:30:00'

    bwa_mem:
        ntask: '{threads}'
        

##### Things to notice/think about:

Expand statement
* Compare this examples' expand statement with the previous examples' verbose output definitions

Environment management
* We're using [biocontainers](https://biocontainers.pro/#/) along with snakemake's `singularity` directive. Consider granularity of software needs

Cluster configuration
* The \_\_default\_\_ configuration
* Rule-specific configurations override the defaults
* We added more threads to the process - will be requested at job submission and used during execution



##### Running the workflow on the GreatLakes cluster:

I'll give the whole command-line invocation first, and then explain below (it may seem complex at first glance)

    snakemake --snakefile rnaseq.snakefile --configfile rnaseq_config.yml --use-singularity --jobs 144 --cluster-config greatlakes_config.yml --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time}'
    
Breaking down the new additions to the command line invocation:
* `--use-singularity` tells snakemake to make use of the `singularity:` blocks within the rules
* `--jobs` must be used for cluster submission; it's the number of concurrent jobs sent to the scheduler
* `--cluster-config` provides the name of the file where job submission configuration details are defined
* `--cluster` is a template string which is filled-in with the values from greatlakes_config.yml, in order to produce the job submission commands for the cluster. This allows the flexibility of working with other resource manager/schedulers by modifying the template.


## You have reached the end of example 3b ✅

</details>

Finally, in example 3c we extend the workflow to produce bigwig files - display tracks which can be opened in a genome browser to visualize the alignment results. If you're inclined, try to achieve this on your own before expanding the example below.

Hint: deeptools bamCoverage is a great tool for producing the bigwig. Intermediate steps are required, though.

<details><summary>Expand - Ex. 3c</summary>

rnaseq.snakefile contents:

    rule all:
        input:
            expand("coverage/{sample}.bigwig", sample=config['samples'])

    rule bwa_index:
        input:
            config['ref_fasta']
        output:
            expand(config['ref_prefix'] + ".{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
        params:
            ref_prefix = config['ref_prefix']
        singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
        shell:
            "bwa index -p {params.ref_prefix} {input}"

    rule bwa_mem:
        input:
            sa = "reference/chr22.sa",
            r1 = "read_data/{sample}/{sample}_R1.fastq.gz",
            r2 = "read_data/{sample}/{sample}_R2.fastq.gz"
        output:
            "aligned/{sample}.sam"
        threads: 10
        params:
            ref_prefix = config['ref_prefix']
        singularity: "docker://biocontainers/bwa:v0.7.15_cv4"
        shell:
            "bwa mem -t {threads} {params.ref_prefix} {input.r1} {input.r2} > {output}"

    rule sam_to_bam_bai:
        input:
            "aligned/{sample}.sam"
        output:
            bam = "aligned/{sample}.bam",
            bai = "aligned/{sample}.bai"
        singularity: "docker://biocontainers/samtools:v1.9-4-deb_cv1"
        shell:
            "samtools sort -o {output.bam} {input} && "
            "samtools index {output.bam} {output.bai}"

    rule deeptools_bamCoverage:
        input:
            bam = "aligned/{sample}.bam",
            bai = "aligned/{sample}.bai"
        output:
            "coverage/{sample}.bigwig"
        singularity: "docker://quay.io/biocontainers/deeptools:3.0.1--py36_1"
        shell:
            "bamCoverage -b {input.bam} -o {output}"


## You have reached the end of example 3c ✅

</details>
