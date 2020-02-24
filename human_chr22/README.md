# Snakemake RNAseq example

This is an example of a simple bioinformatics workflow. This example covers a few important and useful snakemake features such as environment management and job submission which enable possibilities e.g. running workflows on a high-performance compute cluster with SLURM. This folder contains some simulated RNAseq read data which maps to hg38 (human) chromosome 22. The chr22 reference from ENSEMBL is also included in this repository.

##### Note: The basics of snakemake are covered in the [toy example](../toy_example/).

Main topics covered in this example:
1. Using Singularity for environment management
2. Submitting workflow jobs to the SLURM cluster


Get started by expanding Example 2a:

<details><summary>Expand - Ex. 2a</summary>

Here we'll make a simple workflow for aligning the reads. There are three rules, one for unzipping the fasta reference, one for creating a bwa index from the given fasta, and one to use bwa mem to align the reads to the index.

#### Note: You'll have to complete the first rule in order for this to work!

ex-2a.smk has the following contents:

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


config-ex2.yml has the following contents:

    ref_fasta_gz: Homo_sapiens.GRCh38.dna_sm.chr22.fa.gz
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


Review these files and consider the following topics:
* Named inputs/outputs
* Linking rules with dependencies
* Specifying resources on a per-rule basis


Perform a dry-run

    snakemake --snakefile ex-2a.smk --configfile config-ex2.yml --dry-run

Now try running it (run the previous command without the --dry-run flag).
What happens? Why?

## You have reached the end of example 2a ✅

</details>


Now we'll make a few changes to enable us to run this on the compute cluster

<details><summary>Expand - Ex. 2b</summary>

ex-2b.smk has the following contents:

    rule all:
        input:
            expand("aligned/{sample}.sam", sample=config['samples'])

    rule fasta_gunzip:
        input:
            config['ref_fasta_gz']
        output:
            config['ref_prefix'] + ".fa"
        shell:
            "gunzip -c {input} > {output}"

    rule bwa_index:
        input:
            config['ref_prefix'] + ".fa"
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


config-ex2.yml is unmodified from the previous example

jobsub-config.yml has the following contents:

    __default__:
        name: '{rule}_{wildcards}'
        account: your_account
        partition: standard
        nodes: '1'
        ntask: '1'
        memory: '1024'
        time:   '00:30:00'

    bwa_mem:
        ntask: '{threads}'


Review these files and consider the following:
* Snakemake + Singularity
* Granularity of software needs
* The job-submission configuration
* Compare this examples' expand statement with the previous examples' verbose output definitions


Running the workflow on the GreatLakes cluster:

#### Note: The singularity executable must be available. On GreatLakes, `module load singularity`
#### Note: You'll have to change the account name in `jobsub-config.yml` to your own greatlakes account

I'll give the whole command-line invocation first, and then explain below (it may seem complex at first glance)

    module load singularity
    snakemake --snakefile ex-2b.smk --configfile config-ex2.yml --use-singularity --jobs 144 --cluster-config jobsub-config.yml --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time}'

Breaking down the new additions to the command line invocation:
* `--use-singularity` tells snakemake to make use of the `singularity:` blocks within the rules
* `--jobs` must be used for cluster submission; it's the number of concurrent jobs sent to the scheduler
* `--cluster-config` provides the name of the file where job submission configuration details are defined
* `--cluster` is a template string which is filled-in with the values from jobsub-config.yml, in order to produce the job submission commands for the cluster. This allows the flexibility of working with other resource manager/schedulers by modifying the template.


## You have reached the end of example 2b ✅

</details>

Finally, example 2c is an open-ended exercise, where you can extend the workflow to generate additional targets.

Before expanding example 2c, try creating your own by extending example 2b.

    cp ex-2b.smk ex-2c.smk

Edit ex-2c.smk to achieve the following:
* extend the workflow to produce bigwig files - display tracks which can be opened in a genome browser to visualize the alignment results.
* Hint: deeptools bamCoverage is a great tool for producing the bigwig. Intermediate steps are required, though. It only accepts sorted and indexed bam as input

<details><summary>Expand - Ex. 2c</summary>

Before expanding the contents of ex-2c.smk, try copying ex-2b.smk

ex-2c.smk has the following contents:

    rule all:
        input:
            expand("coverage/{sample}.bigwig", sample=config['samples'])

    rule fasta_gunzip:
        input:
            config['ref_fasta_gz']
        output:
            temp(config['ref_prefix'] + ".fa")
        shell:
            "gunzip -c {input} > {output}"

    rule bwa_index:
        input:
            config['ref_prefix'] + ".fa"
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


Example 2c dry-run:

    snakemake --snakefile ex-2c.smk --configfile config-ex2.yml --dry-run

Running example 2c on GreatLakes:

        snakemake --snakefile ex-2c.smk --configfile config-ex2.yml --use-singularity --jobs 144 --cluster-config jobsub-config.yml --cluster 'sbatch --job-name={cluster.name} --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks-per-node={cluster.ntask} --mem={cluster.memory} --time={cluster.time}'


## You have reached the end of example 2c ✅

</details>
