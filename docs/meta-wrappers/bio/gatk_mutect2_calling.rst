.. _`GATK Variant calling best practice workflow`:

GATK VARIANT CALLING BEST PRACTICE WORKFLOW
===========================================

Call short variants (SNP+INDEL) with GATK's Mutect2:

  +----------------+---------------------------+-------------------------------------------------------------+
  | Step           | Tool                      | Reason                                                      |
  +================+===========================+=============================================================+
  +----------------+---------------------------+-------------------------------------------------------------+
  | Indexing       | Picard                    | Create genome sequence dictionnary                          |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Indexing       | Samtools                  | Index fasta genome sequence                                 |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Grouping       | Picard                    | Add or replace possible missing read groups                 |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Indexing       | Sambamba                  | Index re-grouped BAM-formatted alignments                   |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Calling        | Mutect2                   | Call short variants with Mutect2                            |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Contaminations | GetPileupSummaries        | Tabulates pileup metrics for inferring contamination        |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Contaminations | CalculateContamination    | Estimate cross sample contamination                         |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Orientation    | LearnReadOrientationModel | Search for sequencing artifacts based on read orientation   |
  +----------------+---------------------------+-------------------------------------------------------------+
  | Filtering      | FilterMutectCalls         | Use previously estimated biases and filter variants         |
  +----------------+---------------------------+-------------------------------------------------------------+



Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule create_dict:
        input:
            "genome.fasta",
        output:
            "genome.dict",
        threads: 1
        resources:
            mem_mb=1024,
        log:
            "logs/picard/create_dict.log",
        params:
            extra="",
        wrapper:
            "v6.1.0-15-g668257f4/bio/picard/createsequencedictionary"


    rule samtools_index:
        input:
            "genome.fasta",
        output:
            "genome.fasta.fai",
        log:
            "logs/genome_index.log",
        params:
            extra="",  # optional params string
        wrapper:
            "v6.1.0-15-g668257f4/bio/samtools/faidx"


    rule picard_replace_read_groups:
        input:
            "mapped/{sample}.bam",
        output:
            "picard/{sample}.bam",
        threads: 1
        resources:
            mem_mb=1024,
        log:
            "logs/picard/replace_rg/{sample}.log",
        params:
            # Required for GATK
            extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
        wrapper:
            "v6.1.0-15-g668257f4/bio/picard/addorreplacereadgroups"


    rule sambamba_index_picard_bam:
        input:
            "picard/{sample}.bam",
        output:
            "picard/{sample}.bam.bai",
        threads: 1
        log:
            "logs/sambamba/index/{sample}.log",
        params:
            extra="",
        wrapper:
            "v6.1.0-15-g668257f4/bio/sambamba/index"


    rule mutect2_call:
        input:
            fasta="genome.fasta",
            fasta_dict="genome.dict",
            fasta_fai="genome.fasta.fai",
            map="picard/{sample}.bam",
            map_idx="picard/{sample}.bam.bai",
            intervals="regions.bed",
        output:
            vcf="variant/{sample}.vcf",
            bam="variant/{sample}.bam",
            f1r2="counts/{sample}.f1r2.tar.gz",
        threads: 1
        resources:
            mem_mb=1024,
        params:
            extra=" --tumor-sample {sample} ",
        log:
            "logs/mutect/{sample}.log",
        wrapper:
            "v6.1.0-15-g668257f4/bio/gatk/mutect"


    rule gatk_get_pileup_summaries:
        input:
            bam="picard/{sample}.bam",
            bai_bai="picard/{sample}.bam.bai",
            variants="known.vcf.gz",
            variants_tbi="known.vcf.gz.tbi",
            intervals="regions.bed",
        output:
            temp("summaries/{sample}.table"),
        threads: 1
        resources:
            mem_mb=1024,
        params:
            extra="",
        log:
            "logs/summary/{sample}.log",
        wrapper:
            "v6.1.0-15-g668257f4/bio/gatk/getpileupsummaries"


    rule gatk_calculate_contamination:
        input:
            tumor="summaries/{sample}.table",
        output:
            temp("contamination/{sample}.pileups.table"),
        threads: 1
        resources:
            mem_mb=1024,
        log:
            "logs/contamination/{sample}.log",
        params:
            extra="",
        wrapper:
            "v6.1.0-15-g668257f4/bio/gatk/calculatecontamination"


    rule gatk_learn_read_orientation_model:
        input:
            f1r2="counts/{sample}.f1r2.tar.gz",
        output:
            temp("artifacts_prior/{sample}.tar.gz"),
        threads: 1
        resources:
            mem_mb=1024,
        params:
            extra="",
        log:
            "logs/learnreadorientationbias/{sample}.log",
        wrapper:
            "v6.1.0-15-g668257f4/bio/gatk/learnreadorientationmodel"


    rule filter_mutect_calls:
        input:
            vcf="variant/{sample}.vcf",
            ref="genome.fasta",
            ref_dict="genome.dict",
            ref_fai="genome.fasta.fai",
            bam="picard/{sample}.bam",
            bam_bai="picard/{sample}.bam.bai",
            contamination="contamination/{sample}.pileups.table",
            f1r2="artifacts_prior/{sample}.tar.gz",
        output:
            vcf="variant/{sample}.filtered.vcf.gz",
            vcf_idx="variant/{sample}.filtered.vcf.gz.tbi",
        threads: 1
        resources:
            mem_mb=1024,
        log:
            "logs/gatk/filter/{sample}.log",
        params:
            extra="--create-output-variant-index --min-median-mapping-quality 35 --max-alt-allele-count 3",
            java_opts="",
        wrapper:
            "v6.1.0-15-g668257f4/bio/gatk/filtermutectcalls"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/samtools/faidx`

* :ref:`bio/picard/createsequencedictionary`

* :ref:`bio/sambamba/index`

* :ref:`bio/picard/addorreplacereadgroups`

* :ref:`bio/gatk/mutect`

* :ref:`bio/gatk/getpileupsummaries`

* :ref:`bio/gatk/calculatecontamination`

* :ref:`bio/gatk/learnreadorientationmodel`

* :ref:`bio/gatk/filtermutectcalls`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

