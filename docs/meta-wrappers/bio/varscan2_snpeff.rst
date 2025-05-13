.. _`Varscan2 - SnpEff`:

VARSCAN2 - SNPEFF
=================

This meta-wrapper includes the following steps:

  +----------------+----------+------------------------------------------------+
  | Step           | Tool     | Reason                                         |
  +================+==========+================================================+
  | Get reference  | Ensembl  | Acquire genome sequence                        |
  +----------------+----------+------------------------------------------------+
  | Call variants  | Samtools | Acquire a list of potential variants           |
  +----------------+----------+------------------------------------------------+
  | Call variants  | Varscan2 | Call variants with Varscan2                    |
  +----------------+----------+------------------------------------------------+
  | Annotate calls | SnpEff   | Download variants database                     |
  +----------------+----------+------------------------------------------------+
  | Annotate calls | SnpEff   | Annotate variants with downloaded database     |
  +----------------+----------+------------------------------------------------+



Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule get_genome_fasta:
        output:
            "genome.fasta",
        threads: 1
        log:
            "logs/get_genome_fasta.log",
        params:
            species="saccharomyces_cerevisiae",
            datatype="dna",
            build="R64-1-1",
            release="105",
        wrapper:
            "v6.1.0-15-g668257f4/bio/reference/ensembl-sequence"


    rule samtools_mpileup:
        input:
            bam=expand(
                "{sample}.sorted.bam",
                sample=("a", "b"),
            ),
            reference_genome="genome.fasta",
        output:
            "samtools/mpileup.gz",
        threads: 1
        log:
            "logs/samtools_mpileup.log",
        params:
            extra=" --count-orphans ",
        wrapper:
            "v6.1.0-15-g668257f4/bio/samtools/mpileup"


    rule varscan2_somatic:
        input:
            mpileup="samtools/mpileup.gz",
        output:
            snp="varscan2/snp.vcf",
            indel="varscan2/indel.vcf",
        threads: 1
        log:
            "logs/varscan2_somatic.log",
        params:
            extra=" --strand-filter 1 ",
        wrapper:
            "v6.1.0-15-g668257f4/bio/varscan/somatic"


    rule snpeff_download:
        output:
            directory("resources/snpeff/R64-1-1.105"),
        threads: 1
        log:
            "logs/snpeff_download.log",
        params:
            reference="R64-1-1.105",
        wrapper:
            "v6.1.0-15-g668257f4/bio/snpeff/download"


    rule snpeff_annotate:
        input:
            calls="varscan2/snp.vcf",
            db="resources/snpeff/R64-1-1.105",
        output:
            calls="snpeff/annotated.vcf",
            stats="snpeff/annotated.html",
            csvstats="snpeff/annotated.csv",
        threads: 1
        log:
            "logs/snpeff_annotate.log",
        params:
            extra=" -nodownload -noLog ",
        wrapper:
            "v6.1.0-15-g668257f4/bio/snpeff/annotate"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/reference/ensembl-sequence`

* :ref:`bio/samtools/mpileup`

* :ref:`bio/varscan/somatic`

* :ref:`bio/snpeff/download`

* :ref:`bio/snpeff/annotate`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

