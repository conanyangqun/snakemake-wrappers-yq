.. _`Bowtie2 Sambamba`:

BOWTIE2 SAMBAMBA
================

Map reads with Bowtie 2, and post-process alignments with Sambamba:

  +----------------+----------+-----------------------------------------+
  | Step           | Tool     | Reason                                  |
  +================+==========+=========================================+
  | Indexation     | Bowtie 2 | Create genome sequence index            |
  +----------------+----------+-----------------------------------------+
  | Mapping        | Bowtie 2 | Perform read mapping                    |
  +----------------+----------+-----------------------------------------+
  | Sort           | Sambamba | Perform sort based on mapping position  |
  +----------------+----------+-----------------------------------------+
  | Quality filter | Sambamba | Perform mapping quality filter          |
  +----------------+----------+-----------------------------------------+
  | Deduplication  | Sambamba | Identify possible sequencing duplicates |
  +----------------+----------+-----------------------------------------+
  | Indexation     | Sambamba | Index deduplicated reads                |
  +----------------+----------+-----------------------------------------+



Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule bowtie2_build:
        input:
            ref="genome.fasta",
        output:
            multiext(
                "genome",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        log:
            "logs/bowtie2_build/build.log",
        params:
            extra="",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/bowtie2/build"


    rule bowtie2_alignment:
        input:
            sample=["{sample}.R1.fastq", "{sample}.R2.fastq"],
            idx=multiext(
                "genome",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        output:
            temp("mapped/{sample}.bam"),
        log:
            "logs/bowtie2/{sample}.log",
        params:
            extra=(
                " --rg-id {sample} "
                "--rg 'SM:{sample} LB:FakeLib PU:FakePU.1.{sample} PL:ILLUMINA' "
            ),
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/bowtie2/align"


    rule sambamba_sort:
        input:
            "mapped/{sample}.bam",
        output:
            temp("mapped/{sample}.sorted.bam"),
        params:
            "",
        log:
            "logs/sambamba-sort/{sample}.log",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/sambamba/sort"


    rule sambamba_view:
        input:
            "mapped/{sample}.sorted.bam",
        output:
            temp("mapped/{sample}.filtered.bam"),
        params:
            extra=(
                " --format 'bam' "
                "--filter 'mapping_quality >= 30 and not (unmapped or mate_is_unmapped)' "
            ),
        log:
            "logs/sambamba-view/{sample}.log",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/sambamba/view"


    rule sambamba_markdup:
        input:
            "mapped/{sample}.filtered.bam",
        output:
            "mapped/{sample}.rmdup.bam",
        params:
            extra=" --remove-duplicates ",  # optional parameters
        log:
            "logs/sambamba-markdup/{sample}.log",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/sambamba/markdup"


    rule sambamba_index:
        input:
            "mapped/{sample}.rmdup.bam",
        output:
            "mapped/{sample}.rmdup.bam.bai",
        params:
            extra="",
        log:
            "logs/sambamba-index/{sample}.log",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/sambamba/index"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bowtie2/build`

* :ref:`bio/bowtie2/align`

* :ref:`bio/sambamba/sort`

* :ref:`bio/sambamba/view`

* :ref:`bio/sambamba/markdup`

* :ref:`bio/sambamba/index`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

