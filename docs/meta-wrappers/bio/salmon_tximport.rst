.. _`Salmon Tximport`:

SALMON TXIMPORT
===============

This meta-wrapper includes the following steps:

  +----------------+----------+------------------------------------------------+
  | Step           | Tool     | Reason                                         |
  +================+==========+================================================+
  | Indexation     | Bash     | Identify decoy sequences                       |
  +----------------+----------+------------------------------------------------+
  | Indexation     | Salmon   | Create decoy aware gentrome                    |
  |                |          | (genome + trancriptome) index                  |
  +----------------+----------+------------------------------------------------+
  | Quantification | Salmon   | Quantify sequenced reads                       |
  +----------------+----------+------------------------------------------------+
  | Quantification | Tximport | Import counts and inferential replicates in R  |
  |                |          | as a ready-to-use SummarizedExperiment object. |
  +----------------+----------+------------------------------------------------+



Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule salmon_decoy_sequences:
        input:
            transcriptome="resources/transcriptome.fasta",
            genome="resources/genome.fasta",
        output:
            gentrome=temp("resources/gentrome.fasta"),
            decoys=temp("resources/decoys.txt"),
        threads: 1
        log:
            "decoys.log",
        wrapper:
            "v6.1.0-15-g668257f4/bio/salmon/decoys"


    rule salmon_index_gentrome:
        input:
            sequences="resources/gentrome.fasta",
            decoys="resources/decoys.txt",
        output:
            multiext(
                "salmon/transcriptome_index/",
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
        cache: True
        log:
            "logs/salmon/transcriptome_index.log",
        threads: 2
        params:
            # optional parameters
            extra="",
        wrapper:
            "v6.1.0-15-g668257f4/bio/salmon/index"


    rule salmon_quant_reads:
        input:
            r="reads/{sample}.fastq.gz",
            index=multiext(
                "salmon/transcriptome_index/",
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            ),
            gtf="resources/annotation.gtf",
        output:
            quant=temp("pseudo_mapping/{sample}/quant.sf"),
            quant_gene=temp("pseudo_mapping/{sample}/quant.genes.sf"),
            lib=temp("pseudo_mapping/{sample}/lib_format_counts.json"),
            aux_info=temp(directory("pseudo_mapping/{sample}/aux_info")),
            cmd_info=temp("pseudo_mapping/{sample}/cmd_info.json"),
            libparams=temp(directory("pseudo_mapping/{sample}/libParams")),
            logs=temp(directory("pseudo_mapping/{sample}/logs")),
        log:
            "logs/salmon/{sample}.log",
        params:
            # optional parameters
            libtype="A",
            extra="--numBootstraps 32",
        threads: 2
        wrapper:
            "v6.1.0-15-g668257f4/bio/salmon/quant"


    rule tximport:
        input:
            quant=expand(
                "pseudo_mapping/{sample}/quant.sf", sample=["S1", "S2", "S3", "S4"]
            ),
            lib=expand(
                "pseudo_mapping/{sample}/lib_format_counts.json",
                sample=["S1", "S2", "S3", "S4"],
            ),
            aux_info=expand(
                "pseudo_mapping/{sample}/aux_info", sample=["S1", "S2", "S3", "S4"]
            ),
            cmd_info=expand(
                "pseudo_mapping/{sample}/cmd_info.json", sample=["S1", "S2", "S3", "S4"]
            ),
            libparams=expand(
                "pseudo_mapping/{sample}/libParams", sample=["S1", "S2", "S3", "S4"]
            ),
            logs=expand("pseudo_mapping/{sample}/logs", sample=["S1", "S2", "S3", "S4"]),
            tx_to_gene="resources/tx2gene.tsv",
        output:
            txi="tximport/SummarizedExperimentObject.RDS",
        params:
            extra="type='salmon'",
        log:
            "logs/tximport.log"
        wrapper:
            "v6.1.0-15-g668257f4/bio/tximport"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/salmon/decoys`

* :ref:`bio/salmon/index`

* :ref:`bio/salmon/quant`

* :ref:`bio/tximport`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

