.. _`calc_consensus_reads`:

CALC_CONSENSUS_READS
====================

Performs consensus read calculation on a marked bam and merges single, 
paired and skipped-consensus reads into a single sorted bam file.



Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule calc_consensus_reads:
        input:
            # sorted bam file
            "mapped/{sample}.marked.bam",
        output:
            # non-overlapping consensus read pairs will be written into consensus_r1 and consensus_r2
            consensus_r1=temp("results/consensus_reads/{sample}.1.fq"),
            consensus_r2=temp("results/consensus_reads/{sample}.2.fq"),
            # consensus reads from single end records or overlapping read pairs will be merged into a single end record
            consensus_se=temp("results/consensus_reads/{sample}.se.fq"),
            # skipped reads (soft-clipped or unpropper mapped reads) will be skipped and unmarked
            skipped=temp("results/consensus_reads/{sample}.skipped.bam"),
        params:
            extra="",
        log:
            "logs/consensus/{sample}.log",
        wrapper:
            "v6.1.0-15-g668257f4/bio/rbt/collapse_reads_to_fragments-bam"


    rule map_consensus_reads:
        input:
            reads=lambda wc: expand(
                "results/consensus_reads/{sample}.{read}.fq",
                sample=wc.sample,
                read="se" if wc.read_type == "se" else (1, 2),
            ),
            idx=multiext("resources/genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        output:
            temp("results/consensus_mapped/{sample}.{read_type}.bam"),
        params:
            extra=r"-C -R '@RG\tID:{sample}\tSM:{sample}'",
            index=lambda w, input: os.path.splitext(input.idx[0])[0],
            sort="samtools",
            sort_order="coordinate",
        wildcard_constraints:
            read_type="pe|se",
        log:
            "logs/bwa_mem/{sample}.{read_type}.consensus.log",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/bwa/mem"


    rule sort_skipped_reads:
        input:
            "results/consensus_reads/{sample}.skipped.bam",
        output:
            temp("results/consensus_reads/{sample}.skipped.sorted.bam"),
        params:
            extra="-m 4G",
            tmp_dir="/tmp/",
        log:
            "logs/sort_consensus/{sample}.log",
        # Samtools takes additional threads through its option -@
        threads: 8  # This value - 1 will be sent to -@.
        wrapper:
            "v6.1.0-15-g668257f4/bio/samtools/sort"


    rule mark_duplicates_skipped:
        input:
            bams=["results/consensus_reads/{sample}.skipped.sorted.bam"],
        output:
            bam=temp("results/consensus_dupmarked/{sample}.skipped.marked.bam"),
            metrics="results/consensus_dupmarked/{sample}.skipped.metrics.txt",
        log:
            "logs/picard/marked/{sample}.log",
        params:
            extra="--VALIDATION_STRINGENCY LENIENT --TAG_DUPLICATE_SET_MEMBERS true",
        # optional specification of memory usage of the JVM that snakemake will respect with global
        # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
        # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
        # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
        resources:
            mem_mb=1024,
        wrapper:
            "v6.1.0-15-g668257f4/bio/picard/markduplicates"


    rule merge_consensus_reads:
        input:
            "results/consensus_dupmarked/{sample}.skipped.marked.bam",
            "results/consensus_mapped/{sample}.se.bam",
            "results/consensus_mapped/{sample}.pe.bam",
        output:
            "results/consensus/{sample}.bam",
        log:
            "logs/samtools_merge/{sample}.log",
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/samtools/merge"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/rbt/collapse_reads_to_fragments-bam`

* :ref:`bio/bwa/mem`

* :ref:`bio/bwa/index`

* :ref:`bio/samtools/merge`

* :ref:`bio/samtools/sort`

* :ref:`bio/picard/markduplicates`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


