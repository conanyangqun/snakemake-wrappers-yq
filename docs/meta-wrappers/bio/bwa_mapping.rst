.. _`bwa_mapping`:

BWA_MAPPING
===========

Map reads with bwa-mem and index with samtools index - this is just a test for subworkflows


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule bwa_mem:
        input:
            reads=["reads/{sample}.1.fastq", "reads/{sample}.2.fastq"],
            idx=multiext("genome", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        output:
            "mapped/{sample}.bam"
        log:
            "logs/bwa_mem/{sample}.log"
        params:
            extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
            sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
            sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
            sort_extra=""            # Extra args for samtools/picard.
        threads: 8
        wrapper:
            "v6.1.0-15-g668257f4/bio/bwa/mem"

    rule samtools_index:
        input:
            "mapped/{sample}.bam"
        output:
            "mapped/{sample}.bam.bai"
        log:
            "logs/samtools_index/{sample}.log"
        params:
            "" # optional params string
        wrapper:
            "v6.1.0-15-g668257f4/bio/samtools/index"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bwa/mem`

* :ref:`bio/samtools/index`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Jan Forster

