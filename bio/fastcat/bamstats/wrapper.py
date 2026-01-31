__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# required params.
assert 'bam' in snakemake.input.keys()
assert 'tsv' in snakemake.output.keys()

# optional params.
extras = snakemake.params.get("extras", "")
bamstats = 'bamstats' if 'bamstats' not in snakemake.params else os.path.abspath(snakemake.params.bamstats)
threads = snakemake.threads if snakemake.threads else 1

# process params.
bam = os.path.abspath(snakemake.input.bam)
tsv = os.path.abspath(snakemake.output.tsv)
out_dir = os.path.dirname(tsv)

# histograms.
if '--histograms' not in extras:
    extras += ' --histograms ' + os.path.join(out_dir, 'bamstats-histograms')
else:
    pass

# run.
shell(
    "({bamstats} -t {threads} {extras} {bam} >{tsv}) {log}"
)
