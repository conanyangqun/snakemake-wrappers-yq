# mosdepth wrapper.
__author__ = "yangqun"

import sys
import os
from snakemake.shell import shell

# params
assert 'bam' in snakemake.input.keys()
assert 'bai' in snakemake.input.keys()
assert 'summary' in snakemake.output.keys()

bam = snakemake.input.get('bam')
bai = snakemake.input.get('bai')
summary = snakemake.output.get('summary')
prefix = summary.replace(".mosdepth.summary.txt", "")

# optional
mosdepth = snakemake.params.get('mosdepth', 'mosdepth')
by = snakemake.params.get('by', '')
if by:
    by = '--by ' + by

thresholds = snakemake.params.get('thresholds', '')
if thresholds:
    thresholds = '--thresholds ' + thresholds

extras = snakemake.params.get('extras', '')

# threads
threads = snakemake.threads if snakemake.threads <= 4 else 4

# run
shell(
    "{mosdepth} -t {threads} {by} {thresholds} {extras} {prefix} {bam}"
)
