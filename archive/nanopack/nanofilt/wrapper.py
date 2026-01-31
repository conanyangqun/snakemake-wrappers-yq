# nanofilt wrapper
__author__ = "yangqun"

import sys
from snakemake.shell import shell

# arguments
assert 'fastq' in snakemake.input.keys()
assert 'out' in snakemake.output.keys()
fastq = snakemake.input.get('fastq')
out = snakemake.output.get('out')

# optional arguments
min_len = snakemake.params.get('min_len')
max_len = snakemake.params.get('max_len')
extras = snakemake.params.get('extras', '')
nanofilt = snakemake.params.get('nanofilt', "NanoFilt")

extract_cmd = "zcat " + fastq if fastq.endswith('gz') else "cat " + fastq
nanofilt_options = ""
if min_len:
    nanofilt_options += f" -l {min_len}"

if max_len:
    nanofilt_options += f" --maxlength {max_len}"

nanofilt_options += extras

out_cmd = f" | gzip >{out}" if out.endswith('gz') else f" >{out}"

# run
shell(
    "{extract_cmd} | {nanofilt} {nanofilt_options} {out_cmd}"
)
