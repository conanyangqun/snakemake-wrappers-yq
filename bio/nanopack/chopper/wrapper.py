# chopper from nanopack.
__author__ = "yangqun"
__github__ = "https://github.com/conanyangqun/snakemake-wrappers"

import os
import sys

from snakemake.shell import shell

# params.
if 'fastq' not in snakemake.input.keys():
    sys.exit("fastq must be given.")
if 'out' not in snakemake.output.keys():
    sys.exit('out must be given.')

fastq = snakemake.input.get('fastq')
out = snakemake.output.get('out')

# optional params
chopper = snakemake.params.get('chopper', 'chopper')
quality = snakemake.params.get('quality', '0')
min_len = snakemake.params.get('min_len', '1')
max_len = snakemake.params.get('max_len', '2147483647')
headcrop = snakemake.params.get('headcrop', '0')
tailcrop = snakemake.params.get('tailcrop', '0')
extras = snakemake.params.get('extras', '')

# commands.
input_cmd = f'zcat {fastq}' if fastq.endswith('.gz') else f'cat {fastq}'
out_cmd = f' | gzip >{out}' if out.endswith('.gz') else f' >{out}'
threads_cmd = f'--threads {snakemake.threads}'

shell(
    "{input_cmd} "
    "| {chopper} {threads_cmd} --maxlength {max_len} --minlength {min_len} "
    "--headcrop {headcrop} --tailcrop {tailcrop} --quality {quality} "
    "{out_cmd}"
)
