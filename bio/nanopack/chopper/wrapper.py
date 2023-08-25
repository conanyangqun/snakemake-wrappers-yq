# chopper from nanopack.
__author__ = "yangqun"
__github__ = "https://github.com/conanyangqun/snakemake-wrappers"

import os
import sys

from snakemake.shell import shell

# params.
if 'fastqs' not in snakemake.input.keys():
    sys.exit("fastqs must be given.")
if 'out' not in snakemake.output.keys():
    sys.exit('out must be given.')

fastqs = snakemake.input.get('fastqs')
out = snakemake.output.get('out')

# check input file extension.
exts = set([os.path.splitext(x)[-1] for x in fastqs])
if len(exts) > 1:
    sys.exit('input file must be only 1 file type.')

ext = exts.pop()
ext = ext.lower()

# optional params
chopper = snakemake.params.get('chopper', 'chopper')
quality = snakemake.params.get('quality', '0')
min_len = snakemake.params.get('min_len', '1')
max_len = snakemake.params.get('max_len', '2147483647')
headcrop = snakemake.params.get('headcrop', '0')
tailcrop = snakemake.params.get('tailcrop', '0')
extras = snakemake.params.get('extras', '')

# commands.
if ext == '.gz':
    input_cmd = 'zcat {}'.format(' '.join(fastqs))
else:
    input_cmd = 'cat {}'.format(' '.join(fastqs))

out_cmd = f' | gzip >{out}' if out.endswith('.gz') else f' >{out}'
threads_cmd = f'--threads {snakemake.threads}'

shell(
    "{input_cmd} "
    "| {chopper} {threads_cmd} --maxlength {max_len} --minlength {min_len} "
    "--headcrop {headcrop} --tailcrop {tailcrop} --quality {quality} "
    "{out_cmd}"
)
