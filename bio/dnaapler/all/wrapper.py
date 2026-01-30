__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# required params.
assert 'input_file' in snakemake.input.keys()
assert 'out_file' in snakemake.output.keys()

input_file = snakemake.input.input_file
out_file = snakemake.output.out_file

# optional params.
extras = snakemake.params.get("extras", "")
dnaapler = 'dnaapler' if 'dnaapler' not in snakemake.params else os.path.abspath(snakemake.params.dnaapler)
threads = snakemake.threads if snakemake.threads else 1

# input type and output type must be same.
input_ext = os.path.splitext(input_file)[1]
if input_ext.lower() not in ['.fasta', '.gfa']:
    raise ValueError(f"input_file must be of fasta or gfa, but got {input_ext}")
out_ext = os.path.splitext(out_file)[1]
if out_ext.lower() not in ['.fasta', '.gfa']:
    raise ValueError(f"out_file must be of fasta or gfa, but got {out_ext}")
if input_ext.lower() != out_ext.lower():
    raise ValueError(f"input_ext must be same as out_ext, but got {input_ext} != {out_ext}")

# out_dir and prefix.
out_dir = os.path.dirname(os.path.abspath(out_file))
out_prefix = os.path.basename(out_file).replace(f'_reoriented{out_ext}', '')
extras += f"--prefix {out_prefix} " if '--prefix' or '-p' not in extras else ""

# run.
shell(
    "({dnaapler} all -f -t {threads} {extras} --input {input_file} --output {out_dir}) {log}"
)
