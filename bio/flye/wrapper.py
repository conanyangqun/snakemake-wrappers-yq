__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import os

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
ACCEPT_INPUT_TYPE = [
    "pacbio-raw", 
    "pacbio-corr", 
    "pacbio-hifi", 
    "nano-raw", 
    "nano-corr", 
    "nano-hq",
]

# required params.
assert 'data' in snakemake.input.keys()
assert 'fasta' in snakemake.output.keys()
assert 'input_type' in snakemake.params.keys()

data = snakemake.input.data
fasta = snakemake.output.fasta
input_type = snakemake.params.input_type

if input_type not in ACCEPT_INPUT_TYPE:
    raise ValueError(f"input_type must be of {ACCEPT_INPUT_TYPE}, but got {input_type}")

# optional params.
extras = snakemake.params.get("extras", "")
flye = 'flye' if 'flye' not in snakemake.params else os.path.abspath(snakemake.params.flye)
threads = snakemake.threads if snakemake.threads else 1

# process params.
input_cmd = ""
if ',' in data:
    input_cmd = f"--{input_type} {' '.join(data.split(','))}"
elif isinstance(data, str):
    input_cmd = f"--{input_type} {data}"
else:
    raise ValueError(f"data must be of list or str, but got {type(data)}")

out_dir = os.path.dirname(os.path.abspath(fasta))

# run.
shell(
    "({flye} -t {threads} {extras} {input_cmd} --out-dir {out_dir}) {log}"
)
