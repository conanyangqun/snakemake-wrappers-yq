# fastcat wrapper.

__author__ = "yangqun"
__email__ = "yangqunxinli@163.com"
__verison__ = "1.0.0"

import sys
import os
from snakemake.shell import shell

data = snakemake.input.get('fq_or_path')
out_fq = snakemake.output.get('out_fq')
file_summary = snakemake.output.get('file_summary')
read_summary = snakemake.output.get('read_summary')
extras = snakemake.params.get('extras', '')
fastcat = snakemake.params.get('fastcat', 'fastcat')

# input data
if not data:
    sys.exit("Please give fastq files or path.")

if os.path.isdir(data):
    extras += ' -x'

# output data
out_cmd = ''
if out_fq:
    if out_fq.endswith('gz'):
        # compressed fastq.
        out_cmd = " | gzip >{}".format(out_fq)
    else:
        # fastq.
        out_cmd = " >{}".format(out_fq)
else:
    # do not output fastq.
    out_cmd = " >/dev/null"

cmds = ""
# summary
if file_summary:
    cmds += ' -f ' + file_summary
if read_summary:
    cmds += ' -r ' + read_summary

# run
shell(
    "{fastcat} {cmds} {extras} {data} {out_cmd}"
)
