# minimap2 wrapper, modified from snakemake wrapper
__author__ = "yangqun"
__email__ = "yangqunxinli@163.com"
__verison__ = "1.0.0"

import sys
import os
from snakemake.shell import shell

assert 'out' in snakemake.output.keys()
out = snakemake.output.get('out')
out_name, out_ext = os.path.splitext(out)
out_ext = out_ext[1:].upper()

extras = snakemake.params.get('extras', '')

rg = snakemake.params.get('rg', '')
rg = "" if not rg else f'-R {rg}'

sort_type = snakemake.params.get('sort_type', 'coordinate')
sort_extras = snakemake.params.get('sort_extras', '')

pipe_cmd = ""
index_cmd = ""

if out_ext not in ['SAM', 'BAM', 'PAF']:
    sys.exit("out file extension must be .sam, .bam, .paf")

if sort_type not in ['none', 'coordinate', 'queryname']:
    sys.exit('samtools sort must be none, coordinate, queryname. got {}'.format(sort_type))

if out_ext != "PAF":
    # sam/bam/cram
    extras += " -a"
    if sort_type == "none":
        # 不排序
        out_cmd = "| samtools view -h --output-fmt {} -o {}".format(out_ext, out)
    elif sort_type == "queryname":
        out_cmd = "| samtools sort -n --output-fmt {} -o {}".format(out_ext, out)
    else:
        out_cmd = "| samtools sort --output-fmt {} -o {}".format(out_ext, out)
else:
    out_cmd = "-o {}".format(out)

if out_ext in ['BAM'] and sort_type == 'coordinate':
    index_cmd = "&& samtools index {}".format(out)

shell(
    "minimap2 -t {snakemake.threads} {extras} {rg} {snakemake.input.ref_fa} {snakemake.input.fqs_fa} "
    "{out_cmd} {index_cmd}"
)
