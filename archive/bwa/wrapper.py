# bwa wrapper
__author__ = "yangqun"

import os
import sys
from snakemake.shell import shell

# input and output
assert 'ref_fa' in snakemake.input.keys()
assert 'fqs' in snakemake.input.keys()
assert 'out' in snakemake.output.keys()

# arguments
ref_fa = snakemake.input.get('ref_fa')
fqs = snakemake.input.get('fqs')
out = snakemake.output.get('out')

bwa_cmd = snakemake.params.get('bwa_cmd', 'mem -K 100000000 -v 3 -Y')
bwa_path = snakemake.params.get('bwa_path', '/usr/gitc/bwa')

rg = snakemake.params.get('rg', '')
rg = f'-R {rg}' if rg else ''

samtools = snakemake.params.get('samtools', 'samtools')
sort_type = snakemake.params.get('sort_type', 'coordinate')
sort_extras = snakemake.params.get('sort_extras', '')

bwa_threads = "" if snakemake.threads <= 1 else '-t {}'.format(snakemake.threads)
samtools_threads = "" if snakemake.threads <= 1 else '-@ {}'.format(snakemake.threads)

# check arguments.
if len(fqs) > 2:
    sys.exit("get {} files. At most 2 files.".format(len(fqs)))

out_name, out_ext = os.path.splitext(out)
out_ext = out_ext[1:].upper()
if out_ext not in ['SAM', 'BAM']:
    sys.exit("out file extension must be sam|bam, got {}".format(out))

if sort_type not in ['none', 'coordinate', 'queryname']:
    sys.exit('samtools sort must be none, coordinate, queryname. got {}'.format(sort_type))

# sort result
if sort_type == "none":
    out_cmd = f" | {samtools} view {samtools_threads} -h --output-fmt {out_ext} -o {out}"
elif sort_type == 'queryname':
    out_cmd = f" | {samtools} sort {samtools_threads} -n --output-fmt {out_ext} -o {out}"
elif sort_type == 'coordinate':
    out_cmd = f" | {samtools} sort {samtools_threads} --output-fmt {out_ext} -o {out}"
else:
    sys.exit("samtools sort type does not support. got {}".format(sort_type))

# index
if sort_type == "coordinate" and out_ext == "BAM":
    index_cmd = f" && samtools index {samtools_threads} {out}"

# run
shell(
    "{bwa_path} {bwa_cmd} {bwa_threads} {rg} {ref_fa} {fqs}"
    "{out_cmd} {index_cmd}"
)
