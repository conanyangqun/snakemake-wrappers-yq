# bwa wrapper
__author__ = "yangqun"
__email__ = "yangqun@jabrehoo.com"
__verison__ = "1.0.0"

import sys
from snakemake.shell import shell

# arguments.
assert 'ref_fa' in snakemake.input.keys()
assert 'fqs' in snakemake.input.keys()
assert 'out' in snakemake.output.keys()

ref_fa = snakemake.input.get('ref_fa')
fqs = snakemake.input.get('fqs')
out = snakemake.output.get('out')

bwa_path = snakemake.params.get('bwa_path', '/usr/gitc/bwa')

default_cmd = "mem -K 100000000 -v 3 -Y"
bwa_cmd = snakemake.params.get('bwa_cmd', default_cmd)

rg = snakemake.params.get('rg')
rg_cmd = f"-R {rg}" if rg else ""

threads = "" if snakemake.threads <= 1 else "-t {}".format(snakemake.threads)

# check arguments.
if len(fqs) > 2:
    sys.exit("get {} fastq files.".format(len(fqs)))

# run
shell(
    "{bwa_path} {bwa_cmd} {threads} {rg_cmd} {ref_fa} {fqs}"
)
