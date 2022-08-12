# trimmomatic wrapper
__author__ = "yangqun"
__email__ = "yangqun@jabrehoo.com"
__verison__ = "1.0.0"

import sys
import os
from snakemake.shell import shell

# arguments.
for tag in ['data', 'trimmomatic', 'adapters']:
    assert tag in snakemake.input.keys()

data = snakemake.input.get('data')
trimmomatic = snakemake.input.get('trimmomatic')
adapters = snakemake.input.get('adapters')

java = snakemake.input.get('java', 'java')
extras = snakemake.params.get('extras', '')

assert 'outfiles' in snakemake.output.keys()
outfiles = snakemake.output.get('outfiles')

# PE or SE
if len(data) not in [1, 2]:
    sys.exit("input data must be SE or PE. got {} files.".format(len(data)))
elif len(data) == 1:
    if len(outfiles) != 1:
        sys.exit("Got SE data, but outfiles is no 1.")
    else:
        mode = "SE"
elif len(data) == 2:
    if len(outfiles) != 4:
        sys.exit("Got PE data, but outfiles is no 4.")
    else:
        mode = "PE"

shell(
    "{java} -jar {trimmomatic} {mode} {data} {outfiles} {extras} "
    "ILLUMINACLIP:{adapters}:2:30:10:2:True SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36"
)
