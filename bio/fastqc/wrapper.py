# fastqc wrapper
__author__ = "yangqun"
__email__ = "yangqun@jabrehoo.com"
__verison__ = "1.0.0"

import sys
import os
from snakemake.shell import shell


# arguments
assert "data" in snakemake.input.keys()
input_data = snakemake.input['data']
fastqc = snakemake.input.get('fastqc', 'fastqc')
java = snakemake.input.get('java', '')
outdir = snakemake.params.get('outdir', '')
extras = snakemake.params.get('extras', '')

# program
java_path = f"-j {java}" if java else ""
out_cmd = f"-o {outdir}" if outdir else ""

shell(
    "{fastqc} {java_path} {out_cmd} {extras} {input_data}"
)
