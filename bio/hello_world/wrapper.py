__author__ = "yangqun"
__email__ = "yangqun@jabrehoo.com"

import sys
from snakemake.shell import shell

name = snakemake.params.get('name', '')
name = name if name else "world"

shell(
    "echo \"Hello {name}\" >{snakemake.output}"
)
