# kraken2  taxonomic sequence classification.
__author__ = "yangqun"
__github__ = "https://github.com/conanyangqun/snakemake-wrappers"

import os
import sys

from snakemake.shell import shell

# required params.
for tag in ['fastq', 'db_hash']:
    if tag not in snakemake.input.keys():
        sys.exit(f"{tag} must be given.")

for tag in ['out', 'report']:
    if tag not in snakemake.output.keys():
        sys.exit(f"{tag} must be given.")

fastq = snakemake.input.get('fastq')
fastq_ext = os.path.splitext(fastq)[-1].lower()
db_hash = snakemake.input.get('db_hash')
db_path = os.path.dirname(os.path.abspath(db_hash))

out = snakemake.output.get('out')
report = snakemake.output.get('report')

# optional params.
kraken2 = snakemake.params.get('kraken2', 'kraken2')
extras = snakemake.params.get('extras', '')
unclassified_out = snakemake.output.get('unclassified_out', '')
classified_out = snakemake.output.get('classified_out', '')

if "--paired" in extras:
    if '#' not in unclassified_out or '#' not in classified_out:
        sys.exit("you specified --paired, but no # sign in unclassified or classified out name.")

cmd_classified = ""
if classified_out:
    cmd_classified += f" --classified_out {classified_out}"
if unclassified_out:
    cmd_classified += f" --unclassified-out {unclassified_out}"

if fastq_ext == 'gz':
    extras += " --gzip-compressed"

cmd_threads = " --threads 1" if snakemake.threads <= 1 else f" --threads {snakemake.threads}"

# run.
shell(
    "{kraken2} {cmd_threads} --db {db_path} {extras} {cmd_classified} "
    "--output {out} --report {report} {fastq}"
)
