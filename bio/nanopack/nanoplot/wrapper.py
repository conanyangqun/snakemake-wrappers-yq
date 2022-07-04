__author__ = "yangqun"
__email__ = "yangqun@jabrehoo.com"
__verison__ = "1.0"

import sys
import os
from snakemake.shell import shell

fastqs = snakemake.input.get('fastqs', '')
bams = snakemake.input.get('bams', '')
input_files = ""

if fastqs and bams:
    # 两者都有
    sys.exit("Either provide fastqs or bams, not both.")
elif not fastqs and not bams:
    # 两者都空
    sys.exit("Either provide fastqs or bams.")
elif fastqs and not bams:
    # fastqs
    input_files = "--fastq " + "--fastq ".join(fastqs)
elif bams and not fastqs:
    # bams
    input_files = "--bam " + "--bam ".join(bams)

out_html = snakemake.output.get('html')
out_prefix = os.path.basename(out_html).replace("NanoPlot-report.html", "")
out_dir = os.path.dirname(out_html)

title = snakemake.params.get('title', '')
title = "" if not title else f"--title {title}"

threads = "" if snakemake.threads <= 1 else "-t {}".format(snakemake.threads)

shell(
    "NanoPlot {threads} -o {out_dir} -p {out_prefix} --tsv_stats --N50 {title} {input_files}"
)
