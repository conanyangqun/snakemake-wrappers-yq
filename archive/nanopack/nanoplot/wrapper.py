__author__ = "yangqun"

import sys
import os
from snakemake.shell import shell

fastqs = snakemake.input.get('fastqs', '')
bams = snakemake.input.get('bams', '')
out_html = snakemake.output.get('html')
title = snakemake.params.get('title', '')
extras = snakemake.params.get('extras', '')

if not out_html.endswith(".NanoPlot-report.html"):
    sys.exit(
        "html report must be endswith .NanoPlot-report.html"
    )

# 准备输入参数
input_files = ""

if fastqs and bams:
    # 两者都有
    sys.exit("Either provide fastqs or bams, not both.")
elif not fastqs and not bams:
    # 两者都空
    sys.exit("Either provide fastqs or bams.")
elif fastqs and not bams:
    # fastqs
    input_files = "--fastq " + " --fastq ".join(fastqs)
elif bams and not fastqs:
    # bams
    input_files = "--bam " + " --bam ".join(bams)

# 准备输出
out_prefix = os.path.basename(out_html).replace("NanoPlot-report.html", "")
out_dir = os.path.dirname(out_html)

title = "" if not title else f"--title {title}"

threads = "" if snakemake.threads <= 1 else "-t {}".format(snakemake.threads)

shell(
    "NanoPlot {threads} -o {out_dir} -p {out_prefix} {extras} {title} {input_files}"
)
