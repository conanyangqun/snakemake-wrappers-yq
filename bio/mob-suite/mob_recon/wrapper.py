__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import sys
import os
from snakemake.shell import shell


# 处理输入参数
assert "fasta" in snakemake.input.keys()
input_fasta = snakemake.input['fasta']

# 处理输出参数
assert "contig_report" in snakemake.output.keys()
contig_report = snakemake.output['contig_report']

# 从输出文件路径中提取目录
output_dir = os.path.dirname(os.path.abspath(contig_report))

# 处理参数
mob_recon = snakemake.params.get('mob_recon', 'mob_recon')
extras = snakemake.params.get('extras', '')
threads = snakemake.threads if snakemake.threads else 1

# 构建线程参数
threads_cmd = f"--num_threads {threads}"

# 处理日志
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# 执行命令
shell(
    "{mob_recon} --force {threads_cmd} {extras} -i {input_fasta} -o {output_dir} {log}"
)
