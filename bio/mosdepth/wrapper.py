__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import sys
import os
from snakemake.shell import shell


# 处理输入参数
assert "bam" in snakemake.input.keys()
input_bam = snakemake.input['bam']

# 处理输出参数
assert "summary" in snakemake.output.keys()
summary_file = snakemake.output['summary']

# 从输出文件路径中提取目录和前缀名
output_dir = os.path.dirname(summary_file)
output_prefix = os.path.basename(summary_file).replace('.mosdepth.summary.txt', '')
full_prefix = os.path.join(output_dir, output_prefix)

# 处理参数
mosdepth = snakemake.params.get('mosdepth', 'mosdepth')
extras = snakemake.params.get('extras', '')
threads = snakemake.threads if snakemake.threads else 1

# 构建线程参数
threads_cmd = f"-t {threads}"

# 处理日志
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# 执行命令
shell(
    "{mosdepth} {threads_cmd} {extras} {full_prefix} {input_bam} {log}"
)
