__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import sys
import os
from snakemake.shell import shell

# 处理输入参数
assert "contigs" in snakemake.input.keys()
input_contigs = snakemake.input['contigs']

# 处理多个输入文件
if isinstance(input_contigs, list):
    input_contigs_str = " ".join(input_contigs)
else:
    input_contigs_str = input_contigs

# 处理输出参数
assert "report" in snakemake.output.keys()
output_report = snakemake.output['report']

# 检查输出文件格式
output_ext = os.path.splitext(output_report)[1].lower()
assert output_ext in ['.tsv', '.csv'], f"Output report must be .tsv or .csv format, got {output_ext}"

# 处理json输出
json_output = snakemake.output.get('json', None)
json_cmd = f"--json {json_output}" if json_output else ""

# 处理参数
mlst = snakemake.params.get('mlst', 'mlst')
extras = snakemake.params.get('extras', '')
threads = snakemake.threads if snakemake.threads else 1

# 构建线程参数
threads_cmd = f"--threads {threads}"

# 构建输出格式参数
if output_ext == '.csv':
    format_cmd = "--csv"
else:
    format_cmd = ""

# 处理日志
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# 执行命令
shell(
    "{mlst} {threads_cmd} {format_cmd} {json_cmd} {extras} {input_contigs_str} >{output_report} "
    "{log}"
)
