__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import sys
import os
from snakemake.shell import shell

# 处理输入参数
assert "data" in snakemake.input.keys()
input_data = snakemake.input['data']

# 处理输出参数
assert "sig" in snakemake.output.keys()
sig_file = snakemake.output['sig']

# 处理参数
sourmash = snakemake.params.get('sourmash', 'sourmash')
subcmd = snakemake.params.get('subcmd', 'dna')
extras = snakemake.params.get('extras', '')
# threads = snakemake.threads if snakemake.threads else 1

# 构建线程参数
# threads_cmd = f"--threads {threads}"

# 处理日志
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# 执行命令
shell(
    "{sourmash} sketch {subcmd} {extras} --output {sig_file} {input_data} "
    "{log}"
)