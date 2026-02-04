__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import sys
import os
from snakemake.shell import shell

# 处理输入参数
assert "query" in snakemake.input.keys()
query = snakemake.input['query']

assert "db" in snakemake.input.keys()
db = snakemake.input['db']

# 处理输出参数
assert "csv" in snakemake.output.keys()
csv_file = snakemake.output['csv']

# 处理参数
sourmash = snakemake.params.get('sourmash', 'sourmash')
extras = snakemake.params.get('extras', '')
# threads = snakemake.threads if snakemake.threads else 1

# 构建线程参数
# threads_cmd = f"--threads {threads}"

# 处理日志
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# 执行命令
shell(
    "{sourmash} gather {extras} {query} {db} --output {csv_file} "
    "{log}"
)