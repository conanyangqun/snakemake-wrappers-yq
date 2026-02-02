__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import os
from snakemake.shell import shell

# 获取输入文件
fasta = snakemake.input.fasta
db = snakemake.input.db

# 获取输出文件
json_output = snakemake.output.json

# 获取参数
extras = snakemake.params.get("extras", "")
bakta_cmd = snakemake.params.get("bakta", "bakta")

# 从输出文件路径推断输出目录和前缀
output_dir = os.path.dirname(json_output)
prefix = os.path.basename(json_output).replace(".json", "")

# 构建命令行
cmd = [
    bakta_cmd,
    "--db", db,
    "--output", output_dir,
    "--prefix", prefix,
    "--json",
    extras,
    fasta
]

# 执行命令
shell(" ".join(cmd))
