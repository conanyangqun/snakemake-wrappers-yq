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

assert "db_version_json" in snakemake.input.keys()
db_version_json = os.path.abspath(snakemake.input['db_version_json'])

# 从数据库版本json文件的父目录获取数据库路径
db_path = os.path.dirname(db_version_json)
assert db_path, "Failed to get parent directory from db_version_json file"

# 处理输出参数
assert "json" in snakemake.output.keys()
json_file = snakemake.output['json']

# 从输出文件路径中提取目录和前缀名
output_dir = os.path.dirname(os.path.abspath(json_file))
output_prefix = os.path.basename(json_file).replace('.json', '')

# 处理参数
bakta = snakemake.params.get('bakta', 'bakta')
extras = snakemake.params.get('extras', '')
threads = snakemake.threads if snakemake.threads else 1

# 构建线程参数
threads_cmd = f"--threads {threads}"

# 处理日志
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# 特殊处理
amrfinderplus_db = os.path.join(db_path, "amrfinderplus-db")
pre_cmd = f"amrfinder_update --force_update --database {amrfinderplus_db}"

# 执行命令
shell(
    "({pre_cmd} && "
    "{bakta} --force {threads_cmd} {extras} --db {db_path} --output {output_dir} --prefix {output_prefix} {input_fasta}) "
    "{log}"
)
