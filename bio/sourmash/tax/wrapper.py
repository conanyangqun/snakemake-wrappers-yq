__author__ = "Yang Qun"
__copyright__ = "Copyright 2026, Yang Qun"
__email__ = "yangqunxinli@163.com"
__license__ = "MIT"

import sys
import os
from snakemake.shell import shell

# 处理输入参数
assert "gather_csv" in snakemake.input.keys()
gather_csv = os.path.abspath(snakemake.input['gather_csv'])

assert "taxonomy_csv" in snakemake.input.keys()
taxonomy_csv = snakemake.input['taxonomy_csv']
# 处理taxonomy_csv可能是列表的情况
if isinstance(taxonomy_csv, list):
    taxonomy_csv = " ".join([os.path.abspath(t) for t in taxonomy_csv])
else:
    taxonomy_csv = os.path.abspath(taxonomy_csv)

# 处理输出参数
assert "tax_csv" in snakemake.output.keys()
tax_csv = snakemake.output['tax_csv']

# 从输出文件路径中提取目录和basename
tax_csv_abspath = os.path.abspath(tax_csv)
output_dir = os.path.dirname(tax_csv_abspath)
tax_csv_basename = os.path.basename(tax_csv_abspath)

# 校验tax_csv的basename必须以.with-lineages.csv结尾
if not tax_csv_basename.endswith('.with-lineages.csv'):
    raise ValueError(f"tax_csv basename must end with '.with-lineages.csv', got {tax_csv_basename}")

# 处理输入参数
gather_csv_basename = os.path.basename(gather_csv)

# 校验gather_csv的basename不能与tax_csv的basename（去掉.with-lineages.csv后缀后）相同
tax_csv_basename_no_suffix = tax_csv_basename.replace('.with-lineages.csv', '')
if gather_csv_basename == tax_csv_basename_no_suffix:
    raise ValueError(f"gather_csv basename must not be same as tax_csv basename without .with-lineages.csv suffix. gather_csv basename: {gather_csv_basename}, tax_csv basename without suffix: {tax_csv_basename_no_suffix}")

# 处理参数
sourmash = snakemake.params.get('sourmash', 'sourmash')
extras = snakemake.params.get('extras', '')

# 处理日志
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# 构建临时gather_csv文件路径
temp_gather_csv_basename = f"{tax_csv_basename_no_suffix}.csv"
temp_gather_csv = os.path.join(output_dir, temp_gather_csv_basename)

# 拷贝gather_csv到输出目录
shell(
    "cp {gather_csv} {temp_gather_csv} {log}"
)

# 执行命令
shell(
    "{sourmash} tax annotate {extras} "
    "--gather-csv {temp_gather_csv} "
    "--taxonomy-csv {taxonomy_csv} "
    "--output-dir {output_dir} "
    "{log}"
)
