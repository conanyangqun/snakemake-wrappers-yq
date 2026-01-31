# split_hp wrapper
__author__ = "yangqun"


import sys
import uuid

from snakemake.shell import shell

s_id = str(uuid.uuid4()) # 避免临时文件冲突
haplotagged_bam = snakemake.input.get('bam')
hp1_bam = snakemake.output.get("hp1_bam")
hp2_bam = snakemake.output.get("hp2_bam")
nohp_bam = snakemake.output.get("nohp_bam", '')

if hp1_bam and hp2_bam:
    pass
else:
    sys.exit("must give hp1_bam and hp2_bam")

all_cmd = ""

hp1_cmd = f"""
samtools view -@ {snakemake.threads} -d HP:1 {haplotagged_bam} | cut -f1 >{s_id}.hp1.qnames.txt \
&& samtools view -@ {snakemake.threads} -b -o {hp1_bam} -N {s_id}.hp1.qnames.txt {haplotagged_bam} \
&& samtools index -@ {snakemake.threads} {hp1_bam} \
&& rm {s_id}.hp1.qnames.txt
"""

hp2_cmd = f"""
samtools view -@ {snakemake.threads} -d HP:2 {haplotagged_bam} | cut -f1 >{s_id}.hp2.qnames.txt \
&& samtools view -@ {snakemake.threads} -b -o {hp2_bam} -N {s_id}.hp2.qnames.txt {haplotagged_bam} \
&& samtools index -@ {snakemake.threads} {hp2_bam} \
&& rm {s_id}.hp2.qnames.txt
"""

nohp_cmd = f"""
samtools view -@ {snakemake.threads} -d HP {haplotagged_bam} | cut -f1 >{s_id}.hp.qnames.txt \
&& samtools view -@ {snakemake.threads} {haplotagged_bam} | cut -f1 >{s_id}.qnames.txt \
&& grep -f {s_id}.hp.qnames.txt -v {s_id}.qnames.txt >{s_id}.nohp.qnames.txt \
&& samtools view -@ {snakemake.threads} -b -o {nohp_bam} -N {s_id}.nohp.qnames.txt {haplotagged_bam} \
&& samtools index -@ {snakemake.threads} {nohp_bam} \
&& rm {s_id}.hp.qnames.txt {s_id}.qnames.txt {s_id}.nohp.qnames.txt
"""

# grep -f use too many memory, if on a large bam. must find a alternative method.

all_cmd = hp1_cmd + '\n' + hp2_cmd + '\n'
if nohp_bam:
    all_cmd += nohp_cmd

shell(
    "{all_cmd}"
)
