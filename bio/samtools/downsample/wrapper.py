# split_hp wrapper
__author__ = "yangqun"
__email__ = "yangqun@jabrehoo.com"
__verison__ = "1.0.0"

from snakemake.shell import shell

input_bam = snakemake.input.get('bam')
output_bam = snakemake.output.get('bam')
output_bai = snakemake.output.get('bai', '')
extras = snakemake.params.get('extras', '')

threads = 1 if snakemake.threads <= 1 else snakemake.threads

# index
pipe_cmd = ""
if output_bai:
    pipe_cmd = f"&& samtools index -@ {threads} {output_bam}"
else:
    pass

shell(
    "samtools view -@ {threads} -O BAM -o {output_bam} {extras} {input_bam} {pipe_cmd}"
)
