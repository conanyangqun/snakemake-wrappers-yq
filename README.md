# snakemake-wrappers-yq

wrappers是snakemake实现模块化最精细的级别，可实现程序在不同流程间的重用。官方维护了一个[snakemake wrappers仓库](https://snakemake-wrappers.readthedocs.io/en/stable/)，可看作一个协作项目。

snakemake-wrappers-yq是我根据官方的wrappers仓库实现的一个仓库，主要有以下两个目的：
- 实现一些官方没有/不合适的wrapper。
- 学习官方wrapper的实现方式、python测试、打包、发布等软件开发、发布知识。

**特别注意，本仓库并未严格遵守官方wrappers仓库的指南来组织，因此官网文档中的部分功能并不适用本仓库**。例如：
- 本仓库只实现了基本的wrappers，未纳入测试等内容。因此无法像官方一样执行CI测试（后续考虑增加）。
- 本仓库每个wrapper的`environment.yaml`为空。这主要是因为我在日常开发中主要用docker/singularity容器打包软件环境，对conda用的较少。

## wrappers-list

当前已经实现的wrappers列表。

名称|功能|软件|版本
---|---|---|---
bwa|将short reads用bwa mem比对到ref上|bwa|
fastcat|合并fastq文件并生成统计信息|fastcat|
fastqc|对fastq、bam、sam进行QC|fastqc|
hello_world|示例|-|-
kraken2|用kraken2执行物种分类|kraken2|
minimap2|将long reads比对到ref上|minimap2|
mosdepth|对BAM/CRAM文件进行QC|mosdepth|
nanopack/chopper|对fastq数据进行过滤|chopper|
nanopack/nanofilt|对fastq数据进行过滤|nanofilt|
nanopack/nanoplot|对fastq数据进行QC并绘图|nanoplot|
samtools/downsample|对bam文件进行向下抽样|samtools|
samtools/split_hp|根据HP tag将bam中的reads进行拆分｜samtools|
trimmomatic|对Illumina数据进行过滤|trimmomatic|
