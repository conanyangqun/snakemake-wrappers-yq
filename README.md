# snakemake-wrappers-yq

wrappers是snakemake实现模块化最精细的级别，可在不同流程间重用。官方维护了一个[snakemake wrappers仓库](https://snakemake-wrappers.readthedocs.io/en/stable/)，可看作一个协作项目。

snakemake-wrappers-yq是我仿照官方的仓库实现的一个类似仓库，主要功能如下：
- 实现一些官方没有/不合适的wrapper。
- 更新`test/Snakefile`，支持自动化测试、日常调用。

**特别注意，本仓库优化了官方仓库的部分代码，因此官网文档中的部分功能可能不适用本仓库**。例如：
- ~~本仓库只实现了基本的wrappers，未纳入测试等内容。因此无法像官方一样执行CI测试（后续考虑增加）~~。新开发的wrapper已经加入基于`pytest`的测试和demo数据。
- ~~本仓库每个wrapper的`environment.yaml`为空。这主要是因为我在日常开发中主要用docker/singularity容器打包软件环境，对conda用的较少~~。新开发的wrapper已经加入了基于`conda`的环境配置文件。

**请注意，本项目正处于快速迭代中，尚未发布 (v0.y.z) 版本号**。

## 基本用法

### 作为wrapper使用

与官方wrapper相同，在`Snakefile`中引入wrapper即可。

### 作为workflow使用

官方仓库中的`test/Snakefile`文件展示了wrapper的基本用法，但其本身只能用于固定场景的测试，无法直接用于日常数据分析。对于临时需要执行某个分析的场景，重新编写一个workflow来调用该wrapper显得过于繁琐。为此，本仓库中的`test/Snakefile`进行了优化，同时支持自动化测试、日常调用。

安装snakemake并激活。

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

克隆本仓库到本地某个位置，例如`/Users/yangqun/src/snakemake-wrappers-yq/`下。

```bash
git clone https://github.com/conanyangqun/snakemake-wrappers-yq.git
```

创建测试目录。

```bash
mkdir -p /Users/yangqun/test
cd /Users/yangqun/test
```

每个`test/Snakefile`的参数不同，具体查看源代码中首部。以`Flye`为例，展示如何使用。基本语法如下：

```bash
snakemake \
    -s /Users/yangqun/src/snakemake-wrappers-yq/bio/flye/test/Snakefile \
    --wrapper-prefix /Users/yangqun/src/snakemake-wrappers-yq \
    --cores <n> \
    [--use-conda] \
    [--printshellcmds] \
    [--show-failed-logs] \
    [--config fastq=xxx s_id=xxx data_type=xxx sample_info=xxx out_dir=xxx]
```

分析demo数据。程序会自动下载demo数据，并执行分析，适合快速测试。命令如下：

```bash
snakemake \
    -s /Users/yangqun/src/snakemake-wrappers-yq/bio/flye/test/Snakefile \
    --wrapper-prefix /Users/yangqun/src/snakemake-wrappers-yq \
    --cores 1 \
    --use-conda \
    --printshellcmds \
    --show-failed-logs \
```

分析用户的单个数据。命令如下：

```bash
snakemake \
    -s /Users/yangqun/src/snakemake-wrappers-yq/bio/flye/test/Snakefile \
    --wrapper-prefix /Users/yangqun/src/snakemake-wrappers-yq \
    --cores 1 \
    --use-conda \
    --printshellcmds \
    --show-failed-logs
```

分析用户的多个数据。创建一个tsv文件，包含`s_id, data_type, R1`三列，分别表示样本ID、`Flye`支持的数据类型、fastq文件路径。命令如下：

```bash
snakemake \
    -s /Users/yangqun/src/snakemake-wrappers-yq/bio/flye/test/Snakefile \
    --wrapper-prefix /Users/yangqun/src/snakemake-wrappers-yq \
    --cores 1 \
    --use-conda \
    --printshellcmds \
    --show-failed-logs \
    --config sample_info=xxx.tsv out_dir=xxx
```

TODO: 后续需要为每个wrapper的mini-workflow创建单独的README文件，将mini-workflow的使用方法放入该文件中。、

## wrappers-list

见[文档](https://conanyangqun.github.io/snakemake-wrappers-yq/index.html)。
