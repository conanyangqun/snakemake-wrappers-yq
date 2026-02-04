## wrapper开发流程

AI模型按照以下流程开发某个wrapper。
- 要开发的wrapper称为目标wrapper。
- 阅读`bio/template`目录下的目录结构，了解wrapper的组织形式。其中`meta.yaml`文件包含了软件的元信息，`env/environment.yaml`文件包含了软件的conda依赖环境，`wrapper.py`文件为程序逻辑，`test/Snakefile`文件为mini-workflow，可用于执行测试用例。
- 阅读`bakta/bakta`的内容，了解良好wrapper的编码形式。
- 阅读目标wrapper目录下的`meta.yaml`文件，了解要开发的wrapper的基本信息。
- 根据目标wrapper的URL，了解软件的基本用法。
- 开发`environment.yaml`文件的内容，以`major.minor`形式指定软件的版本，这里不需要安装和测试。
- 开发`wrapper.py`文件的内容，注意，如果软件支持多线程、日志等，需要增加相关参数。
- 在`test/Snakefile`文件中编写mini-workflow，包含测试用例。

AI模型按照以下步骤测试目标wrapper。
- 在项目根目录下创建`<wrapper_name>_test`目录，作为测试wrapper的工作目录，进入该目录。
- 用`conda env create -p $PWD/conda_env_<wrapper_name> -f ../bio/<wrapper_name>/environment.yaml`创建conda环境。将`<wrapper_name>`替换为目标wrapper的名称。
- 用`conda activate $PWD/conda_env_<wrapper_name>`激活conda环境。将`<wrapper_name>`替换为目标wrapper的名称。
- 用`conda install -c conda-forge -c bioconda snakemake`在wrapper的conda环境中安装`snakemake`。
- 执行测试用例，测试wrapper是否可以成功运行。注意，不需要添加`--use-conda`参数。
- 测试成功后，用`conda deactivate`退出conda环境。
- 切换到项目目录下，删除`$PWD/conda_env_<wrapper_name>`目录。将`<wrapper_name>`替换为目标wrapper的名称。

