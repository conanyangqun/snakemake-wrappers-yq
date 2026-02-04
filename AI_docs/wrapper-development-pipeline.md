## wrapper开发流程

阅读`AI_docs/project-description.md`文件，了解项目结构。

AI模型按照以下流程开发某个wrapper。
- 要开发的wrapper称为目标wrapper。
- 阅读`bio/template`目录下的目录和文件结构，了解wrapper的组织形式。其中`meta.yaml`文件包含了软件的元信息，`env/environment.yaml`文件包含了软件的conda依赖环境，`wrapper.py`文件为程序逻辑，`test/Snakefile`文件为mini-workflow，可用于执行测试用例。
- 阅读`bakta/bakta`的内容，了解良好wrapper的编码形式。
- 阅读目标wrapper目录下的`meta.yaml`文件，了解要开发的wrapper的基本信息。
- 根据目标wrapper的URL，了解软件的详细用法和参数。
- 开发`environment.yaml`文件的内容，注意，要以`major.minor`形式指定软件的版本，不需要安装和测试。
- 开发`wrapper.py`文件的内容，注意，如果软件支持多线程、日志等，需要增加相关参数。
- 在`test/Snakefile`文件中编写mini-workflow，包含测试用例。
- 在`test_wrappers.py`文件中编写测试函数。
- 开发完毕后，等待用户审核代码，不需要执行测试步骤。

AI模型按照以下步骤测试目标wrapper。
- 注意需要将以下命令中的`<wrapper_name>`替换为目标wrapper的名称。
- 如果某命令包含了`conda`，需要单独运行此命令，不要用`&&`与其他命令拼接。
- 在项目根目录下创建`<wrapper_name>_test`目录作为测试wrapper的工作目录，进入该目录。
- 读取wrapper的`environment.yaml`文件，查看wrapper依赖的软件。
- 单独运行此命令，不要与其他命令拼接。用`conda install -c conda-forge -c bioconda -y -p $PWD/conda_env_<wrapper_name> snakemake`创建测试环境并安装`snakemake`。在此命令中补充上wrapper需要的软件，同时安装。
- 单独运行此命令，不要与其他命令拼接。用`conda activate $PWD/conda_env_<wrapper_name>`激活conda环境。
- 执行测试用例，测试wrapper是否可以成功运行。注意，不需要添加`--use-conda`参数。
- 如果mini-workflow本身流程运行错误，查看标准输出定位出错的rule。如果该rule有日志，查看日志文件，查找可能出错的原因。
- 测试成功后，用`conda deactivate`退出conda环境。
- 切换到项目目录下，删除`$PWD/conda_env_<wrapper_name>`目录。
