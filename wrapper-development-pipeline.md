## wrapper开发流程

AI模型按照以下流程开发某个wrapper。
- 要开发的wrapper称为目标wrapper。
- 阅读`mosdepth`的内容，了解wrapper的组织形式。
- 阅读目标wrapper目录下的`meta.yaml`文件，了解要开发的wrapper的基本信息。
- 开发`environment.yaml`文件的内容，以`major.minor`形式指定软件的版本，这里不需要安装和测试。
- 开发`wrapper.py`文件的内容，注意，如果软件支持多线程、日志等，需要增加相关参数。
- 在`test/Snakefile`文件中编写测试用例。
- 不需要测试wrapper是否可以成功运行。
- 用`conda activate`激活项目根目录下的`snakemake-wrappers-development`环境。注意，此环境为匿名环境，需要使用绝对路径激活，不能使用单纯的名称。
- 在项目根目录下创建`<wrapper_name>_test`目录，作为测试wrapper的工作目录，进入该目录。
- 测试目标wrapper是否可以成功运行。
- 测试成功后，删除`<wrapper_name>_test`目录。
