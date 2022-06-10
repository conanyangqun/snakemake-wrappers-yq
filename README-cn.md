# snakemake-wrappers

wrappers是snakemake中重用程序最精细的单元，是实现程序重用最基本的形式。官方维护了一个[wrappers仓库](https://github.com/snakemake/snakemake-wrappers)，实现了一些常用软件的包装。

本仓库存在的目的，在于实现一些官方没有/不合适的包装器。**特别注意，目前本仓库没有严格遵守官方的wrappers仓库来组织，因此可能有些功能并不适用**。

## 向仓库贡献代码的方式

*此处指往官方wrappers仓库贡献代码的形式*。

基本流程为：

-   fork原仓库为自己的仓库。、
-   clone仓库到本地。
-   创建新的分支。
-   提交改动到新的分支，并提交到远程仓库。
-   创建一个PR。

代码的编写风格如下：

-   提供`meta.yaml`文件，描述包装器。
-   提供`environment.yaml`列出所有需要的软件包，遵守[最佳实践](https://stackoverflow.com/a/64594513/2352071)。

-   添加`wrapper.py`或`wrapper.R`文件，可以处理任意的输入/输出路径。

-   在`test`文件夹提供小型测试用例，提供一个实例`Snakefile`展示如何使用wrapper，rule名字应该具有描述性且以[snake_case](https://en.wikipedia.org/wiki/Snake_case)形式编写，一些小型测试数据，在`test.py`中实现调用。
-   确保python文件具有一致的格式，snakefiles lint。

更多信息参考：[贡献代码](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html)。

