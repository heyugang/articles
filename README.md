# 代码说明

---

该部分包含文章所用核心代码和代码提交脚本。

## 核心代码

文章核心代码，包括SNP选择和品种分类：

- **R代码提交后台运行的代码文件**
    - sbatch_Rscript_fst.sh，FST方法对应的代码提交文件
    - sbatch_Rscript_relf.sh，Relief-F方法对应的代码提交文件
    - sbatch_Rscript_mrmr.sh，mrmr方法对应的代码提交文件

\
- **在服务器的conda环境下，运行R的R代码文件**
    - fst_loop_script.R，FST方法运行5个训练集的R代码文件
    - relf_loop_script.R，Relief-F方法运行第 1 个训练集的R代码文件  
    - mrmr_loop_script.R，运行第 1 个训练集
