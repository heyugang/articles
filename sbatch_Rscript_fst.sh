#!/bin/bash
#SBATCH --job-name=my_R_job       # 作业名称
#SBATCH --ntasks=1                # 运行一个任务
#SBATCH --cpus-per-task=1         # 每个任务一个CPU
#SBATCH --time=01:00:00           # 最大运行时间（hh:mm:ss）
#SBATCH --mem=1G                  # 使用的内存量
#SBATCH --output=my_R_job.out     # 标准输出和错误的文件名
##sbatch -N 1 -c 4 --output=fst_R_job.out --time=24:00:00 -J fst_class sbatch_Rscript.sh

#####################################################
# 这里加载 Conda 环境的命令应该和你在命令行中输入的一样
source /public/home/06025/WORK/Software/anaconda3_202205/etc/profile.d/conda.sh  # 替换成你conda安装的实际路径，这是学校服务器上的路径
conda activate lianghui  # 激活你的 Conda 环境

# 运行 R 脚本
Rscript /public/home/06025/WORK/lianghui/sh_file/fst_loop_10_script.R  # 替换成你的 R 脚本的实际路径
