#!/bin/bash
#BSUB -J seurat_analysis  # 作业名称
#BSUB -q mpi  # 使用 mpi 队列
#BSUB -n 24  # 请求 24 个核心，必须是 24 的整数倍
#BSUB -o seurat_analysis.out  # 标准输出文件
#BSUB -e seurat_analysis.err  # 标准错误文件

# 运行 R 脚本
Rscript pretreatment.R