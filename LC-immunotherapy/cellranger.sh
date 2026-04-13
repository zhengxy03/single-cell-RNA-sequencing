#!/bin/bash
# 定义样本列表
samples=(
  "2018jz4"
  "2018jz11"
  "2018jz12"
  "2018jz37"
  "2018jz38"
  "2018jz46"
  "AK658"
  "AK2834"
  "AK3274"
  "AK4295"
)

# 定义公共参数
transcriptome="/share/home/wangq/zxy/NSCLC/spatial/sc/refgenome/refdata-gex-GRCh38-2024-A"
fastqs="/share/home/wangq/zxy/NSCLC/spatial/sc"
nosecondary="--nosecondary"
create_bam="--create-bam=false"
localcores=24  # 设置使用的核心数

# 循环处理每个样本
for sample in "${samples[@]}"; do
    job_script="job_${sample}.sh"
    
    cat << EOF > $job_script
#!/bin/bash
#BSUB -q mpi
#BSUB -n $localcores  # 使用指定数量的核心
#BSUB -o ${sample}_%J.out
#BSUB -e ${sample}_%J.err
cellranger count --id=$sample \\
                 --transcriptome=$transcriptome \\
                 --fastqs=$fastqs \\
                 --sample=$sample \\
                 $nosecondary \\
                 $create_bam \\
                 --localcores=$localcores
EOF

    # 提交作业
    bsub < $job_script
done