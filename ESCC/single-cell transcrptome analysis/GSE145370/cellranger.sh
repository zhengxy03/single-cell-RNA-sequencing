#!/bin/bash
# 定义样本列表
samples=(
  "SRR11094243"
  "SRR11094244"
  "SRR11094245"
  "SRR11094246"
  "SRR11094247"
  "SRR11094248"
  "SRR11094249"
  "SRR11094250"
  "SRR11094251"
  "SRR11094252"
  "SRR11094253"
  "SRR11094254"
  "SRR11094255"
  "SRR11094256"
  "SRR11094257"
  "SRR11094258"
  "SRR11094259"
)

# 定义公共参数
transcriptome="/share/home/wangq/zxy/ESCC/refgenome/refdata-gex-GRCh38-2024-A"
fastqs="/share/home/wangq/zxy/ESCC/GSE145370"
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