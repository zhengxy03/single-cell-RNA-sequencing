#!/bin/bash

# 定义样本 ID 数组
samples=(
  "2018jz4"
  "2018jz11"
  "2018jz12"
  "2018jz38"
  "2018jz46"
  "AK658"
  "AK2834"
  "AK3274"
  "AK4295"
)

# 源目录和目标目录的基础路径
source_base="/share/home/wangq/zxy/NSCLC/spatial/sc"
target_base="/share/home/wangq/zxy/NSCLC/spatial/sc/outs"

# 遍历样本 ID 数组
for sample in "${samples[@]}"; do
    # 创建目标文件夹
    mkdir -p "$target_base/$sample"
    # 移动文件
    mv "$source_base/$sample/outs/filtered_feature_bc_matrix/"* "$target_base/$sample"
done




#copy
#!/bin/bash

# 定义要复制的文件夹列表
folders=(
    "SRR15093491"
    "SRR15093535"
    "SRR15093530"
    "SRR15093589"
    "SRR15093585"
    "SRR15093592"
    "SRR15093591"
    "SRR15093537"
    "SRR15093533"
)

# 定义源目录和目标目录
source_dir="/share/home/wangq/zxy/ESCC/matrix/GSE160269"
target_dir="/share/home/wangq/zxy/ESCC/treatment"

# 遍历文件夹列表并复制
for folder in "${folders[@]}"; do
    cp -r "$source_dir/$folder" "$target_dir"
done