#!/bin/bash

# 定义样本 ID 数组
samples=(
    "SRR11094242"
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

# 源目录和目标目录的基础路径
source_base="/share/home/wangq/zxy/ESCC/GSE145370"
target_base="/share/home/wangq/zxy/ESCC/matrix/GSE145370"

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