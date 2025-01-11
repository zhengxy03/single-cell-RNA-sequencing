#!/bin/bash

#download cellranger
#mkdir apps
#cd apps
#curl -o cellranger-9.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.0.tar.gz?Expires=1736450109&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=lkSVisSxcbwvfnyuUiAy4AN7z1o7Qv6tXZJUVT0s6qiEJTHoleldbO-dKkIZjIa7hGsAW1yzhvVMKMAQd7V~Y0n27Izat6okFL2M~9l4XTr~RU4Bq7nBrlGrGb7PRsAEjroJjuEz42v5apO2tkcT1~~WzvCiBb6BjpxQerhfw0eMrhjvjNjs9lhpoToxbf-72b5wE-nXRYySFQ6M6zIq3e~l6H9BtPAgaHGeVbtj~ettiPCbbSkDsM34rIT5ivIWu9uzmqmqRrn2yFt4pFX3bNO4pnuinyhpqvG5JTOaxZfvISIokWjiChmS3p--MK1qhjpEQN2txDhAN39K~zbMpQ__"
#tar xvfz cellranger-9.0.0.tar.gz
#cd cellranger-9.0.0
#export PATH="$(pwd):$PATH"

#cellranger

#download reference genome
curl -O "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
REFERENCE="/home/zxy0303/project/project/ESCC/refgenome/refdata-gex-GRCh38-2024-A"

#rename for cellranger
FASTQ_DIR="/home/zxy0303/project/project/ESCC/GSE203115"
#<SampleID>_S<SampleNumber>_L<LaneNumber>_R<ReadNumber>_<SetNumber>.fastq.gz
#<SampleID>：样本的唯一标识符（如 SampleA）。
#<SampleNumber>：样本的编号（如 S1）。
#<LaneNumber>：测序 Lane 的编号（如 L001）。
#<ReadNumber>：测序读段的编号（R1 表示正向读段，R2 表示反向读段）。
#<SetNumber>：文件集合编号（通常为 001）。

#GSE203115
sample_number=1

for file in *_1.fastq.gz; do
    sample_id=$(echo "$file" | sed 's/_1\.fastq\.gz//')
    sample_number_formatted=$(printf "S%02d" "$sample_number")
    mv "${sample_id}_1.fastq.gz" "${sample_id}_${sample_number_formatted}_L001_R1_001.fastq.gz"
    mv "${sample_id}_2.fastq.gz" "${sample_id}_${sample_number_formatted}_L001_R2_001.fastq.gz"
    sample_number=$((sample_number + 1))
done

#cellranger count
#one sample
cellranger count --id=SRR19226083 \
                 --transcriptome=/home/zxy0303/project/project/ESCC/refgenome/refdata-gex-GRCh38-2024-A \
                 --fastqs=/home/zxy0303/project/project/ESCC/GSE203115 \
                 --sample=SRR19226083 \
                 --nosecondary \
                 --create-bam=false

#batch process
declare -a samples=(
  "SRR19226085"
  "SRR19226088"
  "SRR19226089"
  "SRR19226084"
  "SRR19226086"
  "SRR19226087"
  "SRR19226090"
  "SRR19226091"
  "SRR19226092"
  "SRR19226093"
  "SRR19226094"
)

for sample_prefix in "${samples[@]}"; do
    echo "now processing $sample_prefix"
    cellranger count --id="$sample_prefix" \
                     --transcriptome=/home/zxy0303/project/project/ESCC/refgenome/refdata-gex-GRCh38-2024-A \
                     --fastqs=/home/zxy0303/project/project/ESCC/GSE203115 \
                     --sample="$sample_prefix" \
                     --nosecondary \
                     --create-bam=false
done

#Integration of results
cat <<EOF > aggregation.csv
SRR19226083,/home/zxy0303/project/project/ESCC/GSE203115/SRR19226083/outs/molecule_info.h5
SRR19226084,
SRR19226085,
SRR19226086,
EOF
