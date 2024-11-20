# single-cell-RNA-sequencing
# 1 download data
## 1.1 experiment data
GSE144240
```
mkdir -p ~/project/scRNA
cd ~/project/scRNA
mkdir genome sequence
cd sequence
nohup prefetch SRR11832836 SRR11832837 -O . &

parallel -j 4 "
  fastq-dump --split-files --gzip {1}
" ::: $(ls *.sra)

rm *.sra
```
## 1.2 ref data
ensembl-human
```
cd ../genome
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```
## 1.3 annotation data
```
aria2c -d ./ -Z https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
```
# 2 quality control and trimming
* fastqc
```
cd ../sequence
fastqc --threads 3 *.fastq.gz
```
* trimmomaric
```
parallel -j 4 "
  java -jar ~/biosoft/Trimmomatic-0.39/trimmomatic-0.39.jar \
  PE -phred33 {1} ../trim/{1} \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:5:15 MINLEN:30 \
" ::: $(ls *.gz)
```
# 3 seq alignment
* hisat2
```
cd ~/project/rat/genome
mkdir index
cd index

hisat2-build -p 6 ../rn7.chr1.fa rn7.chr1
