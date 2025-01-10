#!/bin/bash

#download anchr via cargo
#cargo install --git https://github.com/wang-q/anchr --branch main

#anchr help

#fetching data
#data_source:GSE160269

mkdir -p GSE160269
cd GSE160269
cat <<EOF > source.csv

EOF

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml
mlr --itsv --ocsv cat ena_info.tsv > ena_info.csv

mlr --icsv --omd cat ena_info.csv

aria2c -x 4 -s 2 -c -j 2 -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt


#GSE203115
cat <<EOF > source_GSE203115.csv
SRX15289982,,Illumina NovaSeq 6000
SRX15289978,,Illumina NovaSeq 6000
SRX15289977,,Illumina NovaSeq 6000
SRX15289984,,Illumina NovaSeq 6000
SRX15289983,,Illumina NovaSeq 6000
SRX15289981,,Illumina NovaSeq 6000
SRX15289979,,Illumina NovaSeq 6000
SRX15289980,,Illumina NovaSeq 6000
SRX15289976,,Illumina NovaSeq 6000
SRX15289975,,Illumina NovaSeq 6000
SRX15289974,,Illumina NovaSeq 6000
SRX15289973,,Illumina NovaSeq 6000
EOF

#GSE196756
cat <<EOF > source.csv
SRX14189304,,Illumina NovaSeq 6000
SRX14189303,,Illumina NovaSeq 6000
SRX14189302,,Illumina NovaSeq 6000
SRX14189301,,Illumina NovaSeq 6000
SRX14189300,,Illumina NovaSeq 6000
SRX14189299,,Illumina NovaSeq 6000
EOF

#GSE145370
cat <<EOF > source.csv
SRX7733065,,Illumina NovaSeq 6000
SRX7733066,,Illumina NovaSeq 6000
SRX7733067,,Illumina NovaSeq 6000
SRX7733068,,Illumina NovaSeq 6000
SRX7733069,,Illumina NovaSeq 6000
SRX7733069,,Illumina NovaSeq 6000
SRX7733069,,Illumina NovaSeq 6000
SRX7733069,,Illumina NovaSeq 6000
SRX7733069,,Illumina NovaSeq 6000
SRX7733070,,Illumina NovaSeq 6000
SRX7733071,,Illumina NovaSeq 6000
SRX7733072,,Illumina NovaSeq 6000
SRX7733073,,Illumina NovaSeq 6000
SRX7733074,,Illumina NovaSeq 6000
SRX7733075,,Illumina NovaSeq 6000
SRX7733076,,Illumina NovaSeq 6000
SRX7733077,,Illumina NovaSeq 6000
SRX7733078,,Illumina NovaSeq 6000
EOF