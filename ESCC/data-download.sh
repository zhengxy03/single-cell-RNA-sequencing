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