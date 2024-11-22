parallel -k -j 4 "
    hisat2 -p 8 -x ~/project/scRNA/genome/chr1_index/chr1_index \
    -U {1}.fastq.gz \
    -S ../align/chr1_align/{1}.sam 2>../align/chr1_align/{1}.log
" ::: $(ls *.gz | perl -p -e 's/.fastq.gz$//')


parallel -j 4 "
    htseq-count -f bam -s yes {1}.sort.bam ../../annotation/hg.gtf \
        > ../HTseq/{1}.count 2>../HTseq/{1}.log
" ::: $(ls *.sort.bam | perl -p -e 's/\.sort\.bam$//')