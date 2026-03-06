# adaptor trimmed
for i in control test
do
 cutadapt -j 20 --discard-untrimmed \
     -a GTTTTAGAGC -g TAAGTAGAGGCTTTATATATCTTGTGGAAAGGACGAAACACC \
     -A GGTGTTTCGT -G CCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC \
     -m 15 --times 2 \
     -o ./${i}_1_trim.fastq.gz -p ./${i}_2_trim.fastq.gz \
     ../1.raw-data/${i}_R1.fq.gz ../1.raw-data/${i}_R2.fq.gz
done

# count analysis
mageck count -l library.csv \
  --fastq control_1_trim.fastq.gz ../2.trim-data/treat_1_trim.fastq.gz \
  --fastq-2 control_2_trim.fastq.gz ../2.trim-data/treat_2_trim.fastq.gz \
  --norm-method control \
  --control-sgrna control_sg.txt \
  --sample-label control,treat \
  -n crispr

# diff test
mageck test -k crispr.count.txt \
  -t treat \
  -c control \
  -n treat_vs_Control \
  --norm-method control \
  --control-sgrna control_sg.txt \
  --gene-lfc-method alphamedian \
  --output-prefix difftest
