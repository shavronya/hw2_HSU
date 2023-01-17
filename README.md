# hw2
mkdir ./HW2_HUS
cd ./HW2_HUS
mkdir ./raw
mkdir ./output
cd ./raw
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz 
gunzip SRR292678sub_S1_L001_R1_001.fastq
gunzip SRR292678sub_S1_L001_R2_001.fastq
cd ../
fastqc ./raw/SRR292678sub_S1_L001_R1_001.fastq ./raw/SRR292678sub_S1_L001_R2_001.fastq  -o ./output 
