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
cd ./raw
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R1_001.fastq.gz
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R2_001.fastq.gz
gunzip SRR292862_S2_L001_R1_001.fastq
gunzip SRR292862_S2_L001_R2_001.fastq
cd ../
fastqc ./raw/SRR292862_S2_L001_R1_001.fastq ./raw/SRR292862_S2_L001_R2_001.fastq  -o ./output
cd ./raw
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R1_001.fastq.gz 
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R2_001.fastq.gz 
gunzip SRR292770_S1_L001_R1_001.fastq
gunzip SRR292770_S1_L001_R2_001.fastq
cd ../
fastqc ./raw/SRR292770_S1_L001_R1_001.fastq ./raw/SRR292770_S1_L001_R2_001.fastq  -o ./output
conda install -c bioconda kmer-jellyfish
N = (M*L)/(L-K+1)
Genome_size = T/N
(N: Depth of coverage, M: Kmer peak, K: Kmer-size, L: avg read length T: Total bases)
You can use k-mer sizes of 31, count kmers with jellyfish count and make a histogram file with jellyfish histo.
