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
N = 125 * 
You can use k-mer sizes of 31, count kmers with jellyfish count and make a histogram file with jellyfish histo.
https://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html

 jellyfish count -t 2 -C -m 31 -s 10M -o ./output/SRR292678_31mer ./raw/SRR292678sub_S1_L001_R1_001.fastq  ./raw/SRR292678sub_S1_L001_R2_001.fastq 
 
-t 2
specifies the number of threads to be used. 
-C
specifies the both strands are considered. 
-m 31
specified that now you are counting for 31 mer (i.e., k=31)
-s 10M
is some kind of magical number specification of hash size. This should be as high as the physical memory allows. The higher the faster, but exceeding the available memory leads to failure or extremely slow counting.

jellyfish histo -o ./output/SRR292678_31mer.histo ./output/SRR292678_31mer 

spec_31 <- read.table("SRR292678_31mer.histo")

plot(spec_31[5:200,],type="l")

points(spec_31[16:200,])

sum(as.numeric(spec_31[,1]*spec_31[,2]))

#659921520

spec_31

#125 peak

sum(as.numeric(spec_31[,1]*spec_31[,2]))/125

#5279372

tar -xzf SPAdes-3.15.4-Linux.tar.gz 

./SPAdes-3.15.4-Linux/bin/spades.py --test 

./SPAdes-3.15.4-Linux/bin/spades.py -1 ./raw/SRR292678sub_S1_L001_R1_001.fastq -2 ./raw/SRR292678sub_S1_L001_R2_001.fastq -o ./output/spades

 wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz
 
 tar -xzf quast-5.2.0.tar.gz
 
 cd quast-5.2.0

./quast-5.2.0/quast.py ./output/spades/contigs.fasta ./output/spades/scaffolds.fasta

./SPAdes-3.15.4-Linux/bin/spades.py --pe1-1 ./raw/SRR292678sub_S1_L001_R1_001.fastq --pe1-2 ./raw/SRR292678sub_S1_L001_R2_001.fastq --mp1-1 ./raw/SRR292770_S1_L001_R1_001.fastq --mp1-2 ./raw/SRR292770_S1_L001_R2_001.fastq --mp2-1 ./raw/SRR292862_S2_L001_R1_001.fastq --mp2-2 ./raw/SRR292862_S2_L001_R2_001.fastq -o ./output/spades2 

./quast-5.2.0/quast.py ./output/spades/contigs.fasta ./output/spades/scaffolds.fasta ./output/spades2/contigs.fasta ./output/spades2/scaffolds.fasta

conda install -c "bioconda/label/cf201901" prokka

 prokka --outdir prokka --prefix E.coli_X ./output/spades2/scaffolds.fasta  --centre X --compliant 

conda install -c "bioconda/label/cf201901" barrnap

barrnap rrna.fa < ./output/spades2/scaffolds.fasta > rrna.gff 

https://www.java.com/en/download/manual.jsp

https://darlinglab.org/mauve/download.html

barrnap ./output/spades/scaffolds.fasta --incseq > rrna.fa

https://blast.ncbi.nlm.nih.gov/Blast.cgi

Filed database Ref Seq Genome Database (refseq_genomes)

Organism filed E. coli (taxid:562)

Enterz Query field 1900/01/01:2011/01/01[PDAT]

https://blast.ncbi.nlm.nih.gov/Blast.cgi#alnHdr_218693476

Escherichia coli 55989



