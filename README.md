HW2
2011 Germany E. coli O104:H4 outbreak
SCAMT Bioinformatics course 

##Preparation

###create a directory and folders<p>
<code>mkdir ./HW2_HUS</code>
 
<code>cd ./HW2_HUS</code>

<code>mkdir ./raw</code>

<code>mkdir ./output</code>

<code>cd ./raw</code>
 

###upload and upzip data

<code>wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz</code>
 
<code>wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz </code>
 
<code>gunzip SRR292678sub_S1_L001_R1_001.fastq</code>
 
<code>gunzip SRR292678sub_S1_L001_R2_001.fastq</code>

<code>wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R1_001.fastq.gz</code>
 
<code>wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R2_001.fastq.gz</code>
 
<code>gunzip SRR292862_S2_L001_R1_001.fastq</code>
 
<code>gunzip SRR292862_S2_L001_R2_001.fastq</code>
 
<code>wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R1_001.fastq.gz </code>
 
<code>wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R2_001.fastq.gz </code>
 
<code>gunzip SRR292770_S1_L001_R1_001.fastq</code>
 
<code>gunzip SRR292770_S1_L001_R2_001.fastq</code>
 
 
##Exploring the dataset
 
SRR292678 - paired end, insert size 470 bp (forward reads, reverse reads, 400 Mb each)<p>
SRR292862 – mate pair, insert size 2 kb (forward reads, reverse reads, 200 Mb each)<p>
SRR292770 – mate pair, insert size 6 kb (forward reads, reverse reads, 200 Mb each)<p>
 
<code>cd ../</code>
 
<code>fastqc ./raw/SRR292678sub_S1_L001_R1_001.fastq ./raw/SRR292678sub_S1_L001_R2_001.fastq  -o ./output </code>
 
<code>fastqc ./raw/SRR292862_S2_L001_R1_001.fastq ./raw/SRR292862_S2_L001_R2_001.fastq  -o ./output</code>
 
<code>fastqc ./raw/SRR292770_S1_L001_R1_001.fastq ./raw/SRR292770_S1_L001_R2_001.fastq  -o ./output</code>
 
According to the FastQC reports per base sequence quality is excellent. No trimming is needed

##K-mer profile and genome size estimation



<code></code>
<code></code>
<code></code>
<code></code>


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

https://www.ncbi.nlm.nih.gov/nuccore/NC_011748.1?report=genbank

Escherichia coli 55989

Java installation https://www.java.com/en/download/manual.jsp

To compare the E. coli X with the reference genome, first install Mauve
(http://darlinglab.org/mauve/download.html) on your computer. Open “Mauve” and
select “File” → “Align with progressiveMauve...”. Press “Add sequences” and select the
reference genome, then the annotated E. coli X genome (“scaffolds.gbk”), and start the
alignment.

stxB stxA

https://cge.food.dtu.dk/services/ResFinder/

https://cge.food.dtu.dk//cgi-bin/webface.fcgi?jobid=63C8471E00005413E38574C2

bla

jellyfish count -t 2 -C -m 31 -s 10M -o ./output/31mer ./output/spades2/scaffolds.fasta 

jellyfish histo -o ./output/31mer.histo ./output/31mer
