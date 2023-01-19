## HW2
## 2011 Germany E. coli O104:H4 outbreak
### SCAMT Bioinformatics course 

## Preparation

create a directory<p>

<code>mkdir ./HW2_HUS</code>
 
<code>cd ./HW2_HUS</code>

<code>mkdir ./raw</code>

<code>mkdir ./output</code>
 

upload and upzip raw data
 
<code>cd ./raw</code>

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
 
 
## Exploring the dataset
 
SRR292678 - paired end, insert size 470 bp (forward reads, reverse reads, 400 Mb each)<p>
SRR292862 – mate pair, insert size 2 kb (forward reads, reverse reads, 200 Mb each)<p>
SRR292770 – mate pair, insert size 6 kb (forward reads, reverse reads, 200 Mb each)<p>
 
<code>cd ../</code>
 
<code>fastqc ./raw/SRR292678sub_S1_L001_R1_001.fastq ./raw/SRR292678sub_S1_L001_R2_001.fastq  -o ./output </code>
 
<code>fastqc ./raw/SRR292862_S2_L001_R1_001.fastq ./raw/SRR292862_S2_L001_R2_001.fastq  -o ./output</code>
 
<code>fastqc ./raw/SRR292770_S1_L001_R1_001.fastq ./raw/SRR292770_S1_L001_R2_001.fastq  -o ./output</code>
 
According to the FastQC reports per base sequence quality is excellent. No trimming is needed

## K-mer profile and genome size estimation

install jellyfish

<code>conda install -c bioconda kmer-jellyfish</code>

check CPU
 
<code>lscpu</code>
 
run jellyfish

<code>jellyfish count -t 2 -C -m 31 -s 128M -o ./output/SRR292678_31mer ./raw/SRR292678sub_S1_L001_R1_001.fastq  ./raw/SRR292678sub_S1_L001_R2_001.fastq </code>

* <code>-t 2</code> specifies the number of threads to be used <p>
* <code>-C</code> specifies the both strands are considered <p>
* <code>-m 31</code> specifies that now you are counting for 31 mer (i.e., k=31)<p>
* <code>-s 128M</code> is some kind of magical number specification of hash size. This should be as high as the physical memory allows<p>
* <code>-o </code> specifies the prefix of output file names<p>

create histogram 
 
<code>jellyfish histo -o ./output/SRR292678_31mer.histo ./output/SRR292678_31mer</code>

visualize histogram in R and define a peak

<code>spec_31 <- read.table("SRR292678_31mer.histo")</code><p>
<code>plot(spec_31[5:200,],type="l")</code><p>
<code>points(spec_31[16:200,])</code><p>
<code>sum(as.numeric(spec_31[,1]*spec_31[,2]))</code><p>
the total number of k-mer in the distribution is 659921520<p>
<code>spec_31</code><p>
peak is 125 <p>
<code>sum(as.numeric(spec_31[,1]*spec_31[,2]))/125</code><p>
the genome size is 5279372

calculation <p>
<code>N = (M*L)/(L-K+1)</code> <p>
<code>Genome_size = T/N </code><p>
<code>N</code> Depth of coverage <p>
<code>M</code> Kmer peak <p>
<code>K</code> Kmer-size <p>
<code>L</code> avg read length <p>
<code>T</code> Total bases <p>

N = 125 * 90 / (90 - 31 + 1) = 187.5
 
Genome_size = 5499346 * 2 * 90 / 187.5 = 5279372

Use <code>http://qb.cshl.edu/genomescope/</code> to build the graph
 
![image](https://user-images.githubusercontent.com/114799897/213418960-3d59d76a-9046-4499-b167-e47af4669c2b.png)


## Assembling E. coli X genome from paired reads

download SPAdes (read correction and assembly)

<code>https://cab.spbu.ru/software/spades/ </code>

upzip SPAdes

<code>tar -xzf SPAdes-3.15.4-Linux.tar.gz</code>

test SPAdes

<code>./SPAdes-3.15.4-Linux/bin/spades.py --test</code>

run SPAdes for SRR292678

<code>./SPAdes-3.15.4-Linux/bin/spades.py -1 ./raw/SRR292678sub_S1_L001_R1_001.fastq -2 ./raw/SRR292678sub_S1_L001_R2_001.fastq -o ./output/spades</code>

get QUAST (evaluates genome assemblies by computing various metrics)

<code>wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz</code>

unzip QUAST

<code> tar -xzf quast-5.2.0.tar.gz</code>

run QUAST

<code>./quast-5.2.0/quast.py ./output/spades/contigs.fasta ./output/spades/scaffolds.fasta</code>

## Effect of read correction


<code>jellyfish count -t 2 -C -m 31 -s 128M -o ./output/scaffolds_31mer ./output/spades/scaffolds.fasta </code>

<code> jellyfish histo -o ./output/scaffolds_31mer.histo ./output/scaffolds_31mer </code>

Use <code>http://qb.cshl.edu/genomescope/</code> to build the graph

## Impact of reads with large insert size

run SPAdes for SRR292678 as a paired ends, SRR292862 and SRR292770 as a mate pairs

<code>./SPAdes-3.15.4-Linux/bin/spades.py --pe1-1 ./raw/SRR292678sub_S1_L001_R1_001.fastq --pe1-2 ./raw/SRR292678sub_S1_L001_R2_001.fastq --mp1-1 ./raw/SRR292770_S1_L001_R1_001.fastq --mp1-2 ./raw/SRR292770_S1_L001_R2_001.fastq --mp2-1 ./raw/SRR292862_S2_L001_R1_001.fastq --mp2-2 ./raw/SRR292862_S2_L001_R2_001.fastq -o ./output/spades3</code>

run QUAST

<code>./quast-5.2.0/quast.py ./output/spades/contigs.fasta ./output/spades/scaffolds.fasta ./output/spades3/contigs.fasta ./output/spades3/scaffolds.fasta</code>

Save the QUAST report data in your lab journal. In your week report you should provide the
main metrics (N50 and number of contigs) for single-library and three-library assemblies.
Also answer, how did the quality of the assembly improve compared to the previous run of
SPAdes, and why.

## Genome Annotation

get Prokka

<code>conda install -c "bioconda/label/cf201901" prokka</code>

run Prokka
 
<code> prokka --outdir prokka --prefix E.coli_X ./output/spades3/scaffolds.fasta</code>

## Finding the closest relative of E. coli X

install Barrnap (rRNA genes prediction tool)

<code>conda install -c "bioconda/label/cf201901" barrnap</code>

run Barrnap 

<code>barrnap ./output/spades3/scaffolds.fasta --incseq > rrna.fa </code>

go to 
 
 <code>https://blast.ncbi.nlm.nih.gov/Blast.cgi -> nucleotide BLAST</code>

choose:

* RefSeq Genome Database (refseq_genomes)

* Organism fE. coli (taxid:562)

* Enterz Query 1900/01/01:2011/01/01[PDAT]

The closest one is seq NC_011748.1	Escherichia coli 55989

Escherichia coli 55989, complete sequence

<code>https://www.ncbi.nlm.nih.gov/nuccore/NC_011748.1?report=genbank</code>

Donwload FASTA and save it as “55989.fasta”

## What is the genetic cause of HUS?

Install Java

<code>https://www.java.com/en/download/manual.jsp</code>

Install Mauve

<code>https://darlinglab.org/mauve/download.html</code>

Open “Mauve” and select “File” → “Align with progressiveMauve...”. Press “Add sequences” and select the reference genome, then the annotated E. coli X genome (“scaffolds.gbk”), and start the alignment.

Shiga toxin-related genes
stxB (3483605 - 3483874, 270)
stxA (3483886 - 3484845, 960)

![image](https://user-images.githubusercontent.com/114799897/213423412-7962b59e-cea8-4e31-bd57-0cde5cfd5697.png)


## Tracing the source of toxin genes in E. coli X

download annotation

<code> https://disk.yandex.ru/d/aWhOBLVIXR7Oaw </code>

based on this annotation, what is the origin of these toxin genes in E.coli X?

## Antibiotic resistance detection

go to ResFinder
 
<code>https://cge.food.dtu.dk/services/ResFinder/</code>

upload the “scaffolds.fasta” file from the SPAdes output (select “Acquired antimicrobial resistance genes”, select “E. coli”)

results for E. coli 55989
 
<code>https://cge.food.dtu.dk//cgi-bin/webface.fcgi?jobid=63C9164B00001F6BCB7B53E9</code>

results for E. coli X
 
<code>https://cge.food.dtu.dk//cgi-bin/webface.fcgi?jobid=63C9159D00001D302B006644</code>

## Antibiotic resistance mechanism

Search for these enzymes in the same way that we looked for the toxin genes: by using the Sequence Navigator in Mauve. Determine by looking at neighboring genes how E. coli X obtained these genes.

beta-lactamase

bla (4758802-4759935, 1134)
![image](https://user-images.githubusercontent.com/114799897/213423286-2384dec8-3937-4709-a1e6-f915532b9266.png)

