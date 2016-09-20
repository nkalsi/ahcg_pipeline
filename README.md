# ahcg_pipeline
Variant calling pipeline for genomic data analysis

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
To get the pipeline
```{sh}
git clone https://github.com/shashidhar22/ahcg_pipeline
```
Downloading and installing modules 
```{sh}
git pull origin master
```

Downloading test data: 
```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```
Gunzipping downloaded NIST files
```{sh}
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```
Downloading the reference genome 
```{sh}
- wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz
- tar -xvzf resources.tar.gz 
```

Building the bowtie index
```{sh}
./lib/bowtie2-2.2.9/bowtie2-build ./resources/genome/hg19.fa hg19
```
Installing samtools and reating fasta index file using samtools faidx
sudo apt-get install samtools
```{sh}
- samtools faidx hg19.fa
```
Installing java version 1.8
```{sh}
sudo apt-get install software-properties-common python-software-properties
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
sudo apt-get install oracle-java8-installer
```
Building sequence dictionaries
```{sh}
java -jar ./lib/picard.jar CreateSequenceDictionary R=./resources/genome/hg19.fa O=hg19-2.dict
```

Running the pipeline
```{sh}
python3 ./ahcg_pipeline.py -t /home/vagrant/nkalsi/ahcg_pipeline/lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b /home/vagrant/nkalsi/ahcg_pipeline/lib/bowtie2-2.2.9/bowtie2 -p /home/vagrant/nkalsi/ahcg_pipeline/lib/picard.jar -g /home/vagrant/nkalsi/ahcg_pipeline/lib/GenomeAnalysisTK.jar -i /home/vagrant/nkalsi/ahcg_pipeline/test*.fastq -w /home/vagrant/nkalsi/ahcg_pipeline/hg19 -d /home/vagrant/nkalsi/ahcg_pipeline/resources/dbsnp/dbsnp_138.hg19.vcf -r /home/vagrant/nkalsi/ahcg_pipeline/resources/genome/hg19.fa -a /home/vagrant/nkalsi/ahcg_pipeline/lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa -o out_dir 
```

Change the remote url for the GIT repository 
```{sh}
vi .git/config
url = https://github.com/nkalsi/ahcg_pipeline
```
Pushing files to github
```{sh}
git add README
git commit -m "readme"
git push -u origin master
```
Adding files to .gitignore
```{sh}
vi .gitignore
#enter the directory name in .gitignore,then 
git add .gitignore
git commit -m "updated"
git push origin master
```
Updating README in github
```{sh}
git add README
git commit -m "updated"
git push origin master
```

Retrieving the sequence for BRCA1
```{sh}
wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
#extract the six BRCA1 entries from the file 
grep 'BRCA1' hg19_refGene.txt > BRCA1.txt
grep 'NM_007294' BRCA1.txt > BRCA1.txt
```

Converting to bed file
```{sh}
./convert_to_bed.sh BRCA1.txt BRCA1.bed
```
```{sh}
gene	ncbi reference	OMIM
BRCA1	NM_007298.3	113705
BRCA2	NM_000059.3	600185
MLH1	NM_000249.3	120436
MSH2	NM_000251.2	609309
MSH6	NM_000179.2	600678
PMS2	NM_000535.5  	
EPCAM	NM_002354.2
TP53	NM_000546.5	191170
PTEN	NM_000314.4  	601728
STK11	NM_000455.4	602216
CDH1	NM_004360.3	192090
PALB2	NM_024675.3
CHEK2	NM_001005735.1
ATM	NM_000051.3	607585
NBN	NM_002485.4	602667
BARD1	NM_000465.3	601593
BRIP1	NM_032043.2	605882
RAD51C	NM_002876.3
RAD51D	NM_001142571.1
AR	NM_000044.3	313700
CASP8	NM_001080124.1	601763
CHEK2 	NM_001005735.1	604373
DIRAS3	NM_004675.2	605193
ERBB2	NM_001005862.1	164870
PALB2	NM_024675.3	601355
RAD50	NM_005732.3	604040
RAD51A	NM_001164269.1	179617  
TGFB1	NM_000660.4	190180
```


#Course : Applied Human Computational Genomics (BIOL8803F)
