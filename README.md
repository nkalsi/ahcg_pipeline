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


Verifying the variant calls from the pipeline with the GIAB variant calls for NA12878:

1. Generating a bed file containing the coordinates (exonic level) for all the genes from color genomics and otogenetics, breast and ovarian cancer panel
The gene list can be found in br_ov/gene_list.txt
Getting all nmids from the breast and ovarian cancer gene list.
```{sh}
awk '{print "\\<" $2 "\\>" }' breastcancer_genes.txt > genelist.txt
```
Getting all nmids from the breast and ovarian cancer gene list.
```{sh}
grep -f gene_nmids.txt hg19_refGene.txt > genes_bed.txt
```
Converting to bed file and adding 20 bp flanks on both ends.
```{sh}
./bedconverter.py -i genes_bed.txt -o ovbr.bed
```
Remove the last 4 lines for the PMS2 gene from the BED file, since exons 12-15 aren't analyzed.

2.Extract variants from your vcf files, falling within the regions listed in the bed file for breast and ovarian cancer panel.
```{sh}
bedtools intersect -header -a giab-with_chr.vcf -b ov_br.bed > giab-vcf_intersect.vcf
```

3. Determine the number of overlapping calls with GIAB variant calls for NA12878.
The GIAB vcf file has a different chromosomal notation.
```{sh}
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' no_chr.vcf > with_chr.vcf 
```
Extract variants from the GIAB vcf file that fall within the regions listed in the bed file for breast and ovarian cancer panel.
```{sh}
bedtools intersect -header -a giab-with_chr.vcf -b ov_br.bed > giab-vcf_intersect.vcf
```
Calculating the number of overlaps
```{sh}
bedtools intersect -header -a giab-vcf_intersect.vcf -b vcf_intersect.vcf > final-output.vcf
```



#Course : Applied Human Computational Genomics (BIOL8803F)
