#!/bin/bash

#Namrata Kalsi
#Pipeline to run variant calling, calculate depth and generate a report for user-specified genes

#USAGE : $0 input.bam gene_list.bed clinvar_file output_directoy_name 

if [[ $# -eq 0 ]] ; then
    echo 'Pipeline to run variant calling, calculate depth and generate a report for user-specified genes'
    echo 'USAGE : ./pipeline_variantCaller.sh input.bam gene_list.bed clinvar_file output_directoy_name'
    exit 0
fi

#Tools used : GATK HaplotypeCaller, GATK VariantRecalibrator, GATK ApplyRecalibratio, bedtools intersect, bedtools genomecov, samtools view 
#scripts involved : cov.py, parse_clnsig.py, draw_script.R
#Python : version 2.7.3 ; modules required : click, pyvcf, logging, gzip
#R : version 3.2.5 ; modules required : ggplot2 
#Make sure the clinvar file has the same naming convention (for eg. chr1 instead of 1) 

#SETTING VARIABLES: 
INPUT_BAM=$1
GENE=$2
BASE="/home/vagrant/nkalsi/ahcg_pipeline"
GATK=$BASE/lib/GenomeAnalysisTK.jar
REF=$BASE/resources/genome/hg19.fa
HAPMAP=$BASE/recalibrator/hapmap_3.3.hg19.sites.vcf.gz
OMNI=$BASE/recalibrator/1000G_omni2.5.hg19.sites.vcf.gz
PHASE=$BASE/recalibrator/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
DBSNP=$BASE/resources/dbsnp/dbsnp_138.hg19.vcf

#MAKING OUTPUT DIRECTORY
mkdir -p $BASE/$4
echo "CALLING VARIANTS"
#HAPLOTYPE CALLER
java -jar $GATK -T HaplotypeCaller -R $REF -I $INPUT_BAM --dbsnp $DBSNP -o $BASE/$4/variants.vcf -nct 1 -gt_mode DISCOVERY

VARIANTS=$BASE/$4/variants.vcf

echo "VARIANT RECALIBRATOR"
#VARIANT RECALIBRATOR
java -Xmx4g -jar $GATK -T VariantRecalibrator -R $REF -input $VARIANTS -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -recalFile $BASE/$4/output.recal -tranchesFile $BASE/$4/output.tranches -rscriptFile $BASE/$4/output.plots.R

RECAL=$BASE/$4/output.recal
TRANCH=$BASE/$4/output.tranches

echo "APPLY RECALIBRATION"
#APPLY RECALIBRATION
java -jar $GATK -T ApplyRecalibration -R $REF -input $VARIANTS -mode SNP --ts_filter_level 99.0 -recalFile $RECAL --tranches_file $TRANCH -o $BASE/$4/variants_recal.vcf

VARIANT_RECAL=$BASE/$4/variants_recal.vcf

echo "CALCULATING COVERAGE" 
#CALCULATING COVERAGE 
samtools view -L $GENE $INPUT_BAM -b > $BASE/$4/genes_of_interest.bam

$BASE/bedtools2/bin/bedtools genomecov -ibam $BASE/$4/genes_of_interest.bam -bga > $BASE/$4/depth_of_coverage.bga.bed

$BASE//bedtools2/bin/bedtools intersect -loj -F 0.10 -a $GENE -b $BASE/$4/depth_of_coverage.bga.bed -bed > $BASE/$4/coverage_compare.bed

awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$6,$7,$4,$8)}' $BASE/$4/coverage_compare.bed > $BASE/$4/depth_of_coverage.final.bed

for gene in $(cut -f4 $GENE | sort -u | xargs)
do
	grep $gene $BASE/$4/depth_of_coverage.final.bed > $BASE/$4/${gene}.txt
	python $BASE/cov.py $BASE/$4/${gene}.txt $BASE/$4/${gene}_coverage.txt 
	#Plotting graph
	./draw_depth.R $BASE/$4/${gene}_coverage.txt $BASE/$4/${gene}_coverage.png
done 

echo "CROSS-REFERNCING WITH CLINVAR"

#CROSS-REFERNCING WITH CLINVAR
#extract gene entries from the clinvar file using genelist
bedtools intersect -a $3 -b $GENE -header > $BASE/$4/clinvar_dcm.vcf
#extract your genes from variants-recal 
bedtools intersect -a $VARIANT_RECAL -b $GENE -header > $BASE/$4/variants_dcm.vcf
#match variants to clinvar 
bedtools intersect -b $BASE/$4/variants_dcm.vcf -a $BASE/$4/clinvar_dcm.vcf -header > $BASE/$4/variants_clinvar.vcf
#report_generation
python ./parse_clnsig.py -i $BASE/$4/variants_clinvar.vcf 2>&1 | tee $BASE/$4/variants_clinvar_temp.txt
name=$(basename $1)
echo "Report for $name" > $BASE/$4/variants_clinvar_table.txt
cut -c 24- $BASE/$4/variants_clinvar_temp.txt  >> $BASE/$4/variants_clinvar_table.txt

convert $BASE/$4/*.png $BASE/$4/all_genes.pdf
pandoc $BASE/$4/variants_clinvar_table.txt -o $BASE/$4/variants_clinvar_table.pdf --no-wrap
pdftk $BASE/$4/variants_clinvar_table.pdf $BASE/$4/all_genes.pdf output $BASE/$4/final_report.pdf

echo "DONE. The final report can be found in $BASE/$4/output final_report.pdf"
#DONE.
