#!/bin/bash

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files
BWA_INDEX="/home/administrator/lifecode/genomes/bwa_hg19_noalt"
REF_GENOME="/home/administrator/lifecode/genomes/bwa_hg19_noalt/hg19.noalt.fa"
REF_GENOME_dict="/home/administrator/lifecode/genomes/Picard_SDict/hg19_noalt/hg19.noalt.fa"
DBSNP_VCF="/home/administrator/lifecode/genomes/dbsnp_hg19/dbsnp_hg19_00-All.vcf.gz"
CLINVAR_VCF="/home/administrator/lifecode/genomes/clnvar_hg19/2025_01_21/clinvar.RENAMED.vcf.gz"
GNOMAD_VCF="/home/administrator/lifecode/genomes/gnomad_hg19/renamed_AF/af-only-gnomad.raw.sites.grch37.CHR.RENAMED_AF.vcf.gz"
COSMIC_VCF="/home/administrator/lifecode/genomes/cosmic_hg19/Cosmic_GenomeScreensMutant_v101_GRCh37.chr.vcf.gz"
TARGETS="/home/administrator/lifecode/genomes/bed_files/WES_HG19/S33266340_Covered.adj.bed"
SNPEFF_JAR="/home/administrator/snpeff/snpEff/snpEff.jar"
HUMANDB="/home/administrator/lifecode/genomes/intervar_humandb/humandb"
INTERVARDB="/home/administrator/lifecode/genomes/intervar_humandb/intervar"

# Function to process a single sample
process_sample() {
	local sample=$1
	local fastq1=$2
	local fastq2=$3

echo "Processing sample: $sample"

# Alignment
#bwa mem -R "@RG\tID:${sample}\tLB:solid_lib\tPL:MGISEQ\tPU:unit1\tSM:${sample}" -t 32 \
#	$REF_GENOME $fastq1 $fastq2 | \
#	samtools view -@ 32 -bS | \
#	samtools sort -@ 32 -o ${sample}_aligned_rg.bam

# Mark Dublicates
#gatk MarkDuplicates \
#	-I ${sample}_aligned_rg.bam \
#	-O ${sample}_aligned_marked.bam \
#	-M ${sample}_output.metrics.txt \
#	--ASSUME_SORT_ORDER coordinate \
#	--CREATE_INDEX true \
#	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

#samtools index ${sample}_aligned_marked.bam

#rm ${sample}_aligned_rg.bam*

###########################################################

# variant calling
#gatk HaplotypeCaller \
#	-R $REF_GENOME_dict \
#	-I ${sample}_aligned_marked.bam \
#	-O ${sample}_variants.vcf.gz \
#	--native-pair-hmm-threads 32

###########################################################

# Extract SNPs
#gatk SelectVariants \
#	-V ${sample}_variants.vcf.gz \
#	-select-type SNP \
#	-O ${sample}_snps.vcf.gz

# Extract INDELs
#gatk SelectVariants \
#	-V ${sample}_variants.vcf.gz \
#	-select-type INDEL \
#	-O ${sample}_indels.vcf.gz

# variant filtering SNPS
#gatk VariantFiltration \
#	-V ${sample}_snps.vcf.gz \
#	-O ${sample}_snps.filtered.vcf.gz \
#	--filter-name "QD_filter" --filter-expression "QD < 2.0" \
#	--filter-name "FS_filter" --filter-expression "FS > 60.0" \
#	--filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
#	--filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
#	--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
#	--filter-name "QUAL_filter" --filter-expression "QUAL < 100.0" \
#	--filter-name "DP_filter" --filter-expression "DP < 10.0" \
#	--filter-name "SOR_filter" --filter-expression "SOR > 3.0"

#rm ${sample}_snps.vcf.gz*

# variant filtering INDELS
#gatk VariantFiltration \
#	-V ${sample}_indels.vcf.gz \
#	-O ${sample}_indels.filtered.vcf.gz \
#	--filter-name "QD_filter" --filter-expression "QD < 2.0" \
#	--filter-name "FS_filter" --filter-expression "FS > 200.0" \
#	--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -20.0" \
#	--filter-name "SOR_filter" --filter-expression "SOR > 10.0"

#rm ${sample}_indels.vcf.gz*

# Merge filtered SNPS and INDELS
#gatk MergeVcfs \
#	-I ${sample}_snps.filtered.vcf.gz \
#	-I ${sample}_indels.filtered.vcf.gz \
#	-O ${sample}_variants.filtered.vcf.gz

#rm ${sample}_snps.filtered.vcf.gz* ${sample}_indels.filtered.vcf.gz*

###########################################################

# Annotate with ClinVar
#echo "Annotate clnvar"

#bcftools annotate --threads 32 -a $CLINVAR_VCF \
#	-c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT,INFO/CLNHGVS,INFO/CLNDISDB,INFO/CLNSIGCONF -O z \
#	-o ${sample}_variants.filtered.clnvar.vcf.gz ${sample}_variants.filtered.vcf.gz

# Index the ClinVar annotated VCF
#bcftools index ${sample}_variants.filtered.clnvar.vcf.gz

#rm ${sample}_variants.filtered.vcf.gz*

###########################################################

# Annotate with dbSNP
#echo "annotate dbsnp"

#bcftools annotate --threads 32 -a $DBSNP_VCF -c ID \
#	${sample}_variants.filtered.clnvar.vcf.gz \
#	-o ${sample}_variants.filtered.clnvar.dbsnp.vcf.gz

# Index the dbSNP annotated VCF
#bcftools index ${sample}_variants.filtered.clnvar.dbsnp.vcf.gz

#rm ${sample}_variants.filtered.clnvar.vcf.gz*

###########################################################

# Annotate with SnpEff
#echo "annotate snpeff" 
#java -jar $SNPEFF_JAR ann -v hg19 -hgvs ${sample}_variants.filtered.clnvar.dbsnp.vcf.gz | \
#bcftools view --threads 32 -Ob -o ${sample}_variants.filtered.clnvar.dbsnp.snpeff.bcf
#bcftools index ${sample}_variants.filtered.clnvar.dbsnp.snpeff.bcf

#rm ${sample}_variants.filtered.clnvar.dbsnp.vcf.gz*

###########################################################

# Annotate with gnomAD
#echo "annotate gnomAD"

#bcftools annotate --threads 32 \
#	-a $GNOMAD_VCF \
#	-c INFO/gnomAD_AF \
#	-O z \
#	-o ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.bcf \
#	${sample}_variants.filtered.clnvar.dbsnp.snpeff.bcf

#bcftools index ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.bcf

#rm ${sample}_variants.filtered.clnvar.dbsnp.snpeff.bcf*

###########################################################

# Annotate with COSMIC

#echo "annotate COSMIC"

#bcftools annotate --threads 32 \
#	-a $COSMIC_VCF \
#	-c INFO/LEGACY_ID \
#	-O z \
#	-o ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.bcf \
#	${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.bcf

#bcftools view ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.bcf > ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.vcf

#rm ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.bcf*

###########################################################

# Annotate Intervar
#echo "annotate Intervar-Assignment of ACGM"

# create link for perl scripts
#ln -s /home/administrator/Annovar/annovar/*.pl .

#Intervar.py -b hg19 \
#	-i ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.vcf --input_type=VCF \
#	-o ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar \
#	-t $INTERVARDB \
#	-d $HUMANDB

#rm *.pl

###########################################################

# Rename intervar chr prefix output
python bio_WES_add_chr_prefix_v0.1.py ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.hg19_multianno.txt.intervar ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.hg19_multianno.txt.chrPrefix.intervar

# Integrate ACGM infor into the annotated vcf
python bio_WES_integrate_ACGM_v0.1.py ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.vcf ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.hg19_multianno.txt.chrPrefix.intervar ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.vcf

###########################################################

echo "Finished processing sample: $sample"

# VEP annotation
#vep -i ${sample}_variants_filtered.vcf.gz -o ${sample}_variants.filtered.vep.vcf.gz --vcf --offline --everything --assembly GRCh37 --fasta $REF_GENOME

}

# Main script
for fastq1 in $FASTQ_DIR/*_1.fq.gz; do
    fastq2="${fastq1/_1.fq.gz/_2.fq.gz}"
    if [ -f "$fastq2" ]; then
        sample=$(basename "$fastq1" | sed 's/_1.fq.gz//')
        process_sample "$sample" "$fastq1" "$fastq2"
    else
        echo "Warning: No matching read 2 file found for $fastq1"
    fi
done

echo "All samples processed"

###########################################################
###########################################################
###########################################################

# Initiate bio_prior script
python bio_WES_prior_v0.1.py ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.vcf ${sample}.xlsx

# Create html report
python bio_WES_report_v0.1.py ${sample}.xlsx ${sample}.html

# convert vcf to bed for creating IGV report
vcf2bed < ${sample}.vcf > ${sample}.bed

# Create IGV report
create_report ${sample}.bed --genome hg19 --flanking 1000 --tracks ${sample}_aligned_marked.bam --output ${sample}.IGV.html

# clean up
#rm ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.hg19_multianno.txt*
#rm ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.vcf*

echo " "
echo " "
echo "PIPELINE COMPLETE! Results saved to ${sample}.xlsx"

