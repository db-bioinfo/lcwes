#!/bin/bash

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files
THREADS=32
REF_GENOME="/home/administrator/lifecode/genomes/databases/bwa_hg19/hg19.fa"
TARGETS="/home/administrator/lifecode/genomes/bed_files/WES_HG19/S33266340_Covered.adj.bed"
DBSNP="/home/administrator/lifecode/genomes/databases/dbsnp_hg19/dbsnp_hg19_00-All.vcf.gz"
CLNVAR="/home/administrator/lifecode/genomes/databases/clnvar_hg19/clinvar.chr.vcf.gz"
COSMIC="/home/administrator/lifecode/genomes/databases/cosmic_hg19/Cosmic_GenomeScreensMutant_Normal_v101_GRCh37.chr.vcf.gz"
SNPEFF="/home/administrator/SNPEFF/snpEff/snpEff.jar"
INTERVARDB="/home/administrator/lifecode/genomes/intervar_humandb/intervar"
HUMANDB="/home/administrator/lifecode/genomes/intervar_humandb/humandb"
FREEBAYES_REGIONS="/home/administrator/lifecode/genomes/freebayes_hg19/hg19_regions.txt"

# Function to process a single sample
process_sample() {
	local sample=$1
	local fastq1=$2
	local fastq2=$3

echo "Processing sample: $sample"

###########################################################

# Alignment
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running alignment..."
bwa mem -R "@RG\tID:${sample}\tLB:exome_lib\tPL:MGISEQ\tPU:unit1\tSM:${sample}" -t $THREADS \
	$REF_GENOME $fastq1 $fastq2 | \
	samtools view -@ $THREADS -bS | \
	samtools sort -@ $THREADS -o ${sample}_aligned_rg.bam

# Mark Dublicates
echo "[$(date '+%Y-%m-%d %H:%M:%S')] MarkDuplicates..."
gatk MarkDuplicates \
	-I ${sample}_aligned_rg.bam \
	-O ${sample}_aligned_marked.bam \
	-M ${sample}_output.metrics.txt \
	--ASSUME_SORT_ORDER coordinate \
	--CREATE_INDEX true \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

samtools index ${sample}_aligned_marked.bam

rm ${sample}_aligned_rg.bam*

###########################################################

# Variant calling GATK
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Callng Variants..."
gatk HaplotypeCaller \
	-R $REF_GENOME \
	-I ${sample}_aligned_marked.bam \
	-O ${sample}_variants.vcf.gz \
	--native-pair-hmm-threads $THREADS

###########################################################

# Extract SNPs
gatk SelectVariants \
	-V ${sample}_variants.vcf.gz \
	-select-type SNP \
	-O ${sample}_snps.vcf.gz

# Extract INDELs
gatk SelectVariants \
	-V ${sample}_variants.vcf.gz \
	-select-type INDEL \
	-O ${sample}_indels.vcf.gz

# variant filtering SNPS
gatk VariantFiltration \
	-V ${sample}_snps.vcf.gz \
	-O ${sample}_snps.filtered.vcf.gz \
	--filter-name "QD_filter" --filter-expression "QD < 2.0" \
	--filter-name "FS_filter" --filter-expression "FS > 60.0" \
	--filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
	--filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
	--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "QUAL_filter" --filter-expression "QUAL < 100.0" \
	--filter-name "DP_filter" --filter-expression "DP < 10.0" \
	--filter-name "SOR_filter" --filter-expression "SOR > 3.0"

rm ${sample}_snps.vcf.gz*

# variant filtering INDELS
gatk VariantFiltration \
	-V ${sample}_indels.vcf.gz \
	-O ${sample}_indels.filtered.vcf.gz \
	--filter-name "QD_filter" --filter-expression "QD < 2.0" \
	--filter-name "FS_filter" --filter-expression "FS > 200.0" \
	--filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -20.0" \
	--filter-name "SOR_filter" --filter-expression "SOR > 10.0"

rm ${sample}_indels.vcf.gz*

# Merge filtered SNPS and INDELS
gatk MergeVcfs \
	-I ${sample}_snps.filtered.vcf.gz \
	-I ${sample}_indels.filtered.vcf.gz \
	-O ${sample}_GATK.filtered.vcf.gz

rm ${sample}_snps.filtered.vcf.gz* ${sample}_indels.filtered.vcf.gz*

###########################################################

# Annotate with dbSNP
echo "annotate dbsnp"
bcftools annotate --threads 32 \
	-a $DBSNP \
	-c ID \
	-O z \
	${sample}_GATK.filtered.vcf.gz \
	-o ${sample}_GATK.filtered.dbsnp.vcf.gz

# Index the dbSNP annotated VCF
bcftools index ${sample}_GATK.filtered.dbsnp.vcf.gz

rm ${sample}_GATK.filtered.vcf.gz*

###########################################################

# Annotate with ClinVar
echo "Annotate ClinVar"
bcftools annotate --threads 32 -a $CLNVAR \
	-c INFO/CLNSIG,INFO/CLNDN,INFO/CLNHGVS,INFO/CLNSIGCONF -O z \
	-o ${sample}_GATK.filtered.dbsnp.clnvar.vcf.gz ${sample}_GATK.filtered.dbsnp.vcf.gz

# Index the ClinVar annotated VCF
bcftools index ${sample}_GATK.filtered.dbsnp.clnvar.vcf.gz

rm ${sample}_GATK.filtered.dbsnp.vcf.gz*

###########################################################

# Annotate with COSMIC
echo "annotate COSMIC"
bcftools annotate --threads 32 \
	-a $COSMIC \
	-c INFO/LEGACY_ID \
	-O z \
	-o ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.vcf.gz \
	${sample}_GATK.filtered.dbsnp.clnvar.vcf.gz

bcftools index ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.vcf.gz

rm ${sample}_GATK.filtered.dbsnp.clnvar.vcf.gz*

###########################################################

# Annotate with SnpEff
echo "annotate snpeff"
java -jar $SNPEFF ann -v hg19 -hgvs -canon -lof ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.vcf.gz | \
bcftools view --threads 32 -Ob -o ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz
bcftools index ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz

rm ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.vcf.gz*

###########################################################

# create link for perl scripts
ln -s /home/administrator/Annovar/annovar/*.pl .

Intervar.py -b hg19 \
	-i ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz --input_type=VCF \
	-o ${sample}_GATK.intervar \
	-t $INTERVARDB \
	-d $HUMANDB
rm *.pl

###########################################################

# Variant Calling Freebayes
echo "Annotate Freebayes"
freebayes-parallel $FREEBAYES_REGIONS 5 -f $REF_GENOME \
	--min-mapping-quality 0 \
	--min-alternate-fraction 0.05 \
	--min-alternate-count 2 \
	--pooled-continuous \
	--use-duplicate-reads \
	${sample}_aligned_marked.bam > ${sample}_freebayes_variants.vcf

bcftools filter -i 'INFO/DP >= 10 && INFO/AO >= 3 && (INFO/SAF >= 1 || INFO/SAR >= 1) && (INFO/AO/INFO/DP >= 0.05) && ((QUAL >= 20) || (INFO/AO >= 15 && INFO/DP >= 20))' -m '+' -s 'LOWQUAL' ${sample}_freebayes_variants.vcf > ${sample}_freebayes.filtered.vcf

bgzip ${sample}_freebayes.filtered.vcf
tabix -p vcf ${sample}_freebayes.filtered.vcf.gz

rm ${sample}_freebayes_variants.vcf

# Annotate with dbSNP
echo "annotate dbsnp"
bcftools annotate --threads 32 \
	-a $DBSNP \
	-c ID \
	-O z \
	${sample}_freebayes.filtered.vcf.gz \
	-o ${sample}_freebayes.filtered.dbsnp.vcf.gz

# Index the dbSNP annotated VCF
bcftools index ${sample}_freebayes.filtered.dbsnp.vcf.gz

rm ${sample}_freebayes.filtered.vcf.gz*

# Annotate with ClinVar
echo "Annotate ClinVar"
bcftools annotate --threads 32 -a $CLNVAR \
	-c INFO/CLNSIG,INFO/CLNDN,INFO/CLNHGVS,INFO/CLNSIGCONF -O z \
	-o ${sample}_freebayes.filtered.dbsnp.clnvar.vcf.gz ${sample}_freebayes.filtered.dbsnp.vcf.gz

# Index the ClinVar annotated VCF
bcftools index ${sample}_freebayes.filtered.dbsnp.clnvar.vcf.gz

rm ${sample}_freebayes.filtered.dbsnp.vcf.gz*

###########################################################

# Annotate with COSMIC
echo "annotate COSMIC"
bcftools annotate --threads 32 \
	-a $COSMIC \
	-c INFO/LEGACY_ID \
	-O z \
	-o ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.vcf.gz \
	${sample}_freebayes.filtered.dbsnp.clnvar.vcf.gz

bcftools index ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.vcf.gz

rm ${sample}_freebayes.filtered.dbsnp.clnvar.vcf.gz*

###########################################################

# Annotate with SnpEff
echo "annotate snpeff"
java -jar $SNPEFF ann -v hg19 -hgvs -canon -lof ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.vcf.gz | \
bcftools view --threads 32 -Ob -o ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz
bcftools index ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz

rm ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.vcf.gz*

###########################################################

# create link for perl scripts
ln -s /home/administrator/Annovar/annovar/*.pl .

Intervar.py -b hg19 \
	-i ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz --input_type=VCF \
	-o ${sample}_freebayes.intervar \
	-t $INTERVARDB \
	-d $HUMANDB
rm *.pl

###########################################################

# convert 1 -> chr1
python lcwesint2chr.py ${sample}_GATK.intervar.hg19_multianno.txt.intervar ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar
python lcwesint2chr.py ${sample}_freebayes.intervar.hg19_multianno.txt.intervar ${sample}_freebayes.intervar.hg19_multianno.txt.chr.intervar

# Keep only nessacary columns from .txt.intervar file
python lcwesExtCol1.py ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar ${sample}_GATK.intervar.txt
python lcwesExtCol1.py ${sample}_freebayes.intervar.hg19_multianno.txt.chr.intervar ${sample}_freebayes.intervar.txt

# Keep only nessecary columns from multianno file
python lcwesExtCol2.py ${sample}_GATK.intervar.hg19_multianno.txt ${sample}_GATK.multianno.txt
python lcwesExtCol2.py ${sample}_freebayes.intervar.hg19_multianno.txt ${sample}_freebayes.multianno.txt

# Merge intervar & multianno
python lcwesmerge.py ${sample}_GATK.intervar.txt ${sample}_GATK.multianno.txt ${sample}_GATK.var.txt
python lcwesmerge.py ${sample}_freebayes.intervar.txt ${sample}_freebayes.multianno.txt ${sample}_freebayes.var.txt

# Re-order
cut -f1,2,3,4,5,22,6,7,8,9,10,11,12,13,14,15,16,17,23,24,25,26,27,28,29,30,18,19,20,21,31,32,33 ${sample}_GATK.var.txt > ${sample}_GATK.variants.txt
cut -f1,2,3,4,5,22,6,7,8,9,10,11,12,13,14,15,16,17,23,24,25,26,27,28,29,30,18,19,20,21,31,32,33 ${sample}_freebayes.var.txt > ${sample}_freebayes.variants.txt

###########################################################

# Merge GATK & Freebaes variants
python lcwesmerge2.py \
	--gatk ${sample}_GATK.variants.txt \
	--freebayes ${sample}_freebayes.variants.txt \
	--output ${sample}_variants.txt

rm ${sample}_GATK.variants.txt ${sample}_freebayes.variants.txt

###########################################################

# Variant Prioritization
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Prioritize Variants..."
python prioritize_variants.py ${sample}_variants.txt ${sample}_variants.prioritized.tsv

###########################################################

# First run with GATK variants
python lcmatchmulti.py --tsv ${sample}_variants.prioritized.tsv \
	--gatk_vcf ${sample}_GATK.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz \
	--freebayes_vcf ${sample}_freebayes.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz \
	--output ${sample}_variants_semistructured.tmp

###########################################################

# Re-order columns
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $4, $5, $7, $44, $38, $46, $39, $40, $11, $41, $42, $31, $12, $43, $47, $21, $22, $18, $19, $20}' ${sample}_variants_semistructured.tmp > ${sample}_reorganized_variants.tsv

# Split intervar-evidence
python lcsplitInt.py ${sample}_reorganized_variants.tsv ${sample}_out.tsv

# clean up
rm 327427_GATK.intervar.* 327427_GATK.multianno.txt 327427_GATK.var* 327427_reorganized_variants.tsv 327427_variants_semistructured.tmp variant_merge_20250228_210208.log 327427_GATK.filtered.dbsnp.clnvar.cosmic.snpeff.vcf.gz.avinput

###########################################################

# Create Report
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Creating Report File..."
python lcrep.py ${sample}_out.tsv ${sample}.html

vcf2bed < ${sample}.vcf > ${sample}.bed

# Create IGV report
create_report ${sample}.bed --genome hg19 --flanking 1000 --tracks ${sample}_aligned_marked.bam --output ${sample}.IGV.html

###########################################################

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

echo " "
echo " "
echo "ANALYSIS COMPLETED SUCCESSFULLY!"

