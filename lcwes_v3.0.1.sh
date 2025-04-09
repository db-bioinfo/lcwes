#!/bin/bash

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files
# Alignment-Variant Calling
REF_GENOME="/home/administrator/lifecode/genomes/databases/bwa_hg19/hg19.fa"
TARGETS="/home/administrator/lifecode/genomes/bed_files/WES_HG19/S33266340_Covered.adj.bed"
# Annotation
SNPEFF_JAR="/home/administrator/snpeff/snpEff/snpEff.jar"
CLINVAR_VCF="/home/administrator/lifecode/genomes/databases/clnvar_hg19/clinvar.chr.vcf.gz"
DBSNP_VCF="/home/administrator/lifecode/genomes/databases/dbsnp_hg19/dbsnp_hg19_00-All.vcf.gz"
# Downstream analysis
INTERVARDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/hg19/intervar"
HUMANDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/hg19/humandb"
FREEBAYES_REGIONS="/home/administrator/lifecode/genomes/databases/freebayes_regions_hg19/hg19_regions.txt"
# Computation
THREADS=32

# Function to process a single sample
process_sample() {
	local sample=$1
	local fastq1=$2
	local fastq2=$3

#-------------------------- Filtering & Alignment ---------------------------#

# Quality Control / Trimming
conda run -n FASTP fastp --in1 $fastq1 --in2 $fastq2 \
	--out1 ${sample}_t1.fq.gz --out2 ${sample}_t2.fq.gz \
	--detect_adapter_for_pe \
	--html report.html \
	--json report.json \
	--thread 16

# Alignment
bwa mem -R "@RG\tID:${sample}\tLB:exome_lib\tPL:MGISEQ\tPU:unit1\tSM:${sample}" -t $THREADS \
	$REF_GENOME $fastq1 $fastq2 | \
	samtools view -@ $THREADS -bS | \
	samtools sort -@ $THREADS -o ${sample}_aligned_rg.bam

# Mark Duplicates
gatk MarkDuplicates \
	-I ${sample}_aligned_rg.bam \
	-O ${sample}_aligned_marked.bam \
	-M ${sample}_output.metrics.txt \
	--ASSUME_SORT_ORDER coordinate \
	--CREATE_INDEX true \
	--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

samtools index ${sample}_aligned_marked.bam

rm ${sample}_aligned_rg.bam*

# Prepare annovar scripts
ln -s /home/administrator/Annovar/annovar/*.pl .

#-------------------------- Merge Variants & Report ---------------------------#

# Launch GATK in the background
process_gatk "$sample" &
GATK_PID=$!

# Launch Freebayes in the background
process_freebayes "$sample" &
FREEBAYES_PID=$!

# Wait for both processes to complete
wait $GATK_PID
GATK_STATUS=$?

wait $FREEBAYES_PID
FREEBAYES_STATUS=$?

# Check if both processes completed successfully
if [ $GATK_STATUS -eq 0 ] && [ $FREEBAYES_STATUS -eq 0 ]; then

	echo "Merging"
	# Merge gatk & top 20 freebayes
	python lcwesmer.py ${sample}_gatk.variants.annotated.prioritized.tsv ${sample}_freebayes.variants.annotated.prioritized.20.tsv > ${sample}_variants.tsv

	# Create Report
	python lcwesrep.py ${sample}_variants.tsv ${sample}.html

	# Create bed file
	awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $6}' ${sample}_variants.tsv | head -n 500 > ${sample}_variants.bed

	# Create IGV report
	create_report ${sample}_variants.bed --genome hg19 --flanking 1000 --tracks ${sample}_aligned_marked.bam --output ${sample}.IGV.html

	rm *.pl

else
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR"
	[ $GATK_STATUS -ne 0 ] && echo "process 1 failed with status $GATK_STATUS"
	[ $FREEBAYES_STATUS -ne 0 ] && echo "process 2 failed with status $FREEBAYES_STATUS"
fi

}

#-------------------------- GATK Variant Calling ---------------------------#

# GATK variant calling
process_gatk() {
	local sample=$1

# Variant calling GATK
gatk HaplotypeCaller \
	-R $REF_GENOME \
	-I ${sample}_aligned_marked.bam \
	-O ${sample}_variants.vcf.gz \
	--native-pair-hmm-threads $THREADS

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

#-------------------------- GATK Variant Annoation ---------------------------#

# Annotate with SnpEff
java -jar $SNPEFF_JAR ann -v hg19 ${sample}_GATK.filtered.vcf.gz | \
bcftools view --threads $THREADS -Oz -o ${sample}_gatk.vcf

rm ${sample}_GATK_filtered.vcf.gz*

# Intervar/Annovar annotation
Intervar.py -b hg19 \
	-i ${sample}_gatk.vcf --input_type=VCF \
	-o ${sample}_gatk.intervar \
	-t $INTERVARDB \
	-d $HUMANDB

# convert 1 -> chr1
python lcwesint2chr.py ${sample}_gatk.intervar.hg19_multianno.txt.intervar ${sample}_gatk.intervar.hg19_multianno.txt.chr.intervar

# Keep only nessacary columns from .txt.intervar file
python lcwesExtCol1.py ${sample}_gatk.intervar.hg19_multianno.txt.chr.intervar ${sample}_gatk.intervar.txt

# Keep only nessecary columns from multianno file
python lcwesExtCol2.py ${sample}_gatk.intervar.hg19_multianno.txt ${sample}_gatk.multianno.txt

# Merge intervar & multianno
python lcwesmerge.py ${sample}_gatk.intervar.txt ${sample}_gatk.multianno.txt ${sample}_gatk.var.txt

# Re-order
cut -f1,2,3,4,5,22,6,7,8,9,10,11,12,13,14,15,16,17,23,24,25,26,27,28,29,30,18,19,20,21,31,32,33 ${sample}_gatk.var.txt > ${sample}_gatk.variants.txt

#-------------------------- ADD Hgvs-c Hgvs-p, Effect, CLNHGVS VAF/AF, PASS/Quality ---------------------------#

# Hgvs-c/Hgvs-p/Effect
# Extract only the highest impact annotation for each variant
conda run -n SNPSIFT SnpSift extractFields ${sample}_gatk.vcf \
  CHROM POS REF ALT \
  "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  >> ${sample}_gatk.snpsift.tsv

# bgzip
bgzip ${sample}_gatk.vcf
tabix -p vcf ${sample}_gatk.vcf.gz

# CLNHGVS id
bcftools annotate --threads $THREADS -a $CLINVAR_VCF \
	-c CLNHGVS,INFO -O z \
	-o ${sample}_gatk.clnvar.vcf.gz ${sample}_gatk.vcf.gz

# Index the ClinVar annotated VCF
tabix -p vcf ${sample}_gatk.clnvar.vcf.gz

# RS id
bcftools annotate --threads $THREADS -a $DBSNP_VCF \
	-c ID \
	-o ${sample}_gatk.clnvar.dbsnp.vcf.gz ${sample}_gatk.clnvar.vcf.gz

# Index file
tabix -p vcf ${sample}_gatk.clnvar.dbsnp.vcf.gz

# Extract qual,af,clnhgvs
python lcwesClnExt.py ${sample}_gatk.clnvar.dbsnp.vcf.gz ${sample}_gatk.clnvar.tsv

# Merge clnvar snpsift
python lcwesmergeSnpSift2Clnvar.py ${sample}_gatk.clnvar.tsv ${sample}_gatk.snpsift.tsv ${sample}_gatk.clnsnp.tsv

# Merge clnsnp to gatk
python lcwesmergeClnSnp2Prior.py --threads 32 ${sample}_gatk.variants.txt ${sample}_gatk.clnsnp.tsv ${sample}_gatk.variants.annotated.tsv

#-------------------------- GATK Variant Prioritization ---------------------------#

# Variant Prioritization
python prioritize_variants.py ${sample}_gatk.variants.annotated.tsv ${sample}_gatk.variants.prioritized.tmp

# Split intervar to ACMG
python lcwessplit.py ${sample}_gatk.variants.prioritized.tmp ${sample}_gatk.variants.annotated.prioritized.tsv

# cleanup
rm ${sample}_gatk.intervar.* ${sample}_gatk.multianno.txt ${sample}_gatk.var.txt ${sample}_gatk.variants.txt
rm ${sample}_gatk.filtered.* variant_merge* ${sample}_gatk_hgvs.exonic_variant_function.reorder
rm ${sample}_gatk.variants.prioritized.tmp ${sample}_gatk.variants.prioritized.2.tmp
rm ${sample}_gatk.avinput ${sample}_gatk.variants.prioritized_stats.json

}

#-------------------------- FREEBAYES VARIANT CALLING ---------------------------#

process_freebayes() {
	local sample=$1

# Variant Calling Freebayes
freebayes-parallel $FREEBAYES_REGIONS $THREADS -f $REF_GENOME \
	--min-mapping-quality 0 \
	--min-alternate-fraction 0.05 \
	--min-alternate-count 2 \
	--pooled-continuous \
	--use-duplicate-reads \
	${sample}_aligned_marked.bam > ${sample}_freebayes_variants.vcf

bcftools filter -i 'INFO/DP >= 10 && INFO/AO >= 3 && ((INFO/SAF >= 1 && INFO/SAR >= 1) || INFO/AO >= 9) && (INFO/AO/INFO/DP >= 0.05 || INFO/AO >= 30) && (QUAL >= 20 || (INFO/AO >= 9 && INFO/DP >= 35))' -m '+' -s 'LOWQUAL' ${sample}_freebayes_variants.vcf > ${sample}_freebayes.filtered.tmp.vcf

# Keep only pass variants
bcftools view -f PASS ${sample}_freebayes.filtered.tmp.vcf -Oz -o ${sample}_freebayes.filtered.vcf

rm ${sample}_freebayes.filtered.tmp.vcf

bgzip ${sample}_freebayes.filtered.vcf
tabix -p vcf ${sample}_freebayes.filtered.vcf.gz

#-------------------------- FREEBAYES VARIANT ANNOTATION ---------------------------#

# Annotate with SnpEff
java -jar $SNPEFF_JAR ann -v hg19 ${sample}_freebayes.filtered.vcf.gz | \
bcftools view --threads $THREADS -Oz -o ${sample}_freebayes.vcf

rm ${sample}_freebayes_filtered.vcf.gz*

# Intervar/Annovar annotation
Intervar.py -b hg19 \
	-i ${sample}_freebayes.vcf --input_type=VCF \
	-o ${sample}_freebayes.intervar \
	-t $INTERVARDB \
	-d $HUMANDB

# convert 1 -> chr1
python lcwesint2chr.py ${sample}_freebayes.intervar.hg19_multianno.txt.intervar ${sample}_freebayes.intervar.hg19_multianno.txt.chr.intervar

# Keep only nessacary columns from .txt.intervar file
python lcwesExtCol1.py ${sample}_freebayes.intervar.hg19_multianno.txt.chr.intervar ${sample}_freebayes.intervar.txt

# Keep only nessecary columns from multianno file
python lcwesExtCol2.py ${sample}_freebayes.intervar.hg19_multianno.txt ${sample}_freebayes.multianno.txt

# Merge intervar & multianno
python lcwesmerge.py ${sample}_freebayes.intervar.txt ${sample}_freebayes.multianno.txt ${sample}_freebayes.var.txt

# Re-order
cut -f1,2,3,4,5,22,6,7,8,9,10,11,12,13,14,15,16,17,23,24,25,26,27,28,29,30,18,19,20,21,31,32,33 ${sample}_freebayes.var.txt > ${sample}_freebayes.variants.txt

#-------------------------- ADD Hgvs-c Hgvs-p, Effect, CLNHGVS VAF/AF, PASS/Quality ---------------------------#

# Hgvs-c/Hgvs-p/Effect
# Extract only the highest impact annotation for each variant
conda run -n SNPSIFT SnpSift extractFields ${sample}_freebayes.vcf \
  CHROM POS REF ALT \
  "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  >> ${sample}_freebayes.snpsift.tsv

# bgzip
bgzip ${sample}_freebayes.vcf
tabix -p vcf ${sample}_freebayes.vcf.gz

# CLNHGVS id
bcftools annotate --threads $THREADS -a $CLINVAR_VCF \
	-c CLNHGVS,INFO -O z \
	-o ${sample}_freebayes.clnvar.vcf.gz ${sample}_freebayes.vcf.gz

# Index the ClinVar annotated VCF
tabix -p vcf ${sample}_freebayes.clnvar.vcf.gz

# RS id
bcftools annotate --threads $THREADS -a $DBSNP_VCF \
	-c ID \
	-o ${sample}_freebayes.clnvar.dbsnp.vcf.gz ${sample}_freebayes.clnvar.vcf.gz

# Index file
tabix -p vcf ${sample}_freebayes.clnvar.dbsnp.vcf.gz

# Extract qual,af,clnhgvs
python lcwesClnExt.py ${sample}_freebayes.clnvar.dbsnp.vcf.gz ${sample}_freebayes.clnvar.tsv

# Merge clnvar snpsift
python lcwesmergeSnpSift2Clnvar.py ${sample}_freebayes.clnvar.tsv ${sample}_freebayes.snpsift.tsv ${sample}_freebayes.clnsnp.tsv

# Merge clnsnp to freebayes
python lcwesmergeClnSnp2Prior.py --threads 32 ${sample}_freebayes.variants.txt ${sample}_freebayes.clnsnp.tsv ${sample}_freebayes.variants.annotated.tsv

#-------------------------- Variant Prioritization ---------------------------#

# Variant Prioritization
python prioritize_variants.py ${sample}_freebayes.variants.annotated.tsv ${sample}_freebayes.variants.prioritized.tmp

# Split intervar to ACMG
python lcwessplit.py ${sample}_freebayes.variants.prioritized.tmp ${sample}_freebayes.variants.annotated.prioritized.tsv

# Keep first 20 prioritized variants
head -n 20 ${sample}_freebayes.variants.annotated.prioritized.tsv > ${sample}_freebayes.variants.annotated.prioritized.20.tsv

# cleanup
rm ${sample}_freebayes.intervar.* ${sample}_freebayes.multianno.txt ${sample}_freebayes.var.txt ${sample}_freebayes.variants.txt
rm ${sample}_freebayes.filtered.* variant_merge* ${sample}_freebayes_hgvs.exonic_variant_function.reorder
rm ${sample}_freebayes.variants.prioritized.tmp ${sample}_freebayes.variants.prioritized.2.tmp
rm ${sample}_freebayes.avinput ${sample}_freebayes.variants.prioritized_stats.json

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

echo "ALL DONE!"
