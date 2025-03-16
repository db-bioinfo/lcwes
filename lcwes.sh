#!/bin/bash

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files
THREADS=32
REF_GENOME="/home/administrator/lifecode/genomes/databases/bwa_hg19/hg19.fa"
TARGETS="/home/administrator/lifecode/genomes/bed_files/WES_HG19/S33266340_Covered.adj.bed"
SNPEFF="/home/administrator/SNPEFF/snpEff/snpEff.jar"
INTERVARDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/intervar"
HUMANDB="/home/administrator/lifecode/genomes/databases/intervar_humandb/humandb"
FREEBAYES_REGIONS="/home/administrator/lifecode/genomes/freebayes_hg19/hg19_regions.txt"

#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------

# Function to process a single sample
process_sample() {
	local sample=$1
	local fastq1=$2
	local fastq2=$3

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

#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------

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
	python lcwesmer.py ${sample}_GATK.variants.prioritized.tsv ${sample}_freebayes.variants.prioritized.20.tsv > ${sample}_variants.tsv

	# Create Report
	python lcwesrep.py ${sample}_variants.tsv ${sample}.html

	# Create bed file
	awk 'NR > 1 {print $1 "\t" $2 "\t" $3 "\t" $6}' ${sample}_variants.tsv | head -n 300 > ${sample}_variants.bed

	# Create IGV report
	create_report ${sample}_variants.bed --genome hg19 --flanking 1000 --tracks ${sample}_aligned_marked.bam --output ${sample}.IGV.html

	rm *.pl

else
	echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR"
	[ $GATK_STATUS -ne 0 ] && echo "process 1 failed with status $GATK_STATUS"
	[ $FREEBAYES_STATUS -ne 0 ] && echo "process 2 failed with status $FREEBAYES_STATUS"
fi

}

#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------

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

#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------

# Convert vcf to Annovar format
convert2annovar.pl -format vcf4 ${sample}_GATK.filtered.vcf.gz > ${sample}_GATK.avinput
# Annovar hgvs annotation
annotate_variation.pl -out ${sample}_GATK_hgvs -build hg19 -hgvs ${sample}_GATK.avinput $HUMANDB
# Re order
./lcwesreor.sh ${sample}_GATK_hgvs.exonic_variant_function > ${sample}_GATK_hgvs.exonic_variant_function.reorder
# cleanup
rm ${sample}_GATK_hgvs.log ${sample}_GATK_hgvs.variant_function ${sample}_GATK_hgvs.exonic_variant_function

# Intervar/Annovar annotation
Intervar.py -b hg19 \
	-i ${sample}_GATK.filtered.vcf.gz --input_type=VCF \
	-o ${sample}_GATK.intervar \
	-t $INTERVARDB \
	-d $HUMANDB

# convert 1 -> chr1
python lcwesint2chr.py ${sample}_GATK.intervar.hg19_multianno.txt.intervar ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar

# Keep only nessacary columns from .txt.intervar file
python lcwesExtCol1.py ${sample}_GATK.intervar.hg19_multianno.txt.chr.intervar ${sample}_GATK.intervar.txt

# Keep only nessecary columns from multianno file
python lcwesExtCol2.py ${sample}_GATK.intervar.hg19_multianno.txt ${sample}_GATK.multianno.txt

# Merge intervar & multianno
python lcwesmerge.py ${sample}_GATK.intervar.txt ${sample}_GATK.multianno.txt ${sample}_GATK.var.txt

# Re-order
cut -f1,2,3,4,5,22,6,7,8,9,10,11,12,13,14,15,16,17,23,24,25,26,27,28,29,30,18,19,20,21,31,32,33 ${sample}_GATK.var.txt > ${sample}_GATK.variants.txt

# Variant Prioritization
python prioritize_variants.py ${sample}_GATK.variants.txt ${sample}_GATK.variants.prioritized.tmp

# Split intervar to ACMG
python lcwessplit.py ${sample}_GATK.variants.prioritized.tmp ${sample}_GATK.variants.prioritized.2.tmp

# Add hgvs to final output
python lcweshgvs.py ${sample}_GATK.variants.prioritized.2.tmp ${sample}_GATK_hgvs.exonic_variant_function.reorder ${sample}_GATK.variants.prioritized.tsv

# cleanup
rm ${sample}_GATK.intervar.* ${sample}_GATK.multianno.txt ${sample}_GATK.var.txt ${sample}_GATK.variants.txt
rm ${sample}_GATK.filtered.* variant_merge* ${sample}_GATK_hgvs.exonic_variant_function.reorder
rm ${sample}_GATK.variants.prioritized.tmp ${sample}_GATK.variants.prioritized.2.tmp
rm ${sample}_GATK.avinput ${sample}_GATK.variants.prioritized_stats.json

}

#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------

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

#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------
#==========#----------#==========#----------#==========#----------#==========#----------

# Convert vcf to Annovar format
convert2annovar.pl -format vcf4 ${sample}_freebayes.filtered.vcf.gz > ${sample}_freebayes.avinput
# Annovar hgvs annotation
annotate_variation.pl -out ${sample}_freebayes_hgvs -build hg19 -hgvs ${sample}_freebayes.avinput $HUMANDB
# Re order
./lcwesreor.sh ${sample}_freebayes_hgvs.exonic_variant_function > ${sample}_freebayes_hgvs.exonic_variant_function.reorder
# cleanup
rm ${sample}_freebayes_hgvs.log ${sample}_freebayes_hgvs.variant_function ${sample}_freebayes_hgvs.exonic_variant_function

# Intervar/Annovar Annotation
Intervar.py -b hg19 \
	-i ${sample}_freebayes.filtered.vcf.gz --input_type=VCF \
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

# Variant Prioritization
python prioritize_variants.py ${sample}_freebayes.variants.txt ${sample}_freebayes.variants.prioritized.tmp

# Split intervar to ACMG
python lcwessplit.py ${sample}_freebayes.variants.prioritized.tmp ${sample}_freebayes.variants.prioritized.2.tmp

# Add hgvs to final output
python lcweshgvs.py ${sample}_freebayes.variants.prioritized.2.tmp ${sample}_freebayes_hgvs.exonic_variant_function.reorder ${sample}_freebayes.variants.prioritized.tsv

head -n 20 ${sample}_freebayes.variants.prioritized.tsv > ${sample}_freebayes.variants.prioritized.20.tsv

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
