#!/bin/bash

# Set the directory where the FASTQ files are located
FASTQ_DIR="."

# Set paths to required tools and reference files
BWA_INDEX="/home/administrator/lifecode/genomes/bwa_hg19_noalt"
REF_GENOME="/home/administrator/lifecode/genomes/bwa_hg19_noalt/hg19.noalt.fa"
REF_GENOME_dict="/home/administrator/lifecode/genomes/Picard_SDict/hg19_noalt/hg19.noalt.fa"
TARGETS="/home/administrator/lifecode/genomes/bed_files/WES_HG19/S33266340_Covered.adj.bed"
VEP_CACHEDIR="/home/administrator/.vep/104_GRCh37/"
VEP_GENOME="/home/administrator/.vep/104_GRCh37/homo_sapiens/104_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
HUMANDB="/home/administrator/lifecode/genomes/intervar_humandb/humandb"
INTERVAR="/home/administrator/lifecode/genomes/intervar_humandb/intervar"
GNOMAD_LATEST="/home/administrator/lifecode/genomes/intervar_humandb/humandb/gnomad"
AUTOPVS1="/home/administrator/D3b-autoPVS1"
CLNVAR="home/administrator/lifecode/genomes/AutoGVP_data/clinvar.vcf.gz"
CLNVAR_PUB="/home/administrator/lifecode/genomes/AutoGVP_data/ClinVar-selected-submissions.tsv"
VARIANT_SUM="/home/administrator/lifecode/genomes/AutoGVP_data/variant_summary.txt.gz"
SUB_SUM="/home/administrator/lifecode/genomes/AutoGVP_data/submission_summary.txt.gz"
CLN_IDS="/home/administrator/lifecode/genomes/AutoGVP_data/clinvar_cpg_concept_ids.txt"

# Function to process a single sample
process_sample() {
	local sample=$1
	local fastq1=$2
	local fastq2=$3

echo "Processing sample: $sample"

# Alignment
bwa mem -R "@RG\tID:${sample}\tLB:solid_lib\tPL:SOLID\tPU:unit1\tSM:${sample}" -t 32 \
	$REF_GENOME $fastq1 $fastq2 | \
	samtools view -@ 32 -bS | \
	samtools sort -@ 32 -o ${sample}_aligned_rg.bam

# Mark Dublicates
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

# variant calling
gatk HaplotypeCaller \
	-R $REF_GENOME_dict \
	-I ${sample}_aligned_marked.bam \
	-O ${sample}_variants.vcf.gz \
	-L $TARGETS \
	--native-pair-hmm-threads 5

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
	-O ${sample}_variants.filtered.vcf.gz

rm ${sample}_snps.filtered.vcf.gz* ${sample}_indels.filtered.vcf.gz*

###########################################################

# Convert vcf compatible to ensemble (chr1->1)
python lcchr2int.py ${sample}_variants.filtered.vcf.gz ${sample}_variants.filtered.chr2int.vcf.gz

gunzip ${sample}_variants.filtered.chr2int.vcf.gz

rm ${sample}_variants.filtered.vcf.gz*

###########################################################

# Annotate with vep
echo "Annotate with vep"
conda run -n VEP_104 vep \
	--offline \
	--cache \
	--dir_cache $VEP_CACHEDIR \
	--fasta $VEP_GENOME \
	--use_given_ref \
	--species homo_sapiens \
	--assembly GRCh37 \
	--fork 28 \
	--buffer_size 100000 \
	--xref_refseq \
	--hgvs \
	--hgvsg \
	--canonical \
	--symbol \
	--distance 0 \
	--exclude_predicted \
	--flag_pick \
	--lookup_ref \
	--force \
	--input_file ${sample}_variants.filtered.chr2int.vcf \
	--output_file ${sample}_vep.vcf \
	--format vcf \
	--vcf \
	--no_stats \
	--numbers

##############bgzip ${sample}_vep.vcf

###############bcftools index ${sample}_vep.vcf.gz

rm ${sample}_variants.filtered.chr2int.vcf

###########################################################

# Convert back to original (1 -> chr1)
python lcint2chr.py ${sample}_vep.vcf ${sample}_vep.int2chr.vcf

#############gunzip ${sample}_vep.int2chr.vcf.gz

rm ${sample}_vep.vcf*

# spilit mulitallelic sites

echo "Split multiallelic sites"
bcftools norm --threads 5 -m "-any" ${sample}_vep.int2chr.vcf | \
vt normalize - -n -r $REF_GENOME | \
bgzip -@ 5 -c > ${sample}_vep.int2chr.norm.vcf.gz  && \
tabix ${sample}_vep.int2chr.norm.vcf.gz

rm ${sample}_vep.int2chr.vcf

###########################################################

# Run Intervar

# create link for perl scripts
ln -s /home/administrator/Annovar/annovar/*.pl .

echo "Annotate with Intervar"
Intervar.py \
	-b hg19 \
	-i ${sample}_vep.int2chr.norm.vcf.gz \
	--input_type=VCF \
	-o ${sample}_vep.int2chr.norm.intervar.vcf \
	-t $INTERVAR \
	-d $HUMANDB

###########################################################

# Run Annovar

echo "Annotate Annovar"
perl table_annovar.pl ${sample}_vep.int2chr.norm.vcf.gz $GNOMAD_LATEST \
	--buildver hg19 \
	--outfile ${sample}_vep.int2chr.norm.anno \
	--protocol gnomad211_exome,gnomad211_genome \
	--operation f,f \
	--remove \
	--vcfinput

rm *.pl

###########################################################

# Run AutoPVS1

echo "Run AutoPVS1"

# Switch to AUTOPVS1 dir
my_dir=$(pwd)
cd $AUTOPVS1

# Run autoPVS1
python autoPVS1_from_VEP_vcf.py --genome_version genome_hg19 --vep_vcf $my_dir/${sample}_vep.int2chr.norm.vcf.gz > ${sample}_PVS1.txt
mv ${sample}_PVS1.txt $my_dir/
cd $my_dir

###########################################################

# Run AutoGVP
gunzip ${sample}_vep.int2chr.norm.vcf.gz
ln -s /home/administrator/lifecode/genomes/AutoGVP_data/* .
mkdir results

echo "Running AutoGVP"
conda run -n bcftools_1.17 run_autogvp.sh \
	--workflow="custom" \
	--vcf=${sample}_vep.int2chr.norm.vcf \
	--clinvar=clinvar.vcf.gz \
	--intervar=${sample}_vep.int2chr.norm.intervar.vcf.hg19_multianno.txt.intervar \
	--multianno=${sample}_vep.int2chr.norm.anno.hg19_multianno.txt \
	--autopvs1=${sample}_PVS1.txt \
	--outdir=results \
	--out="${sample}" \
	--selected_clinvar_submissions=ClinVar-selected-submissions.tsv \
	--variant_summary=variant_summary.txt.gz \
	--submission_summary=submission_summary.txt.gz \
	--conceptIDs=clinvar_cpg_concept_ids.txt \
	--conflict_res="latest"

###########################################################

echo "Finished processing sample: $sample"

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

###########################################################
###########################################################
###########################################################

# Initiate bio_prior script
#python bio_WES_prior_v2.py ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.vcf ${sample}.xlsx

# Create html report
#python bio_WES_report_v0.1.py ${sample}.xlsx ${sample}.html

# convert vcf to bed for creating IGV report
#vcf2bed < ${sample}.vcf > ${sample}.bed

# Create IGV report
#create_report ${sample}.bed --genome hg19 --flanking 1000 --tracks ${sample}_aligned_marked.bam --output ${sample}.IGV.html

# clean up
#rm ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.intervar.hg19_multianno.txt*
#rm ${sample}_variants.filtered.clnvar.dbsnp.snpeff.gnomAD.cosmic.vcf*

echo " "
echo " "
echo "PIPELINE COMPLETE! Results saved to ${sample}.xlsx"

