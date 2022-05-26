#!/bin/bash -e

# authors : Heechul Chung, M.S.



if [ $# -lt 2 ]
then
        echo usage: $0 [prefix] [thread]
        exit 1
fi


prefix=$1
thread=$2

output_home=/home/hcchung/WES/metastasis/mappingresults/
fastqdir=/home/hcchung/WES/metastasis/GCMET_WES_fastq/matched_with_RNAseq/
ref=/data/public/GATK/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
picard=/opt/Yonsei/Picard/2.26.4/
gtf=/data/public/GATK/gcp-public-data--broad-references/hg38/v0/gencode.v27.primary_assembly.annotation.gtf

dbsnp=/data/public/GATK/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
mills=/data/public/GATK/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

output_dir=$output_home/${prefix}
mkdir $output_dir
cd ${output_dir}

echo "BWA part 1"
bwa mem $ref ${fastqdir}/${prefix}_1.fastq.gz ${fastqdir}/${prefix}_2.fastq.gz > ${prefix}_aln-pe.sam -t $thread

echo "SAM to BAM"
samtools view -bS ${prefix}_aln-pe.sam -o ${prefix}_aln-pe.bam



echo "Add read groups"
java -Xmx80G -jar ${picard}/picard.jar AddOrReplaceReadGroups \
        I=${prefix}_aln-pe.bam\
        O=${prefix}_rg_added.bam \
        RGLB=${prefix} \
        RGPL=illumina \
        RGPU=${prefix} \
        RGSM=${prefix} \
        VALIDATION_STRINGENCY=LENIENT


echo "Produce unmapped bam"
java -Xmx80G -jar ${picard}/picard.jar RevertSam \
    I=${prefix}_rg_added.bam \
    O=${prefix}_rg_added.reverse.bam \
    SANITIZE=true \


echo "Merge BAM alignment"
java -Xmx80G -jar ${picard}/picard.jar MergeBamAlignment \
      ALIGNED=${prefix}_rg_added.bam \
      UNMAPPED=${prefix}_rg_added.reverse.bam \
      O=${prefix}_merge_alignment.bam \
      R=$ref



echo "Remove duplicates"
java -Xmx80G -jar ${picard}/picard.jar MarkDuplicates \
    I=${prefix}_merge_alignment.bam \
    O=${prefix}_dedupped.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=output.metrics

echo "SortBam"
java -Xmx80G -jar ${picard}/picard.jar SortSam \
      I=${prefix}_dedupped.bam \
      O=${prefix}_sorted_dedupped.bam \
      SORT_ORDER=coordinate



echo "Generates recalibration table for Base Quality Score Recalibration"
gatk --java-options -Xmx80G BaseRecalibrator \
   -I ${prefix}_sorted_dedupped.bam \
   -R ${ref} \
   --known-sites $dbsnp \
   --known-sites $mills \
   -O recal_data.table

echo "Base quality score recalibration and get ready for variation detection"
gatk --java-options -Xmx80G ApplyBQSR \
   -R $ref \
   -I ${prefix}_sorted_dedupped.bam \
   --bqsr-recal-file recal_data.table \
   -O ${prefix}_variant_detection_ready.bam
                                                                                                                                            66,4          30%
