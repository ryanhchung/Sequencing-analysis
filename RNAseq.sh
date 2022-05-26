################Method 1##############

#!/bin/bash

# path
path='/home/jhlee/workspace/rna_seq/data'
out_path='/home/jhlee/workspace/rna_seq/hisat'
reference='/home/jhlee/workspace/test_1/0.Reference'
feature_path='/home/jhlee/workspace/rna_seq/featurecount'

# read txt and convert it as an array
mapfile -t arr < $path'/dir_name.csv'
# arr='GCMET_005_T1'

for i in ${arr[@]:1:5}
do
        A=$path/$i/$i'_1_val_1.fq'
        B=$path/$i/$i'_2_val_2.fq'
        C=$out_path/$i'.bam'
        D=$out_path/$i'.log'
    E=$feature_path/$i'_count.txt'

        F=$path/$i/$i'_1.fastq'
        G=$path/$i/$i'_2.fastq'

    # Trimming
    trim_galore --paired --fastqc $F $G -o $path/$i
        # Hisat2
        hisat2 -p 50 --rna-strandness RF -x $reference'/GRCh38' -1 $A -2 $B 2> $D| samtools view -@ 8 -Sbo $C
    # FeatureCounts
    featureCounts -T 30 -p -s 2 -t exon -g gene_id -a $reference'/Homo_sapiens.GRCh38.105.gtf' -o $E $C

        echo $i
done


################Method 2##############

#!/bin/bash -e

# authors : Heechul Chung, M.S.


if [ $# -lt 2 ]
then
        echo usage: $0 [prefix] [thread]
        exit 1
fi


prefix=$1
thread=$2

output_home=/home/data/metastasis/RNAseq/metastasis_data/mappingresults/
fastqdir=/home/data/metastasis/RNAseq/metastasis_data/GCMET_RNAseq_fastq/
ref=/home/data/reference/Homo_sapiens_assembly38.fasta
genomeDir=/home/data/metastasis/RNAseq/metastasis_data/STAR_genomedir/
gtf=/home/hcchung/reference/gencode.v39.annotation.gtf





echo "Indexing steps -> to change FASTA reference format suitable for STAR"



output_dir=$output_home/${prefix}
mkdir $output_dir
cd ${output_dir}

genomedir=${output_dir}/hg38_2pass
mkdir $genomedir


STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref --runThreadN $thread

echo "1st step done"

STAR --genomeDir $genomeDir  --readFilesIn ${fastqdir}/${prefix}_1.fastq.gz ${fastqdir}/${prefix}_2.fastq.gz --runThreadN $thread --readFilesCommand zcat

echo "2nd step done"
STAR --runMode genomeGenerate --genomeDir $genomedir --genomeFastaFiles $ref --sjdbFileChrStartEnd ${output_dir}/SJ.out.tab --sjdbOverhang 75 --runThreadN $thread

runDir=${output_dir}/STAR_output
mkdir $runDir
cd $runDir

echo "3rd step done"
STAR --genomeDir $genomedir --readFilesIn ${fastqdir}/${prefix}_1.fastq.gz ${fastqdir}/${prefix}_2.fastq.gz --runThreadN $thread --readFilesCommand zcat


echo "OK, let's convert sam to bam file"
samtools view -bS ${runDir}/Aligned.out.sam -o ${runDir}/Aligned.out.bam

echo "Let's sort and finish converting"
samtools sort ${runDir}/Aligned.out.bam -o ${runDir}/${prefix}_Aligned.out.sorted.bam


countdir=${runDir}/Count_dir
mkdir $countdir
cd $countdir

echo "HTseq start!"
htseq-count --stranded=no -m intersection-nonempty -r pos -f bam ${runDir}/${prefix}_Aligned.out.sorted.bam $gtf > ${prefix}_output.txt

#Once generated, no need to perform twice!"
#echo "Generate exon length from gtf file"
#python /home/hcchung/rnaseq/pipeline/txLength.py $gtf /home/hcchung/rnaseq/metastasis_data/mappingresults/
#Name : geneLength.txt


#cd $runDir
#echo "indexing bam files for idxstats"
#samtools index ${prefix}_Aligned.out.sorted.bam

cd $countdir

echo "Pefrofm idxstats to show chromosome length, mapped reads, and unmapped reads."
samtools idxstats ${runDir}/${prefix}_Aligned.out.sorted.bam  > idx.txt

echo "Print total reads"
python /home/ryanchung/pipeline/get_total_read.py idx.txt #After this, we get totalreadLength.txt file in countdir.

readlength=$(< totalreadlength.txt)
echo $readlength


echo "Calculate RPKM"
python /home/ryanchung/pipeline/readcount_to_rpkm.py ${prefix}_output.txt /home/data/metastasis/RNAseq/metastasis_data/mappingresults/geneLength.txt $readlength /home/data/metastasis/RNAseq/metastasis_data/mappingresults/fpkm/

mv /home/data/metastasis/RNAseq/metastasis_data/mappingresults/fpkm/rpkm.txt /home/data/metastasis/RNAseq/metastasis_data/mappingresults/fpkm/${prefix}_rpkm.txt #generate rpkm matrix with prefix
echo "Quantification Done!"

