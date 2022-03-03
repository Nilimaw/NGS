#!/bin/bash

#SBATCH --array=1-95
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH -J dnaseq
#SBATCH -e dnaseq.err%j
#SBATCH -o dnaseq.out%j
#SBATCH -t 24:00:00
#SBATCH -n 4
#SBATCH -p standard
#SBATCH --mail-type=ALL

#Load required modules
module load bowtie2
module load fastqc/0.11.5
module load java
module load perl
module load samtools
ID=${SLURM_ARRAY_TASK_ID}

echo "ID1= $ID"
declare -a arr1=("N711-N508" "N708-N508" "N710-N507" "N701-N508" "N702-N508" "N707-N507" "N705-N508" "N707-N508" "N709-N508" "N711-N507" "N703-N508" "N709-N502" "N712-N506" "N708-N506" "N709-N506" "N710-N506" "N711-N506" "N710-N502" "N701-N507" "N702-N507" "N703-N507" "N711-N502" "N704-N507" "N705-N507" "N709-N507" "N701-AK394" "N702-AK394" "N703-AK394" "N704-AK394" "N705-AK394" "N706-AK394" "N707-AK394" "N708-AK394" "N709-AK394" "N710-AK394" "N711-AK394" "N712-AK394" "N701-N503" "N702-N503" "N703-N503" "N704-N503" "N705-N503" "N701-N502" "N710-N503" "N706-N503" "N707-N503" "N708-N503" "N709-N503" "N702-N502" "N711-N503" "N712-N503" "N701-N504" "N702-N504" "N703-N504" "N703-N502" "N704-N504" "N704-N505" "N706-N502" "N709-N505" "N705-N505" "N706-N505" "N706-N504" "N704-N502" "N711-N504" "N707-N504" "N708-N504" "N709-N504" "N710-N504" "N705-N502" "N712-N504" "N701-N505" "N702-N505" "N703-N505" "N708-N505" "N707-N502" "N702-N506" "N710-N505" "N711-N505" "N701-N506" "N708-N502" "N707-N506" "N703-N506" "N704-N506" "N706-N506" "N707-N505" "N705-N504" "N710-N508" "N712-N505" "N712-N507" "N705-N506" "N706-N507" "N708-N507" "N712-N502" "N704-N508" "N706-N508")


ID2=${arr1[$ID-1]}
echo "ID2= $ID2"


#Set working directory and data directory
workingdir="/scratch/jfay6_lab/Nilima/aneuploids"
mkdir $workingdir/sam_files_global_nuc $workingdir/bam_files_global_nuc
data="/scratch/jfay6_lab/Nilima/novogene_1_31_20/fastq"
mkdir /scratch/jfay6_lab/Nilima/novogene_1_31_20/fastqc
echo "Working in dir: $workingdir"

#Run FastQC
final="/scratch/jfay6_lab/Nilima/novogene_1_31_20/fastqc"

#for file in $data/*.fq.gz
#do
fastqc -o $final $data/*-${ID2}_H773KCCX2_L1_1.fq.gz
fastqc -o $final $data/*-${ID2}_H773KCCX2_L1_2.fq.gz

#doneq

echo "Fastqc done"

#Set genome fasta and gff for Bowtie
REFERENCEINDEX="/scratch/jfay6_lab/Nilima/reference/hybridS.fasta"
REFERENCEGFF="/scratch/jfay6_lab/Nilima/reference/hybridS.gff"

#Build bowtie index
bowtie2-build $REFERENCEINDEX $REFERENCEINDEX

#Align and combine reads from 4 lanes

bowtie2 --local -x $REFERENCEINDEX -1 $data/*-${ID2}_H773KCCX2_L1_1.fq.gz -2 $data/*-${ID2}_H773KCCX2_L1_2.fq.gz -N 1 -q -p 4 -S $workingdir/sam_files_global_nuc/${ID2}.sam

samtools view -bS $workingdir/sam_files_global_nuc/${ID2}.sam > $workingdir/bam_files_global_nuc/${ID2}.bam

java -jar /software/picard/2.12.0/picard.jar AddOrReplaceReadGroups INPUT=$workingdir/bam_files_global_nuc/${ID2}.bam OUTPUT=$workingdir/bam_files_global_nuc/${ID2}.sorted.bam SORT_ORDER=coordinate RGID=${ID2} RGLB=${ID2} RGPL=illumina RGPU=${ID2} RGSM=${ID2}

java -Xmx8g -jar /software/picard/2.12.0/picard.jar MarkDuplicates INPUT=$workingdir/bam_files_global_nuc/${ID2}.sorted.bam OUTPUT=$workingdir/bam_files_global_nuc/${ID2}.dedup.bam METRICS_FILE=$workingdir/bam_files_global_nuc/${ID2}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true

#samtools merge $workingdir/bam_files_global_nuc/N703-N508_L1_all.dedup.bam $workingdir/bam_files_global_nuc/N703-N508_L1.dedup.bam

samtools view $workingdir/bam_files_global_nuc/${ID2}.dedup.bam > $workingdir/sam_files_global_nuc/${ID2}_all.dedup.sam
samtools view -H $workingdir/bam_files_global_nuc/${ID2}.dedup.bam >$workingdir/sam_files_global_nuc/${ID2}_all.header.sam
perl /scratch/jfay6_lab/Nilima/code/AS_check_050118.pl $workingdir/sam_files_global_nuc/${ID2}_all.dedup.sam $workingdir/sam_files_global_nuc/${ID2}_all.as.sam

echo "Alignment and merge done"

#Annotate with chr pos
cat $workingdir/sam_files_global_nuc/${ID2}_all.header.sam $workingdir/sam_files_global_nuc/${ID2}_all.as.sam >$workingdir/sam_files_global_nuc/${ID2}_all.reheader.sam
samtools view -bS $workingdir/sam_files_global_nuc/${ID2}_all.reheader.sam > $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam

#Index
samtools index $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam

echo "Indexing done"

#Chromosome counts - mapped and unmapped reads to each chromosome
samtools idxstats $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam > $workingdir/bam_files_global_nuc/${ID2}_all.chr.counts

#Depth of coverage at each pos of each chromosome
mkdir $workingdir/bam_files_global_nuc/coverage $workingdir/bam_files_global_nuc/coverage_analysis
samtools depth $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam > $workingdir/bam_files_global_nuc/coverage/${ID2}.coverage

#Classify covergae into chromosome and species
declare -a arr=("chrI" "chrII" "chrIII" "chrIV" "chrV" "chrVI" "chrVII" "chrVIII" "chrIX" "chrX" "chrXI" "chrXII" "chrXIII" "chrXIV" "chrXV" "chrXVI" "Sbay_1" "Sbay_2" "Sbay_3" "Sbay_4" "Sbay_5" "Sbay_6" "Sbay_7" "Sbay_8" "Sbay_9" "Sbay_10" "Sbay_11" "Sbay_12" "Sbay_13" "Sbay_14" "Sbay_15" "Sbay_16" "chrmt" "r23N702_mtDNA_gapped")

for file in "${arr[@]}"
do
echo "$file"
mkdir $workingdir/bam_files_global_nuc/coverage_analysis/$file
awk '$1 == "'$file'"' $workingdir/bam_files_global_nuc/coverage/${ID2}.coverage > $workingdir/bam_files_global_nuc/coverage_analysis/$file/${file}_${ID2}.coverage
done

#Mapping stats
mkdir $workingdir/bam_files_global_nuc/mapstats
samtools flagstat $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam > $workingdir/bam_files_global_nuc/mapstats/${ID2}.mapstats

#Summarize mapping stats in a file
#total reads
o1=$(samtools view -c $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam)

#mapped reads reads 
o2=$(samtools view -F 0x04 -c $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam)

#unmapped reads
o3=$(samtools view -f4 -c $workingdir/bam_files_global_nuc/${ID2}_all.reheader.bam)

echo ${ID2} $o1 $o2 $o3 >> $workingdir/bam_files_global_nuc/mapstats.txt

echo "Done!"



