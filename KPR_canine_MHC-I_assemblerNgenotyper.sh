#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=KPR
#SBATCH --ntasks=4
#SBATCH --time=72:00:00
#SBATCH --mem=60G

##################################################################
####################### Sample Info ##############################
##################################################################

script='' # path to the KPR scripts
data='' # path to the RNAseq paired end fastq files
result='' # path to the output directory
ref='' # path to the KPR reference


fastq1=''
fastq2=''

K=50
N=1000

DLA88_cutoffDepth=15
DLA12_cutoffDepth=38
cutoffLength=270

module load SAMtools/1.9-GCC-8.3.0
module load Trimmomatic/0.39-Java-1.8.0_144
module load BWA/0.7.17-GCC-8.3.0
module load SRA-Toolkit/2.10.8-centos_linux64
module load Clustal-Omega/1.2.4-GCC-8.3.0


##################################################################
######################## Execution ###############################
##################################################################

script_assembler=$script'/KPR_assembler'
script_genotyper=$script'/KPR_genotyper'
result_tmp=$result'/execution'

mkdir $result
mkdir $result_tmp
cd $result_tmp

# Trimmomatics trimming
time java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 \
$data/$fastq1 $data/$fastq2 \
sample_lane1_forward_paired.fq.gz sample_lane1_forward_unpaired.fq.gz \
sample_lane1_reverse_paired.fq.gz sample_lane1_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# DLA-I read extraction
time bwa aln $ref/DLA_Nucleotide_DB_GenBank_AllExons.fa sample_lane1_forward_paired.fq.gz > sample_lane1_forward_paired.sai
time bwa aln $ref/DLA_Nucleotide_DB_GenBank_AllExons.fa sample_lane1_reverse_paired.fq.gz > sample_lane1_reverse_paired.sai

time bwa sampe $ref/DLA_Nucleotide_DB_GenBank_AllExons.fa \
sample_lane1_forward_paired.sai \
sample_lane1_reverse_paired.sai \
sample_lane1_forward_paired.fq.gz \
sample_lane1_reverse_paired.fq.gz \
> sample_DLAdb.sam

samtools view -bS sample_DLAdb.sam > sample_DLAdb.bam
samtools sort sample_DLAdb.bam -o sample_DLAdb-sorted.bam
samtools index sample_DLAdb-sorted.bam

samtools view sample_DLAdb-sorted.bam | grep -v '^@' | awk '$2%16 == 3 && $3 != "*" && $6 !~ /S/ && $6 !~ /H/ {print $0}' > sample_BWA_reads.sam

gawk '(and(16, $2))' sample_BWA_reads.sam > DLA_reverseStrandReads.sam
gawk '(! and(16,$2))' sample_BWA_reads.sam > DLA_forwardStrandReads.sam

cat DLA_reverseStrandReads.sam | awk '$4 <= 445 {print $0}' > DLA_Exon2_reverseStrandReads.sam
cat DLA_forwardStrandReads.sam | awk '$4 <= 445 {print $0}' > DLA_Exon2_forwardStrandReads.sam

cat DLA_reverseStrandReads.sam | awk '$4 <= 725 {print $0}' > DLA_Exon2-3_reverseStrandReads.sam
cat DLA_forwardStrandReads.sam | awk '$4 <= 725 {print $0}' > DLA_Exon2-3_forwardStrandReads.sam

cat DLA_reverseStrandReads.sam | awk '$4 >= 246 && $4 <= 725 {print $0}' > DLA_Exon3_reverseStrandReads.sam
cat DLA_forwardStrandReads.sam | awk '$4 >= 246 && $4 <= 725 {print $0}' > DLA_Exon3_forwardStrandReads.sam

python $script_assembler/Sam_file_Quality_control.py DLA_Exon2_reverseStrandReads.sam 30
python $script_assembler/Sam_file_Quality_control.py DLA_Exon2_forwardStrandReads.sam 30
python $script_assembler/Sam_file_Quality_control.py DLA_Exon2-3_reverseStrandReads.sam 30
python $script_assembler/Sam_file_Quality_control.py DLA_Exon2-3_forwardStrandReads.sam 30
python $script_assembler/Sam_file_Quality_control.py DLA_Exon3_reverseStrandReads.sam 30
python $script_assembler/Sam_file_Quality_control.py DLA_Exon3_forwardStrandReads.sam 30

for i in DLA_Exon*.sam_QualityFiltered; do awk '{OFS="\t"; print ">"$1"\n"$10}' $i > "${i/%.sam_QualityFiltered/.fa}"; done

python $script_assembler/Get_Pairwise_reads.py DLA_Exon2_forwardStrandReads.fa DLA_Exon2_reverseStrandReads.fa
python $script_assembler/Get_Pairwise_reads.py DLA_Exon2-3_forwardStrandReads.fa DLA_Exon2-3_reverseStrandReads.fa
python $script_assembler/Get_Pairwise_reads.py DLA_Exon3_forwardStrandReads.fa DLA_Exon3_reverseStrandReads.fa

######################## Assembly #######################################

mkdir $result_tmp/Test4
cd $result_tmp

## 100bp
# DLA-88 1-73;74-345;346-625
# DLA-12 1-64;65-336;337-616
# DLA-64 1-73;74-345;346-625

### Forwd ###
python $script_assembler/Kmer_assembly4_V6.py DLA_Exon2_forwardStrandReads.fa-pair DLA_Exon2_reverseStrandReads.fa-pair $K $N 3 > Test4/'Kmer_Exon2_contigs-'$K'mer_'$N'.txt'
python $script_assembler/Kmer_assembly4_V6.py DLA_Exon2-3_forwardStrandReads.fa-pair DLA_Exon2-3_reverseStrandReads.fa-pair $K $N 3 > Test4/'Kmer_Exon2-3_contigs-'$K'mer_'$N'.txt'
python $script_assembler/Kmer_assembly4_V6.py DLA_Exon3_forwardStrandReads.fa-pair DLA_Exon3_reverseStrandReads.fa-pair $K $N 3 > Test4/'Kmer_Exon3_contigs-'$K'mer_'$N'.txt'

cd $result_tmp/Test4
mv 'Kmer_Exon2_contigs-'$K'mer_'$N'.txt' 'Kmer_Exon2_contigs-'$K'mer_'$N'.txt_merged'
mv 'Kmer_Exon2-3_contigs-'$K'mer_'$N'.txt' 'Kmer_Exon2-3_contigs-'$K'mer_'$N'.txt_merged'
mv 'Kmer_Exon3_contigs-'$K'mer_'$N'.txt' 'Kmer_Exon3_contigs-'$K'mer_'$N'.txt_merged'

python $script_assembler/Step3_Compare_for_rev_V2.py 'Kmer_Exon2_contigs-'$K'mer_'$N'.txt_merged' $K
python $script_assembler/Step3_Compare_for_rev_V2.py 'Kmer_Exon2-3_contigs-'$K'mer_'$N'.txt_merged' $K
python $script_assembler/Step3_Compare_for_rev_V2.py 'Kmer_Exon3_contigs-'$K'mer_'$N'.txt_merged' $K

### Rev ###

cd $result_tmp
python $script_assembler/Kmer_assembly4_V6.py DLA_Exon2_reverseStrandReads.fa-pair DLA_Exon2_forwardStrandReads.fa-pair $K $N 3 > Test4/'Kmer_Exon2_contigs-'$K'mer_'$N'-REV.txt'
python $script_assembler/Kmer_assembly4_V6.py DLA_Exon2-3_reverseStrandReads.fa-pair DLA_Exon2-3_forwardStrandReads.fa-pair $K $N 3 > Test4/'Kmer_Exon2-3_contigs-'$K'mer_'$N'-REV.txt'
python $script_assembler/Kmer_assembly4_V6.py DLA_Exon3_reverseStrandReads.fa-pair DLA_Exon3_forwardStrandReads.fa-pair $K $N 3 > Test4/'Kmer_Exon3_contigs-'$K'mer_'$N'-REV.txt'

cd $result_tmp/Test4
mv 'Kmer_Exon2_contigs-'$K'mer_'$N'-REV.txt' 'Kmer_Exon2_contigs-'$K'mer_'$N'-REV.txt_merged'
mv 'Kmer_Exon2-3_contigs-'$K'mer_'$N'-REV.txt' 'Kmer_Exon2-3_contigs-'$K'mer_'$N'-REV.txt_merged'
mv 'Kmer_Exon3_contigs-'$K'mer_'$N'-REV.txt' 'Kmer_Exon3_contigs-'$K'mer_'$N'-REV.txt_merged'

python $script_assembler/Step3_Compare_for_rev_V2.py 'Kmer_Exon2_contigs-'$K'mer_'$N'-REV.txt_merged' $K
python $script_assembler/Step3_Compare_for_rev_V2.py 'Kmer_Exon2-3_contigs-'$K'mer_'$N'-REV.txt_merged' $K
python $script_assembler/Step3_Compare_for_rev_V2.py 'Kmer_Exon3_contigs-'$K'mer_'$N'-REV.txt_merged' $K


assembly_merge=$result_tmp'/Test4/Merge2'
mkdir $assembly_merge
cd $assembly_merge

cat $result_tmp/Test4/*merged-Match > sample_contigs.fa # Merge all contigs
python $script_genotyper/Script0.5_merge_and_rename_contigs.py sample_contigs.fa 5 trim # Contigs dedup and trimming

rm -r tmp
cat sample_contigs.fa_merged_and_renamed $ref/DLA_Nucleotide_DB_AllExons.fa > tmp # contig classification (DLA-88*Norm DLA-88*50X DLA-12 DLA-64)

clustalo --force -i tmp -o sample-assembly.clustalO \
 --full --distmat-out=distance_matrix.txt --percent-id
rm tmp

# classify norm and 50X contigs
python $script_genotyper/Analyze_ClustalO_DistMatrix.py distance_matrix.txt sample_contigs.fa_merged_and_renamed

# merge D88-NORM, D12 as Norm
cat DLA88-norm_contigs.txt DLA88-50X_contigs.txt DLA12_contigs.txt > contigs.txt

# BWA mapping and extract read pairs
bwa index -a is sample_contigs.fa_merged_and_renamed

python $script_genotyper/Script9_fasta_to_fastq.py $result_tmp/DLA_Exon2-3_forwardStrandReads.fa $result_tmp/DLA_Exon2-3_reverseStrandReads.fa E

time bwa aln sample_contigs.fa_merged_and_renamed $result_tmp/DLA_Exon2-3_forwardStrandReads.fq > DLA_Exon2-3_forwardStrandReads.sai
time bwa aln sample_contigs.fa_merged_and_renamed $result_tmp/DLA_Exon2-3_reverseStrandReads.fq > DLA_Exon2-3_reverseStrandReads.sai

time bwa sampe sample_contigs.fa_merged_and_renamed  \
DLA_Exon2-3_forwardStrandReads.sai \
DLA_Exon2-3_reverseStrandReads.sai \
$result_tmp/DLA_Exon2-3_forwardStrandReads.fq \
$result_tmp/DLA_Exon2-3_reverseStrandReads.fq \
> sample_contigs.sam

cat sample_contigs.sam | grep -v '^@' | awk '$2%16==3 && $3 != "*" && $6 !~ /H/ && $6 !~ /S/ && $6 !~ /D/ && $6 !~ /I/ {print $0}' | grep -w 'XM:i:0' > sample_contigs_BWA_reads.sam
gawk '(and(16, $2))' sample_contigs_BWA_reads.sam > sample_contigs_BWA_reverseStrandReads.sam
gawk '(! and(16,$2))' sample_contigs_BWA_reads.sam > sample_contigs_BWA_forwardStrandReads.sam

python $script_assembler/Get_Pairwise_reads_V2.py sample_contigs_BWA_forwardStrandReads.sam sample_contigs_BWA_reverseStrandReads.sam 

# individual alignment and completeness checking
rm -r tmp
mkdir tmp

python $script_genotyper/Script5.1_Divide_consensus_sequences.py \
sample_contigs.fa_merged_and_renamed \
$assembly_merge/tmp

cd $assembly_merge/tmp
for i in Align*
do
        cat $i $script_genotyper/Db_Exon23_Nuc_index_reference.fa > foo.fa
        clustalo --force -i foo.fa -o $i'.clustalO'
        rm foo.fa
        python $script_genotyper/Script10_Extract_nuc_from_clustalO2.py $i'.clustalO'
done

for i in Align_*.clustalO_Exon2N3
do
python $script_genotyper/Summerize_coverage.py $i >> complete_contigs.txt
done

while read line
do
python $script_genotyper/classify_contigs.py $line'.clustalO_Exon2N3'
done < complete_contigs.txt

# Adjust mapping positions
cd tmp
ls *_50XDLA | sed 's/_50XDLA//g' > 50XDLA_lst.txt

cd $assembly_merge
python $script_genotyper/Script11_Adjust_sam_positions.py sample_contigs_BWA_reverseStrandReads.sam-pair $assembly_merge/tmp/ tmp/50XDLA_lst.txt
python $script_genotyper/Script11_Adjust_sam_positions.py sample_contigs_BWA_forwardStrandReads.sam-pair $assembly_merge/tmp/ tmp/50XDLA_lst.txt


# complete and incomplete contigs
python $script_genotyper/Script21_complete_incomplete_contig_classification.py \
$assembly_merge/tmp/complete_contigs.txt \
$assembly_merge/tmp \
contigs.txt \
DLA64_contigs.txt


python $script_genotyper/Script5_Contig_translation4.py DLA64_contigs_complete.txt

cat DLA64_contigs_complete.txt_ProteinSeq $script_genotyper/Db_Exon23_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o DLA64_contigs_complete.txt_ProteinSeq.clustalO
rm foo.fa

python $script_genotyper/Script6_Extract_3HVR_from_clustalO3.py DLA64_contigs_complete.txt_ProteinSeq.clustalO

## DLA-79, 12-18-22
# complete and incomplete contigs
python $script_genotyper/Script21_complete_incomplete_contig_classification.py \
$assembly_merge/tmp/complete_contigs.txt \
$assembly_merge/tmp \
contigs.txt \
DLA79_contigs.txt


python $script_genotyper/Script5_Contig_translation4.py DLA79_contigs_complete.txt

cat DLA79_contigs_complete.txt_ProteinSeq $script_genotyper/Db_Exon23_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o DLA79_contigs_complete.txt_ProteinSeq.clustalO
rm foo.fa

python $script_genotyper/Script6_Extract_3HVR_from_clustalO3.py DLA79_contigs_complete.txt_ProteinSeq.clustalO


# analyze normal and 50X seperately
mkdir nonDLA64
cp contigs_incomplete.txt nonDLA64

### normal DLA ###
cd $assembly_merge/nonDLA64

##clustalo --force -i norm_contigs.txt -o sample-assembly.clustalO
clustalo --force -i contigs_incomplete.txt -o sample-assembly.clustalO

# create reference consensus sequence and identify polymorphic sites
python $script_genotyper/Script1_ConsensusSeq_polymorphicSite_identification_ClustalO5.py sample-assembly.clustalO 0.75

# Align reference consensus with DLA database representative alleles
cat sample-assembly.clustalO_ConsensusSequence $script_genotyper/Db_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o sample-assembly.clustalO_ConsensusSequence_index
rm foo.fa

# correct clustalO alignment
python $script_genotyper/Script1.5_Check_alignment_and_correct_ClustalO.py \
sample-assembly.clustalO_ConsensusSequence_index \
sample-assembly.clustalO_ConsensusSequence \
sample-assembly.clustalO

# create reference consensus sequence and identify polymorphic sites
python $script_genotyper/Script1_ConsensusSeq_polymorphicSite_identification_ClustalO5.py sample-assembly.clustalO_corrected 0.75

# Align reference consensus with DLA database representative alleles
cat sample-assembly.clustalO_corrected_ConsensusSequence $script_genotyper/Db_index_reference.fa > foo.fa

if [[ $(grep '>' foo.fa | wc -l) -eq 1 ]]
then
cp foo.fa sample-assembly.clustalO_corrected_ConsensusSequence_index
else
clustalo --force -i foo.fa -o sample-assembly.clustalO_corrected_ConsensusSequence_index
fi

rm foo.fa

python $script_genotyper/Script2_Identify_polymorphicSite_linkages_with_sam11.py \
$assembly_merge/sample_contigs_BWA_forwardStrandReads.sam-pair_adjusted \
$assembly_merge/sample_contigs_BWA_reverseStrandReads.sam-pair_adjusted \
sample-assembly.clustalO_corrected_ConsensusSequence_index \
sample-assembly.clustalO_corrected_PolymorphicSites \
sample-assembly.clustalO_corrected_ContigLinkages \
sample-assembly.clustalO_corrected_ContigGroupsAndOutliers \
60

python $script_genotyper/classify_contigs_completeness.py \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigLinkages \
$assembly_merge/sample_contigs.fa_merged_and_renamed \
$assembly_merge/tmp/complete_contigs.txt \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers

python $script_genotyper/Script3_Merge_contigs_and_generate_consensus3.py \
sample-assembly.clustalO_corrected \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete \
$assembly_merge/sample_contigs_BWA_forwardStrandReads.sam-pair_adjusted \
$assembly_merge/sample_contigs_BWA_reverseStrandReads.sam-pair_adjusted \
sample-assembly.clustalO_corrected_ConsensusSequence_index

python $script_genotyper/Script4_ConsensusSequenceExtension_withPairedReads4.py \
$assembly_merge/sample_contigs_BWA_forwardStrandReads.sam-pair_adjusted \
$assembly_merge/sample_contigs_BWA_reverseStrandReads.sam-pair_adjusted \
sample-assembly.clustalO_corrected_ConsensusSequence_index \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs \
sample-assembly.clustalO_corrected_ConsensusSequence_index_PolymorphicSites_Linkage \
20


########## genotyping ###########
# rename extended contigs
python $script_genotyper/Script3_Rename_extended_contigs.py \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads \
ExtendedAlign

# extract core exons from extended contigs
rm -r tmp
mkdir tmp

python $script_genotyper/Script5.1_Divide_consensus_sequences.py \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed \
tmp/
cd tmp/
for i in ExtendedAlign*
do
        cat $i $script_genotyper/Db_Exon23_Nuc_index_reference.fa > foo.fa
        clustalo --force -i foo.fa -o $i'.clustalO'
        rm foo.fa
        python $script_genotyper/Script10_Extract_nuc_from_clustalO2.py $i'.clustalO'
done

for i in ExtendedAlign_*.clustalO_Exon2N3
do
python $script_genotyper/Summerize_coverage.py $i >> complete_contigs.txt
done

while read line
do
python $script_genotyper/classify_contigs.py $line'.clustalO_Exon2N3'
done < complete_contigs.txt

# extension validation
cd $assembly_merge/nonDLA64

if [ -f *clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed ]
then
if [[ $(grep '>' *clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed | wc -l) -ge 2 ]]
then

mkdir $assembly_merge/nonDLA64/Extended_contigs_validation
cp *clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed $assembly_merge/nonDLA64/Extended_contigs_validation

# BWA mapping
cd $assembly_merge/nonDLA64/Extended_contigs_validation
bwa index -a is sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed

bwa aln sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed $result_tmp/DLA_Exon2-3_forwardStrandReads_QualityFiltered.fq > DLA_Exon2-3_forwardStrandReads_QualityFiltered.sai
bwa aln sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed $result_tmp/DLA_Exon2-3_reverseStrandReads_QualityFiltered.fq > DLA_Exon2-3_reverseStrandReads_QualityFiltered.sai

bwa sampe sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed \
DLA_Exon2-3_forwardStrandReads_QualityFiltered.sai \
DLA_Exon2-3_reverseStrandReads_QualityFiltered.sai \
$result_tmp/DLA_Exon2-3_forwardStrandReads_QualityFiltered.fq \
$result_tmp/DLA_Exon2-3_reverseStrandReads_QualityFiltered.fq \
> sample_contigs.sam

cat sample_contigs.sam | grep -v '^@' | awk '$2%16==3 && $3 != "*" && $6 !~ /H/ && $6 !~ /S/ && $6 !~ /D/ && $6 !~ /I/{print $0}' | grep -w 'XM:i:0' > sample_contigs_BWA_reads.sam
gawk '(and(16, $2))' sample_contigs_BWA_reads.sam > sample_contigs_BWA_reverseStrandReads.sam
gawk '(! and(16,$2))' sample_contigs_BWA_reads.sam > sample_contigs_BWA_forwardStrandReads.sam

python $script_assembler/Get_Pairwise_reads_V2.py sample_contigs_BWA_forwardStrandReads.sam sample_contigs_BWA_reverseStrandReads.sam

# Adjust mapping positions
cd $assembly_merge/nonDLA64/tmp/
ls *_50XDLA | sed 's/_50XDLA//g' > 50XDLA_lst.txt

cd $assembly_merge/nonDLA64/Extended_contigs_validation
python $script_genotyper/Script11_Adjust_sam_positions.py sample_contigs_BWA_forwardStrandReads.sam-pair $assembly_merge/nonDLA64/tmp/ $assembly_merge/nonDLA64/tmp/50XDLA_lst.txt
python $script_genotyper/Script11_Adjust_sam_positions.py sample_contigs_BWA_reverseStrandReads.sam-pair $assembly_merge/nonDLA64/tmp/ $assembly_merge/nonDLA64/tmp/50XDLA_lst.txt

if [[ $(grep '>' sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed | wc -l) -eq 1 ]]
then
cp sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed sample-assembly.clustalO
else
clustalo --force -i sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed -o sample-assembly.clustalO
fi

# create reference consensus sequence and identify polymorphic sites
python $script_genotyper/Script1_ConsensusSeq_polymorphicSite_identification_ClustalO5.py sample-assembly.clustalO 0.75

# Align reference consensus with DLA database representative alleles
cat sample-assembly.clustalO_ConsensusSequence $script_genotyper/Db_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o sample-assembly.clustalO_ConsensusSequence_index
rm foo.fa

# correct clustalO alignment
python $script_genotyper/Script1.5_Check_alignment_and_correct_ClustalO.py \
sample-assembly.clustalO_ConsensusSequence_index \
sample-assembly.clustalO_ConsensusSequence \
sample-assembly.clustalO

# create reference consensus sequence and identify polymorphic sites
python $script_genotyper/Script1_ConsensusSeq_polymorphicSite_identification_ClustalO5.py sample-assembly.clustalO_corrected 0.75

# Align reference consensus with DLA database representative alleles
cat sample-assembly.clustalO_corrected_ConsensusSequence $script_genotyper/Db_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o sample-assembly.clustalO_corrected_ConsensusSequence_index
rm foo.fa

# check if all contigs are identical and no polymorphic site identified
if [ -s sample-assembly.clustalO_corrected_PolymorphicSites ]
then

python $script_genotyper/Script2_Identify_polymorphicSite_linkages_with_sam11.py \
sample_contigs_BWA_forwardStrandReads.sam-pair_adjusted \
sample_contigs_BWA_reverseStrandReads.sam-pair_adjusted \
sample-assembly.clustalO_corrected_ConsensusSequence_index \
sample-assembly.clustalO_corrected_PolymorphicSites \
sample-assembly.clustalO_corrected_ContigLinkages \
sample-assembly.clustalO_corrected_ContigGroupsAndOutliers \
60

python $script_genotyper/Script14_extract_sequences_from_correctedContigLinkages.py \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigLinkages \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed

cp sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed_validated $assembly_merge/nonDLA64
cd $assembly_merge/nonDLA64

else
cd $assembly_merge/nonDLA64
cp sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed_validated

fi

else
cd $assembly_merge/nonDLA64
cp sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed_validated

fi

else
cd $assembly_merge/nonDLA64
fi



# Extract exon 2 n 3
rm -r tmp
mkdir tmp

python $script_assembler/Merge_contigs.py Incomplete_contigs.fa 80

cat sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigGroupsAndOutliers_incomplete_ConsensusSequencesWithUnusedContigs_extendedWithPairedReads_renamed_validated Incomplete_contigs_merged.fa > sample-assembly_extended_merged_contig.fa

python $script_genotyper/Script5.1_Divide_consensus_sequences.py \
sample-assembly_extended_merged_contig.fa \
tmp/

cd tmp/
for i in Align*
do
        cat $i $script_genotyper/Db_Exon23_Nuc_index_reference.fa > foo.fa
        clustalo --force -i foo.fa -o $i'.clustalO'
        rm foo.fa
        python $script_genotyper/Script10_Extract_nuc_from_clustalO2.py $i'.clustalO'
done

for i in Align_*.clustalO_Exon2N3
do
python $script_genotyper/Summerize_coverage.py $i >> complete_contigs.txt
done

while read line
do
python $script_genotyper/classify_contigs.py $line'.clustalO_Exon2N3'
done < complete_contigs.txt




# merge all norm DLA complete contigs
mkdir $assembly_merge/nonDLA64/combined
cd $assembly_merge/nonDLA64

cat Complete_contigs.fa $script_genotyper/Db_Exon23_Nuc_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o Complete_contigs.clustalO
rm foo.fa
python $script_genotyper/Script10_Extract_nuc_from_clustalO2.py Complete_contigs.clustalO



cd $assembly_merge/nonDLA64/combined
cat $assembly_merge/contigs_complete.txt ../tmp/*_normDLA ../tmp/*_50XDLA > sample_nonDLA64_contigs.fa

python $script_genotyper/Script0.5_merge_and_rename_contigs.py sample_nonDLA64_contigs.fa 5 not

if [[ $(grep '>' sample_nonDLA64_contigs.fa_merged_and_renamed | wc -l) -eq 1 ]]
then
cp sample_nonDLA64_contigs.fa_merged_and_renamed sample-assembly.clustalO
else
clustalo --force -i sample_nonDLA64_contigs.fa_merged_and_renamed -o sample-assembly.clustalO
fi

python $script_genotyper/Script1_ConsensusSeq_polymorphicSite_identification_ClustalO5.py sample-assembly.clustalO 0.75

cat sample-assembly.clustalO_ConsensusSequence $script_genotyper/Db_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o sample-assembly.clustalO_ConsensusSequence_index
rm foo.fa


# correct clustalO alignment
python $script_genotyper/Script1.5_Check_alignment_and_correct_ClustalO.py \
sample-assembly.clustalO_ConsensusSequence_index \
sample-assembly.clustalO_ConsensusSequence \
sample-assembly.clustalO

# create reference consensus sequence and identify polymorphic sites
python $script_genotyper/Script1_ConsensusSeq_polymorphicSite_identification_ClustalO5.py sample-assembly.clustalO_corrected 0.75

# Align reference consensus with DLA database representative alleles
cat sample-assembly.clustalO_corrected_ConsensusSequence $script_genotyper/Db_index_reference.fa > foo.fa

if [[ $(grep '>' foo.fa | wc -l) -eq 1 ]]
then
cp foo.fa sample-assembly.clustalO_corrected_ConsensusSequence_index
else
clustalo --force -i foo.fa -o sample-assembly.clustalO_corrected_ConsensusSequence_index
fi

rm foo.fa

python $script_genotyper/Script2_Identify_polymorphicSite_linkages_with_sam11.py \
$assembly_merge/sample_contigs_BWA_forwardStrandReads.sam-pair_adjusted \
$assembly_merge/sample_contigs_BWA_reverseStrandReads.sam-pair_adjusted \
sample-assembly.clustalO_corrected_ConsensusSequence_index \
sample-assembly.clustalO_corrected_PolymorphicSites \
sample-assembly.clustalO_corrected_ContigLinkages \
sample-assembly.clustalO_corrected_ContigGroupsAndOutliers \
60

python $script_genotyper/Script1_extract_sequences_basedOn_linkages2.py \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigLinkages \
sample_nonDLA64_contigs.fa_merged_and_renamed \
Script2_polymorphicSites_nucleotide_frequencies

python $script_genotyper/Script5_Contig_translation4.py Correct_contigs.fa

cat Correct_contigs.fa_ProteinSeq $script_genotyper/Db_Exon23_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o Correct_contigs.fa_ProteinSeq.clustalO
rm foo.fa

python $script_genotyper/Script6_Extract_3HVR_from_clustalO3.py Correct_contigs.fa_ProteinSeq.clustalO

python $script_genotyper/Script2_Genotyping_for_conplete_prot_seq2.py \
Correct_contigs.fa_ProteinSeq.clustalO_3HVR \
Correct_contigs.fa_ProteinSeq.clustalO_Exon2N3 \
$ref/DLA_881264_HVR_Prot-dedup.fasta \
$ref/DLA_881264_Exon23_Prot-dedup.fasta \
Correct_contigs_scores.txt \
sample

# Correct percentage
if [[ $(grep '>' $assembly_merge/nonDLA64/combined/*_contigs.fa_merged_and_renamed | wc -l) -eq 1 ]]
then

cp $assembly_merge/nonDLA64/combined/*_contigs.fa_merged_and_renamed  Correct_contigs.fa

python $script_genotyper/Script5_Contig_translation4.py Correct_contigs.fa
cat Correct_contigs.fa_ProteinSeq $script_genotyper/Db_Exon23_index_reference.fa > foo.fa
clustalo --force -i foo.fa -o Correct_contigs.fa_ProteinSeq.clustalO
rm foo.fa

python $script_genotyper/Script6_Extract_3HVR_from_clustalO3.py Correct_contigs.fa_ProteinSeq.clustalO

## putting DLA-79 input
#

## DLA-79, 12-18-22
python $script_genotyper/Script2_Genotyping_for_complete_prot_seq-noPolymorphicSite.py \
Correct_contigs.fa_ProteinSeq.clustalO_3HVR \
Correct_contigs.fa_ProteinSeq.clustalO_Exon2N3 \
$assembly_merge/DLA64_contigs_complete.txt_ProteinSeq.clustalO_3HVR \
$assembly_merge/DLA64_contigs_complete.txt_ProteinSeq.clustalO_Exon2N3 \
$assembly_merge/DLA79_contigs_complete.txt_ProteinSeq.clustalO_3HVR \
$assembly_merge/DLA79_contigs_complete.txt_ProteinSeq.clustalO_Exon2N3 \
Script2_norm_sequences \
Script2_50X_sequences \
Script2_DLA64_sequences \
Script2_DLA79_sequences \
$ref/DLA_881264_HVR_Prot-dedup.fasta \
$ref/DLA_881264_Exon23_Prot-dedup.fasta \
sample

elif [[ $(grep '>' $assembly_merge/nonDLA64/combined/*_contigs.fa_merged_and_renamed | wc -l) -eq 0 ]]
then

#DLA-79,12-18-22
python $script_genotyper/Script2_Genotyping_for_complete_prot_seq-noPolymorphicSite.py \
Correct_contigs.fa_ProteinSeq.clustalO_3HVR \
Correct_contigs.fa_ProteinSeq.clustalO_Exon2N3 \
$assembly_merge/DLA64_contigs_complete.txt_ProteinSeq.clustalO_3HVR \
$assembly_merge/DLA64_contigs_complete.txt_ProteinSeq.clustalO_Exon2N3 \
$assembly_merge/DLA79_contigs_complete.txt_ProteinSeq.clustalO_3HVR \
$assembly_merge/DLA79_contigs_complete.txt_ProteinSeq.clustalO_Exon2N3 \
Script2_norm_sequences \
Script2_50X_sequences \
Script2_DLA64_sequences \
Script2_DLA79_sequences \
$ref/DLA_881264_HVR_Prot-dedup.fasta \
$ref/DLA_881264_Exon23_Prot-dedup.fasta \
sample


else


python $script_genotyper/Script12_Correct_contig_percentage4.py \
sample-assembly.clustalO_corrected_ConsensusSequence_index_correctedContigLinkages \
Script2_polymorphicSites_nucleotide_frequencies \
Correct_contigs.fa_ProteinSeq \
Script2_norm_sequences \
Script2_50X_sequences

#dla-79, 12-18-22
python $script_genotyper/Script2_Genotyping_for_complete_prot_seq2.py \
Correct_contigs.fa_ProteinSeq.clustalO_3HVR \
Correct_contigs.fa_ProteinSeq.clustalO_Exon2N3 \
$assembly_merge/DLA64_contigs_complete.txt_ProteinSeq.clustalO_3HVR \
$assembly_merge/DLA64_contigs_complete.txt_ProteinSeq.clustalO_Exon2N3 \
$assembly_merge/DLA79_contigs_complete.txt_ProteinSeq.clustalO_3HVR \
$assembly_merge/DLA79_contigs_complete.txt_ProteinSeq.clustalO_Exon2N3 \
Script2_norm_sequences \
Script2_50X_sequences \
Script2_DLA64_sequences \
Script2_DLA79_sequences \
$ref/DLA_881264_HVR_Prot-dedup.fasta \
$ref/DLA_881264_Exon23_Prot-dedup.fasta \
Correct_contigs_percentageScores.txt \
sample

fi



# gene assignment
if [ -f "Genotyping_result.txt" ]
then

python $script_genotyper/Script15_genotypingResult_to_fasta.py Genotyping_result.txt
cat Genotyping_result.txt_newAllele.fa $ref/DLA_881264_Exon23_Prot-dedup.fasta > GeneAssignment.fa
clustalo --force -i GeneAssignment.fa -o GeneAssignment.clustalO --distmat-out=GeneAssignment.dist --percent-id --full
python $script_genotyper/Analyze_ClustalO_DistMatrix_geneAssignment.py GeneAssignment.dist
python $script_genotyper/Script16_combine_geneAssignment_info.py Genotyping_result.txt GeneAssignment.dist_Summary

fi


# contig filtering
result_filter=$assembly_merge'/nonDLA64/combined/Filtering'

mkdir $result_filter
cd $result_filter

python $script_genotyper/Script18_generateContigNuc.py \
../Genotyping_result_withGeneAssignment.txt \
$result_filter/../../../ \
$result_filter

bwa index -a is NonDLA64_WholeNucSeq.fa

time bwa aln NonDLA64_WholeNucSeq.fa $result_tmp/DLA_Exon2-3_forwardStrandReads.fq > DLA_Exon2-3_forwardStrandReads.sai
time bwa aln NonDLA64_WholeNucSeq.fa $result_tmp/DLA_Exon2-3_reverseStrandReads.fq > DLA_Exon2-3_reverseStrandReads.sai

time bwa sampe NonDLA64_WholeNucSeq.fa \
DLA_Exon2-3_forwardStrandReads.sai \
DLA_Exon2-3_reverseStrandReads.sai \
$result_tmp/DLA_Exon2-3_forwardStrandReads.fq \
$result_tmp/DLA_Exon2-3_reverseStrandReads.fq \
-n 1000 \
> sample_contigs.sam

cat sample_contigs.sam | grep -v '^@' | awk '$2%16==3 && $3 != "*" && $6 !~ /H/ && $6 !~ /S/ && $6 !~ /D/ && $6 !~ /I/ {print $0}' | grep -w 'XM:i:0' > sample_contigs_BWA_reads.sam

python $script_genotyper/Script19_polyMorphicSite_identification.py \
NonDLA64_E23NucSeqAlignment.fa

python $script_genotyper/Script19_pairReadsPolymorphicSiteValidation.py \
NonDLA64_E23NucSeqAlignment_PolymorphicSites \
NonDLA64_E23NucSeqAlignment_ContigLinkages \
sample_contigs_BWA_reads.sam

cp ../Genotyping_result_withGeneAssignment.txt .

## add DLA-79
python $script_genotyper/Script20_contig_outputSummary3.py \
Genotyping_result_withGeneAssignment.txt \
$result_filter \
$assembly_merge/DLA64_contigs_complete.txt \
$assembly_merge/DLA79_contigs_complete.txt \
$DLA88_cutoffDepth \
$DLA12_cutoffDepth \
$cutoffLength \
$result_filter/sample_contigs_BWA_reads.sam


# print results
if [ -f "Genotyping_result_withGeneAssignment.txt" ]
then

cp Genotyping_result_withGeneAssignment.txt $result

else

echo "Genotyping failed" > $result/Genotyping_result_withGeneAssignment.txt

fi

if [ -f "Genotyping_result_withGeneAssignment.stat" ]
then

cp Genotyping_result_withGeneAssignment.stat $result

else

echo "Genotyping failed" > $result/Genotyping_result_withGeneAssignment.stat

fi

# checkpoints
wc -l execution/sample_BWA_reads.sam > checkpoint.txt
wc -l execution/Test4/Merge2/nonDLA64/combined/sample_nonDLA64_contigs.fa >> checkpoint.txt


