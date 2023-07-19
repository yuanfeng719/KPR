# KPR
# **A Kmer-based paired-end read (KPR) *de novo* assembler and genotyper to genotype major histocompatibility complex class I (MHC-I) alleles for the dog** 

Kmer-based paired-end read (KPR) *de novo* assembler and genotyper take paired-end RNA-seq reads of an individual as input and output the MHC-I alleles of the individual.  The tools currently only work for the dog, consisting of the following steps: 1) mapping paired-RNA-seq reads to the canine MHC-I (i.e., DLA-I) reference alleles; 2) assemble paired-end RNA-seq reads mapped to the DLA-I regions into contigs *de novo*; and 3) genotyping the assembled contig.  

The software is described in detail at https://www.biorxiv.org/content/10.1101/2020.07.15.205559v2.   

## **Requirements**
The following package/software are required to run the KPR *de novo* assembler and genotyper  
```
  Python 2.7.14  
  SAMtools 1.9  
  Trimmomatic 0.39  
  BWA 0.7.17  
  Clustal-Omega 1.2.4  
```

## **Usage**
The package contains a Unix/Linux shell script (KPR_canine_MHC-I_assemblerNgenotyper.sh), a Python script directory and a DLA-I allele reference directory.  

### Steps to run the package  
1. Downloading the package (Linux/Unix platform is required)  

2. Modify the shell script (KPR_canine_MHC-I_assemblerNgenotyper.sh) by specifying the actual path to the Python script directory, the input RNA-seq fastq file directory, the result directory, and the reference directory, as well as the names of both the forward and reverse RNA-seq fastq files, as below.
```
script=''     # path to the KPR scripts  
data=''       # path to the RNA-seq paired-end fastq files  
result=''     # path to the output directory  
ref=''        # path to the KPR reference  

fastq1=''     # RNA-seq fastq file 1 (forward)  
fastq2=''     # RNA-seq fastq file 2 (reverse)  
```

3.  Adjust parameters: 1) the K-mer length, K; and 2) the number of assembly runs, N, based on the input RNA-seq data. Normally, we would suggest K=half of the RNA-seq read length and N>=1000. Larger N will take more time, but will yield more comprehensive genotyping.  
Other parameters include 3) DLA-88 depth; 4) DLA-12 depth; and 5) Identical length, which are used to identify chimeric and likely-chimeric assemblies. DLA-88 and DLA-12 depths are the minimum number of reads spanning polymorphic site-pairs. Identical length means the tolerable sequence length (nucleotides) shared by two assemblies.  
```
K=50     # K-mer length  
N=1000     # The number of assembly runs  

DLA88_cutoffDepth=15     # DLA-88 depth  
DLA12_cutoffDepth=38     # DLA-12 depth  
cutoffLength=270     # Tolerable identical sequence length shared by any two assemblies
```

4.  Required package/software loading:  the shell script contain lines for loading required package/software that are unique to the UGA (Sapelo2) platform.  If your platform uses a different approach of loading package/software, you will need to make corresponding changes in the script for lines below.  
```
module load SAMtools/1.9-GCC-8.3.0  
module load Trimmomatic/0.39-Java-1.8.0_144  
module load BWA/0.7.17-GCC-8.3.0  
module load SRA-Toolkit/2.10.8-centos_linux64  
module load Clustal-Omega/1.2.4-GCC-8.3.0  
```

5.  Once the shell script execution is finished, the result file, Genotyping_result_withGeneAssignment.stat, will appear in the $result directory you specified, containing the following columns: 
```
[1] Sample  
[2] Contig name  
[3] Relative expression level in the sample  
[4] Allele name (if match a known allele)  
[5] Allele group name (if match a known allele group based on hypervariable regions)  
[6] Gene name  
[7] Exons 2 and 3 protein sequence  
[8] Hypervariable region protein sequence  
[9] Exons 2 and 3 nucleotide sequence  
[10] Allele confidence  
[11] Actually expression level (Exons 2 and 3 reads assigned to this contig)  
[12] Relative expression level in the sample before adjustment (reallocate Chimeric expression levels)  
[13] Actually expression level before adjustment  
[14] Evidence of pairwise polymorphic site linkages (e.g. 89A_132C(116/265) means current contig carries “A” at position 89 and “C” at position 132, a total of 265 reads span these two positions while 116 reads support this A-C linkage)  
[15] Polymorphic site supporting read counts (e.g. 89A:200 means 200 reads carry “A” at position 89)  
[16] Read pair linkages (e.g. 89A;132C;213A;215C;216T;217G;219A(42) means 42 read pairs span polymorphic sites 89-132-213-215-216-217-219, with a linkage of A-C-A-C-T-G-A at these positions, respectively)  
```

### An example output line is provided below:  
```
Dog_ID  NormAlign_1     0.357251481896  DLA-88*006:01   DLA-88*006      DLA-88  GSHSLRYFYTSVSRPGRGDPRFIAVGYVDDTQFVRFDSDAATGRTEPRAPWVEQEGPEYWDPQTRTIKETAQLYRVDLDTLRGYYNQSEAGSHTRQTMYGCDLGPGGRLLRGYSQDAYDGADYIALNEDLRSYTAADTAAQITRRKWEAAGTAEHDRNYLETTCVEWLRRYLEMGKETLQRA     PQTRTIKETAQLYRVDGSHTRQTMYGCDLGPGGRLLRGYSQDTAEHDR        GCTCCCACTCCCTGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGACCCCCGCTTCATCGCCGTCGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCGGCCACTGGGAGGACGGAGCCGCGGGCGCCGTGGGTGGAGCAGGAGGGGCCGGAGTATTGGGACCCGCAGACGCGGACCATCAAGGAGACCGCACAGCTGTACCGAGTGGACCTGGACACCCTGCGCGGCTACTACAACCAGAGCGAGGCCGGGTCTCACACCCGCCAGACCATGTACGGCTGTGACCTGGGGCCCGGCGGGCGCCTCCTCCGCGGGTACAGTCAGGACGCCTACGACGGCGCCGATTACATCGCCCTGAACGAGGACCTGCGCTCCTACACCGCGGCGGACACGGCGGCGCAGATCACCCGGCGCAAGTGGGAAGCGGCAGGTACTGCAGAGCACGATAGGAACTACCTGGAGACGACGTGCGTGGAGTGGCTGCGGAGGTACCTGGAGATGGGGAAGGAGACGCTGCAGCGCGCAG       Known allele    29069.0 0.357251481896  29069.0 89A_132C(116/265),89A_213A(46/107),89A_215C(46/106), …	89A:200,132C:186,213A:157, …	89A(44),89A;132C(71),89A;132C;213A;215C;216T;217G(2), …  
```
