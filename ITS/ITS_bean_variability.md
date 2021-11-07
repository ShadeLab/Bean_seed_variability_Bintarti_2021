# Raw Sequences Analysis of ITS gene from Bean variability
# Date: February 5th 2021
# Author: A. Fina Bintarti

## Analysis of ITS Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME

raw sequence data stored on HPCC
/mnt/research/ShadeLab/Sequence/raw_sequence/Bean_variability/20190828_Amplicon_PE250/

Moved/copy raw sequences (45 samples, 1 mock community sample, 1 negative control and 1 positive control of DNA extraction) to the working space:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/ITS_Bean_variability/rawreads

Renaming the fastq files name in order to make it more simple, i.e.:
300_1_1_1_S23_L001_R1_001.fastq.gz to A11_S23_L001_R1_001.fastq.gz
A11 = plant A, pod # 1, seed # 1

using the loop as follow:
```
for i in 300_1_1_*
do mv "$i" "${i/300_1_1_/A1}"
done
```
## All analysis results are stored on hpcc: 
"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/ITS_Bean_variability/test"

### Quality checking 
```
# count read numbers
for fastq in rawreads/*.fastq
do wc -l $fastq
done > reads_raw.counts

# produce reads quality graphs using FastQC
mkdir stats

cat rawreads/*R1_001.fastq > raw_reads_R1.fastq; cat rawreads/*R2_001.fastq > raw_reads_R2.fastq

module load FastQC/0.11.7-Java-1.8.0_162

fastqc raw_reads_R1.fastq raw_reads_R2.fastq -o stats && rm -rf raw_reads_R1.fastq raw_reads_R2.fastq
```
# Clustering

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

### 1. Merge paired end reads
```
mkdir mergedfastq

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs rawreads/*R1*.fastq -relabel @ -tabbedout merged_tabbed.txt -report merged_summary.txt -fastqout mergedfastq/merged.fastq

# output
10422061  Pairs (10.4M)
   9114049  Merged (9.1M, 87.45%)
   5676521  Alignments with zero diffs (54.47%)
   1211545  Too many diffs (> 5) (11.62%)
         0  Fwd tails Q <= 2 trimmed (0.00%)
         3  Rev tails Q <= 2 trimmed (0.00%)
     96467  No alignment found (0.93%)
         0  Alignment too short (< 16) (0.00%)
     31198  Staggered pairs (0.30%) merged & trimmed
    135.08  Mean alignment length
    363.67  Mean merged length
      0.48  Mean fwd expected errors
      0.83  Mean rev expected errors
      0.23  Mean merged expected errors
```
### 2. Check sequence quality of merged sequences using Usearch and Vsearch
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 mergedfastq/merged.fastq -output stats_eestats2_USEARCH.txt

# output
9114049 reads, max len 467, avg 363.7

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    9071248( 99.5%)    9086647( 99.7%)    9086671( 99.7%)
   100    8893704( 97.6%)    9071369( 99.5%)    9085824( 99.7%)
   150    8740724( 95.9%)    9034750( 99.1%)    9085218( 99.7%)
   200    8529159( 93.6%)    8979547( 98.5%)    9082400( 99.7%)
   250    8436636( 92.6%)    8946909( 98.2%)    9079617( 99.6%)
   300    8212971( 90.1%)    8858516( 97.2%)    9065557( 99.5%)
   350    7951807( 87.2%)    8696132( 95.4%)    8976709( 98.5%)
   400        403(  0.0%)        543(  0.0%)        634(  0.0%)
   450          1(  0.0%)          2(  0.0%)          4(  0.0%)
###############################################################

module load vsearch/2.9.1

vsearch -fastq_stats mergedfastq/merged.fastq -fastq_qmax 42 -log stats_results_VSEARCH.txt
```
### 3. Remove primer and adapters with cutadapt
```
################
CS1-ITS86F (fwd): 5’- GTGAATCATCGAATCTTTGAA – 3’
CS2-ITS4 (rev): 5’- TCCTCCGCTTATTGATATGC – 3’ 
Reverse complement of reverse primer: GCATATCAATAAGCGGAGGA
#################

module load cutadapt (Ver. 2.0)

cutadapt -g GTGAATCATCGAATCTTTGAA -a GCATATCAATAAGCGGAGGA -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o cut_merged.fastq mergedfastq/merged.fastq > cut_adpt_results.txt

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 cut_merged.fastq -output cutdapt_eestats2_USEARCH.txt

# output
9113711 reads, max len 425, avg 322.7

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    9047304( 99.3%)    9084922( 99.7%)    9085302( 99.7%)
   100    8917872( 97.9%)    9072085( 99.5%)    9085218( 99.7%)
   150    8716409( 95.6%)    9026843( 99.0%)    9083865( 99.7%)
   200    8574615( 94.1%)    8994870( 98.7%)    9082405( 99.7%)
   250    8431100( 92.5%)    8942451( 98.1%)    9075005( 99.6%)
   300    8154888( 89.5%)    8786235( 96.4%)    8990591( 98.6%)
   350        447(  0.0%)        581(  0.0%)        672(  0.0%)
   400          1(  0.0%)          2(  0.0%)          5(  0.0%)
```
### 4. Quality filter
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 cut_merged.fastq -output cut_merged.pre_filtered.eestats2.txt -length_cutoffs 100,400,10

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter cut_merged.fastq -fastq_minlen 150 -fastq_maxee 1 -fastaout cut_merged_filtered.fa -fastaout_discarded merged.no_filter.fa -fastqout cut_merged_filtered.fastq

# output
100.0% Filtering, 96.8% passed
   9113711  Reads (9.1M)                    
    262694  Discarded reads with expected errs > 1.00
   8822144  Filtered reads (8.8M, 96.8%)
```
### 5. Find the set of unique sequences (dereplication)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques cut_merged_filtered.fastq -fastaout derep_filtered_cut_merged.fasta -sizeout

#output
8822144 seqs, 599680 uniques, 332168 singletons (55.4%)
```
### 6. Open reference-based OTU picking (using UNITE_v.8.0 at 97% identity treshhold)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -usearch_global derep_filtered_cut_merged.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v.8.0/sh_refs_qiime_ver8_97_02.02.2019.fasta  -strand plus -uc ref_seqs.uc -dbmatched UNITE_reference.fasta -notmatched UNITE_failed_closed.fq

```
### 7. Sorting by size and de novo-based OTU picking on sequences that failed to hit reference
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64  -sortbysize UNITE_failed_closed.fq -fastaout sorted_UNITE_failed_closed.fq

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64  -cluster_otus sorted_UNITE_failed_closed.fq -minsize 2 -otus de_novo_otus.fasta -uparseout uparse_otus.txt -relabel OTU_

# output
100.0% 52 OTUs, 62 chimeras
```
### 8. Combine the seqs of de novo and reference-based OTU picking
```
cat UNITE_reference.fasta de_novo_otus.fasta > REP_seq.fna

# numbering the OTUs for CONSTAX input
/mnt/home/bintarti/python_scripts-master/fasta_number.py REP_seq.fna OTU_ > 01312020_NUMB_REP_seq.fasta
```
### 9. Mapping/Construct OTU table
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -usearch_global mergedfastq/merged.fastq -db 01312020_NUMB_REP_seq.fasta -strand plus -id 0.97 -uc OTU_map.uc -otutabout 01312020_OpenRef_OTU_table.txt

# output
9078959 / 9114049 mapped to OTUs (99.6%)
```
Taxonomic classification using CONSTAX

Please refer to how ‘Running CONSTAX on the MSU HPCC’ on lab guru: https://my.labguru.com/knowledge/documents/330. Here is the publication about CONSTAX tool https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1952-x and how to use it https://github.com/natalie-vandepol/compare_taxonomy. Note: CONSTAX uses python 2.7 to be able to run.
```
# check python version
python --version

# if you have python3 installed and want to swap to python2
conda create -n python2 python=2.7 anaconda
conda activate python2

# download "CONSTAX_hpcc.tar.gz" from the lab guru link above. Put and extract the file on your home directory on hpcc.
# open "CONSTAX_hpcc" directory and follow the instructions from the lab guru link above.
# use the file 'consensus_taxonomy.txt' in the 'outputs' directory as your taxonomy table.
```
### 11. Convert .txt file to .biom file (if you need .biom file!)
```
biom convert -i 01312020_OpenRef_OTU_table.txt -o 01312020_OpenRef_OTU_table.biom --table-type="OTU table" --to-json

# summarize OTU table
biom summarize-table -i 01312020_OpenRef_OTU_table.biom -o 01312020_OpenRef_OTU_table_sum.txt
```










