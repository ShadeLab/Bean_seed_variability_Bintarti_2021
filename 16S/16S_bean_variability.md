# Raw Sequences Analysis of 16S rRNA gene from the bean microbiome variability study

## Analysis of 16S Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME1

raw sequence data stored on HPCC
/mnt/research/ShadeLab/Sequence/raw_sequence/Bean_variability/20190903_16S-V4_PE250/

Moved/copy raw sequences (47 samples, 1 mock community sample, 1 negative control and 1 positive control of DNA extraction) to the working space:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/16S_Bean_variability

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
"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/16S_Bean_variability/16S_97"

# Part I: Clustering

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

## 1) Merge Paired End Reads
```
# decompress the reads
gunzip *.gz

# make directory called "mergedfastq"
mkdir mergedfastq

# merge paired end reads
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs rawseq/*R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -tabbedout mergedfastq/merged.report.txt -alnout mergedfastq/merged_aln.txt

# -fastq_maxdiffs 10: Allow 10 max differences in overlap region

### Output ###

 6511826  Pairs (6.5M)
   5266146  Merged (5.3M, 80.87%)
   1979863  Alignments with zero diffs (30.40%)
   1225635  Too many diffs (> 10) (18.82%)
     20045  No alignment found (0.31%)
         0  Alignment too short (< 16) (0.00%)
         0  Exp.errs. too high (max=1.0) (0.00%)
      4084  Staggered pairs (0.06%) merged & trimmed
    247.35  Mean alignment length
    252.59  Mean merged length
      0.90  Mean fwd expected errors
      1.34  Mean rev expected errors
      0.14  Mean merged expected errors
```
## 2) Check Sequence Quality of Merged Seqs
```
mkdir fastq_info
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged.fq -output fastq_info/eestats.txt

### output ###
5266146 reads, max len 450, avg 252.6

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    5218901( 99.1%)    5262930( 99.9%)    5265661(100.0%)
   100    5172977( 98.2%)    5256742( 99.8%)    5265538(100.0%)
   150    5108536( 97.0%)    5245159( 99.6%)    5264546(100.0%)
   200    5012267( 95.2%)    5224626( 99.2%)    5263511( 99.9%)
   250    4824794( 91.6%)    5156222( 97.9%)    5256694( 99.8%)
   300        388(  0.0%)        947(  0.0%)       1441(  0.0%)
   350        255(  0.0%)        711(  0.0%)       1325(  0.0%)
   400        147(  0.0%)        537(  0.0%)       1172(  0.0%)
```
## 3) Filter and Truncate the Merged Seqs
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter mergedfastq/merged.fq -fastq_maxee 1 -fastq_trunclen 250 -fastqout filtered_merged.fq

### output ###
100.0% Filtering, 97.9% passed
   5266146  Reads (5.3M)                    
      2993  Discarded reads length < 250
    106931  Discarded reads with expected errs > 1.00
   5156222  Filtered reads (5.2M, 97.9%)
```
## 4) Dereplicate Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_merged.fq -fastqout uniques_filtered_merged.fastq -sizeout

### output ###
5156222 seqs, 226310 uniques, 120386 singletons (53.2%)
Min size 1, median 1, max 2335131, avg 22.78
```
## 5) Remove Singeltons
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize uniques_filtered_merged.fastq -fastqout nosigs_uniques_filtered_merged.fastq -minsize 2

### output ###
Sorting 105924 sequences
```
## 6) Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast nosigs_uniques_filtered_merged.fastq -centroids_fastq denoised_nosigs_uniques_filtered_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

### output ###
Seqs  105924 (105.9k)
  Clusters  1344
  Max size  3790435 (3.8M)
  Avg size  3746.9
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  723Mb
      Time  3.00s
Throughput  35.3k seqs/sec.
```
## 7) Closed Reference-based OTU Picking Using SILVA_132 Database
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global denoised_nosigs_uniques_filtered_merged.fastq -id 0.97 -db /mnt/research/ShadeLab/WorkingSpace/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -strand plus -uc ref_seqs.uc -dbmatched SILVA_closed_reference.fasta -notmatchedfq failed_closed.fq

### output ###
100.0% Searching, 58.3% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the home directory.
Produce some output files - ref_seqs.uc (pre-clustered), SILVA_closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking
```
## 8) De novo OTU picking
```
# sort by size
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize failed_closed.fq -fastaout sorted_failed_closed.fq

# cluster de novo
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus sorted_failed_closed.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

### output ###
100.0% 69 OTUs, 172 chimeras
```
## 9) Combine the Rep Sets Between De novo and SILVA Reference-based OTU Picking
```
cat SILVA_closed_reference.fasta denovo_otus.fasta > FULL_REP_SET.fna
```
## 10) Map 'FULL_REP_SET.fna' Back to Pre-dereplicated Sequences and Make OTU Tables
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -usearch_global merged.fq -db FULL_REP_SET.fna -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

### output ###
5259262 / 5266146 mapped to OTUs (99.9%)
OTU_table.txt
OTU_jsn.biom
```
# Part II: Switch to QIIME 1.9.1

## 1) Convert OTU_table.txt to OTU_table.from_txt_json.biom
```
biom convert -i OTU_table.txt -o OTU_table.biom --table-type="OTU table" --to-json

### output ###
OTU_table.biom
```
## 2) Align sequences to the SILVA_132_QIIME_release database with PyNAST 
```
align_seqs.py -i FULL_REP_SET.fna -o alignment -t /mnt/home/bintarti/SILVA_132_QIIME_release/core_alignment/80_core_alignment.fna

### output ###
alignment/FULL_REP_SET_aligned.fasta
alignment/FULL_REP_SET_failures.fasta
alignment/FULL_REP_SET_log.txt
```
## 3) Filter failed alignment from OTU table
```
#Discard all OTUs listed in FULL_REP_SET_failures.fasta from OTU table

filter_otus_from_otu_table.py -i OTU_table.biom -o OTU_filteredfailedalignments.biom -e alignment/FULL_REP_SET_failures.fasta

#from FULL_REP_SET.fna file

filter_fasta.py -f FULL_REP_SET.fna -o FULL_REP_SET_filteredfailedalignments.fa -a alignment/FULL_REP_SET_aligned.fasta

### output ###
OTU_filteredfailedalignments.biom
FULL_REP_SET_filteredfailedalignments.fa
```
## 4) Assign taxonomy to the SILVA_132_QIIME_release database with UCLUST
```
assign_taxonomy.py -i FULL_REP_SET_filteredfailedalignments.fa -o taxonomy -r /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /mnt/home/bintarti/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt
#-r reference -> path to silva db 
# -t taxonomy

### output ###
taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments.txt
taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments.log
```
## 5) Add taxonomy to the OTU table
```
echo "#OTUID"$'\t'"taxonomy"$'\t'"confidence" > templine.txt

cat templine.txt taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments.txt >> taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments_header.txt

biom add-metadata -i OTU_filteredfailedalignments.biom -o OTU_table_tax.biom --observation-metadata-fp=taxonomy/FULL_REP_SET_filteredfailedalignments_tax_assignments_header.txt  --sc-separated=taxonomy --observation-header=OTUID,taxonomy

### output ###
OTU_table_tax.biom
```
## 6) Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i OTU_table_tax.biom -o OTU_table_tax_filt.biom -n D_4__Mitochondria,D_3__Chloroplast,Chlorophyta,Unassigned

#remove same Mito and Chloro sequences from RepSeqs file
filter_fasta.py -f FULL_REP_SET_filteredfailedalignments.fa -o FULL_REP_SET_filteredfailedalignments_rmCM.fa -b OTU_table_tax_filt.biom 

#summarize OTU table

#1.original

biom summarize-table -i OTU_table_tax.biom -o OTU_table_tax_sum.txt

#2.filtered

biom summarize-table -i OTU_table_tax_filt.biom -o OTU_table_tax_filt_sum.txt

#optional sanity check:  count seqs in new fasta, and check that it has fewer than original

#orginal

count_seqs.py -i FULL_REP_SET_filteredfailedalignments.fa

#filtered

count_seqs.py -i FULL_REP_SET_filteredfailedalignments_rmCM.fa

### output ###
OTU_table_tax_filt.biom
FULL_REP_SET_filteredfailedalignments_rmCM.fa
OTU_table_tax_sum.txt
OTU_table_tax_filt_sum.txt
```
## 7) Rarefaction – will not be performed in this study. We will use reads normalization using CSS method from the MetagenomeSeq package on R.
```
single_rarefaction.py -d 11137 -o single_rare.biom -i OTU_table_tax_filt.biom

biom summarize-table -i single_rare.biom -o single_rare_sum.txt

### output ###
single_rare.biom
single_rare_sum.txt
```
## 8) Summarize global taxonomic data
```
summarize_taxa.py -i OTU_table_tax_filt.biom -o taxa_sum

### output ###
taxa_sum/
```
## 9) Make phylogeny with FastTree
```
#First, clean alignment by omitting highly variable regions before tree building - will make tree building more efficient

filter_alignment.py -i alignment/FULL_REP_SET_aligned.fasta -o alignment/filtered_alignment

#make phylogeny and root tree

make_phylogeny.py -i alignment/filtered_alignment/FULL_REP_SET_aligned_pfiltered.fasta -o rep_set.tre -r tree_method_default

### output ###
alignment/filtered_alignment/FULL_REP_SET_aligned_pfiltered.fasta
rep_set.tre
```
## 10) Convert and add taxonomy
```
# 1. filtered OTU table
biom convert -i OTU_table_tax_filt.biom -o OTU_table_tax_filt.txt --header-key taxonomy --to-tsv

### output ###
OTU_table_tax_filt.txt

# 2. original OTU table
biom convert -i OTU_table_tax.biom -o OTU_table_tax.txt --header-key taxonomy --to-tsv


```
# Part III: Switch to R

## 10) Removing DNA contaminant from the OTU_table_tax_filt.txt file was conducted on R using microDecon package (McKnight et al. 2019)
```
### output ###
removing 18 contaminant OTUs = 255 OTUs
```
## 11) Removing OTUs that do not present in sample
```
### output ###
removing 44 zero OTUs = 211 OTUs
```
## 12) Reads normalization using CSS method from the metagenomeSeq package (Paulson et al. 2013)
```
A normalization method avoids biases due to uneven sequencing depth
### output ###
file = otu_norm.txt 
# move the file from the local computer to the hpcc (/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/16S_Bean_variability/test)
```
# Part IV: Switch to QIIME on hpcc

## 13) Convert "otu_norm.txt" to "otu_norm.biom"
```
biom convert -i otu_norm.txt -o otu_norm.biom --table-type="OTU table" --to-json

### output ###
otu_norm.biom
```
## 14) Prune the "rep_set.tre" based on the set of tip names from the "otu_norm.txt" (tips_to_keep.txt)
```
filter_tree.py -i rep_set.tre -t tips_to_keep.txt -o pruned.tre
```
## 15) Calculate PD_whole_tree (alpha diversity)
```
alpha_diversity.py -m PD_whole_tree -i otu_norm.biom -o PD -t pruned.tre
```














