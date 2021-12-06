# Raw Sequences Analysis of ITS gene from Bean variability
# Date: 
# Workflow author: This workflow was made  by Nejc Stopnisek and run as a job on hpcc on his WorkingSpace

## Analysis of ITS Miseq Data
Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). 

raw sequence data stored on HPCC
/mnt/research/ShadeLab/Sequence/raw_sequence/Bean_variability/20190828_Amplicon_PE250/

Moved/copy raw sequences (45 samples, 1 mock community sample, 1 negative control and 1 positive control of DNA extraction) to the working space:
/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/ITS_Bean_variability/rawreads

## All analysis results were moved from Nejc's WorkingSpace to mine on hpcc: 
"/mnt/research/ShadeLab/WorkingSpace/Bintarti/Bean_variability/SeedITS_nejc"

### 1. remove primers with cutadapt

```
files=/mnt/research/ShadeLab/Sequence/raw_sequence/Bean_variability/20190828_Amplicon_PE250_renamed

dir=/mnt/research/ShadeLab/WorkingSpace/Stopnisek/SeedITS

mkdir ${dir}/cut_trim
mkdir ${dir}/mergedFastq
mkdir ${dir}/results 

cd ${files}

for i in *R1_001.fastq
do
cutadapt -g ^GTGAATCATCGAATCTTTGAA -G ^TCCTCCGCTTATTGATATGC --discard-untrimmed -o ${dir}/cut_trim/${i%R1_001.fastq}_cut_R1.fastq -p ${dir}/cut_trim/${i%R1_001.fastq}_cut_R2.fastq ${files}/${i} ${files}/${i%R1_001.fastq}R2_001.fastq >> ${dir}/logFolder/cutadapt_primer_trimming_stats.txt 
done
```

### 2. Merge pairs

```
~/amoA_MiSeq_sequencing/usearch64 -fastq_mergepairs ${dir}/cut_trim/*R1.fastq -relabel @ -fastq_maxee 1 -fastq_maxns 0 -fastqout ${dir}/mergedFastq/merged.fq

```

### 3. Denoise and remove singletons

```
~/amoA_MiSeq_sequencing/usearch64 -fastx_uniques ${dir}/mergedFastq/merged.fq -fastaout ${dir}/results/uniques.fa -sizeout

~/amoA_MiSeq_sequencing/usearch64 -sortbysize ${dir}/results/uniques.fa -fastaout ${dir}/results/uniques_nosig.fa -minsize 2
```

### ZOTU (We used OTU at 97 % of identity threshold, so we skipped this!)
```
~/amoA_MiSeq_sequencing/usearch64 -unoise3 ${dir}/results/uniques_nosig.fa -zotus ${dir}/results/zotus.fa -tabbedout ${dir}/results/unoise_zotus.txt 

#Construct ZOTU table
sed -i 's/Zotu/ZOTU/g' ${dir}/results/zotus.fa

~/amoA_MiSeq_sequencing/usearch64 -otutab ${dir}/mergedFastq/merged.fq -zotus ${dir}/results/zotus.fa -otutabout ${dir}/results/ZOTU_table.txt -biomout ${dir}/results/ZOTU_jsn.biom 
```

### 4. 97% OTUs Clustering
```
~/amoA_MiSeq_sequencing/usearch64 -cluster_otus ${dir}/results/uniques_nosig.fa -minsize 2 -otus ${dir}/results/otus.fasta -uparseout ${dir}/results/uparse_otus.txt -relabel OTU_ --threads 28

~/amoA_MiSeq_sequencing/usearch64 -usearch_global ${dir}/mergedFastq/merged.fq -db ${dir}/results/otus.fasta -strand plus -id 0.97 -otutabout ${dir}/results/otu_table_ITS_UPARSE.txt
```

### 5. Taxonomy Classification using CONSTAX (CONSTAX version 1 only resulted RDP taxonomy so we used that for assigning taxonomy of our OTU table)
```
sh constax.sh

# download "CONSTAX_hpcc.tar.gz" from the lab guru link above. Put and extract the file on your home directory on hpcc.
# open "CONSTAX_hpcc" directory and follow the instructions from the lab guru link above.
# use the file 'consensus_taxonomy.txt' in the 'outputs' directory as your taxonomy table.
```
