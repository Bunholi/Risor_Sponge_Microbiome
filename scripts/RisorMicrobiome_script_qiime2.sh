#!/bin/bash

#16S microbiome Pipeline

## 1 - Setting things up

### Connect to UNH Server

ssh utexasibunholi@ron.sr.unh.edu

### Create a folder to upload your sequences

cd microbiome-risor

wget -r -np -R "index.html*" --http-user=user --http-password=YdVCYKNnMc https://cobb.sr.unh.edu/managed/240410_A01346_0130_BHWKFKDRX3_16Mer040324-IB-P2-risor-microbiome/reads #remember to change address accordingly

### Activate the genomics environment

conda activate genomics

### Rename folder for simplification

mv 240410_A01346_0130_BHWKFKDRX3_16Mer040324-IB-P2-risor-microbiome/ 16S_pl1

cd 16S_pl1

ls reads #most importante folder - all sequences will be there

### Make metadata file in google sheet and validate it using Keemei chrome extension
#### See mdat_skin.tsv as an example 

### input metadata in terminal 
#### Create a new .tsv file in the microbio-risor folder and copy and paste from google sheets
#### Save it by pressing command + s

### Run Joe's python script to clean up placement of files (script is already in the genomics environmen). You will need to have in the same folder: "mdat_skin.tsv" "reads" "report" "stats"

mkmanifest.py #it will generate a manifest.tsv table that will be used to import data into QIIME2

### 2 - Examine read quality with FastQC

mkdir fastqc_raw-reads

fastqc reads*/*_R1_* reads*/*_R2_* -o fastqc_raw-reads #running fastqc program in all the reads that have the pattern matching *_R1_* and *_R2_*

#### If you use Visual Studio Code, you can download each .html file straight from your folder and then visualize it in https://view.qiime2.org/

### 3 - Actiave QIIME 2 environment and importa data into it

conda activate qiime2-2022.2

#### obs 1 - if you are in a different folder, make sure to give the whole path to manifest.tsv
#### obs 2 - make sure your reads names are matching with the manifest.tsv - I had to remove "23_" from my reads name otherwise all of they would have the same name. I wrote a quick bash script to remove it (rmv_charac.sh)

qiime tools import\
   --type 'SampleData[PairedEndSequencesWithQuality]'\
   --input-path manifest.tsv\
   --output-path demux\
   --input-format PairedEndFastqManifestPhred33V2

### 4 - Quality control - summarize

#### Check quality of reads to determine how much to truncate before the paired ends read are joined
#### It will always depend on the length of amplicon and quality of reads

qiime demux summarize\
   --i-data demux.qza\
   --o-visualization demux 

mkdir visuals

mv demux.qzv visuals/

#### If you use Visual Studio Code, you can download your demux.qzv file from the explorer side tab and visualize it at https://view.qiime2.org/

### 5 - Denoising data - DADA2

#### It will fix base errors, trim primers, remove low quality regions, and filter out chimeras
#### Make sure to run dada2 separately for each sequencing run
#### Before choosing trunc, take a look at demux.qzv (interactive quality plot)
#### Trim will be set based on your primer size (in this case 512-926R)

qiime dada2 denoise-paired\
   --i-demultiplexed-seqs /home/users/utexasibunholi/microbio-risor/16S_pl2/demux.qza \
   --p-trunc-len-f 250 \
   --p-trunc-len-r 220 \
   --p-trim-left-f 19 \
   --p-trim-left-r 20 \
   --p-n-threads 16 \
   --o-denoising-stats /home/users/utexasibunholi/microbio-risor/16S_pl1/dada2_250_220/dns \
   --o-table /home/users/utexasibunholi/microbio-risor/16S_pl1/dada2_250_220/table \
   --o-representative-sequences /home/users/utexasibunholi/microbio-risor/16S_pl1/dada2_250_220/rep-seqs

#### Visualizing DADA2 results

##### metadata on denoising
qiime metadata tabulate \
   --m-input-file /home/users/utexasibunholi/microbio-risor/16S_pl2/dns.qza \
   --o-visualization /home/users/utexasibunholi/microbio-risor/16S_pl2/dns

##### unique sequences across all samples
qiime feature-table tabulate-seqs \
   --i-data /home/users/utexasibunholi/microbio-risor/16S_pl2/dada2_240_210/rep-seqs.qza \
   --o-visualization /home/users/utexasibunholi/microbio-risor/16S_pl2/dada2_240_210/rep-seqs

##### generate table of per-sample sequence counts - important file
qiime feature-table summarize\
   --i-table /home/users/utexasibunholi/microbio-risor/16S_pl2/dada2_240_210/table.qza \
   --m-sample-metadata-file /home/users/utexasibunholi/microbio-risor/16S_pl2/mdat.tsv \
   --o-visualization /home/users/utexasibunholi/microbio-risor/16S_pl2/dada2_250_220/table

##### dns.qzv will show the number of reads that passed filtration and how many were merged. Based on that you can evaluate if you need to change trunc options.
##### If relatively few are failing to merge, it suggests that the selected truncation length was good (check what is the best accepted percentage)

### 6 - Taxonomic assignment - using SILVA_databases for 16S

#### vsearch uses a fast heuristic based on strings shared by the query and target sequences in order to quickly identify similar sequences

qiime feature-classifier classify-consensus-vsearch\
   --i-query /home/users/utexasibunholi/microbio-risor/16S_pl2/dada2_250_220/rep-seqs.qza \
   --i-reference-reads /home/share/databases/SILVA_databases/silva-138-99-seqs.qza \
   --i-reference-taxonomy /home/share/databases/SILVA_databases/silva-138-99-tax.qza \
   --p-maxaccepts 5 \
   --p-query-cov 0.4 \
   --p-perc-identity 0.7 \
   --o-classification taxonomy \
   --p-threads 72

#### Visualising vsearch results

##### taxonomy.qza is the main file 
qiime metadata tabulate\
   --m-input-file taxonomy.qza\
   --o-visualization taxonomy.qzv

##### taxa-barplot can be visualize at https://view.qiime2.org/
qiime taxa barplot \
   --i-table table.qza \
   --i-taxonomy taxonomy.qza\
   --o-visualization taxa-barplot\
   --m-metadata-file mdat.tsv

### 7 - Exporting/converting .qza files to .tsv - for following analyses in R

#### table.qza
qiime tools export \
    --input-path /home/users/utexasibunholi/microbio-risor/16S_pl1/dada2_250_220/table.qza \
    --output-path ./ #converting .qza to .biom file

biom convert -i feature-table.biom -o table_250_220_skin.tsv --to-tsv #converting .biom to .tsv file

#### taxonomy .qza
qiime tools export \
    --input-path /home/users/utexasibunholi/microbio-risor/16S_pl1/dada2_250_220/taxonomy.qza \
    --output-path ./ #converting .qza to .biom file

biom convert -i feature-taxonomy.biom -o taxonomy_250_220_skin.tsv --to-tsv #converting .biom to .tsv file