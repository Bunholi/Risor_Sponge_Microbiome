#!/bin/bash

#16S microbiome Pipeline - modified - IVB

## 1 - Setting things up

### Connect to a Server

ssh yourserver@yourserver

### Create a folder to upload your sequences

cd microbiome-risor

### Activate the genomics environment

conda activate genomics

### Make metadata file in google sheet and validate it using Keemei chrome extension

### input metadata in terminal 
#### Create a new .tsv file in the microbio-risor folder and copy and paste from google sheets
#### Save it by pressing command + s

### Run python script to clean up placement of files (script is already in the genomics environmen). You will need to have in the same folder: "mdat_skin.tsv" "reads" "report" "stats"

mkmanifest.py #it will generate a manifest.tsv table that will be used to import data into QIIME2

### 2 - Examine read quality with FastQC

mkdir fastqc_raw-reads

fastqc reads*/*_R1_* reads*/*_R2_* -o fastqc_raw-reads #running fastqc program in all the reads that have the pattern matching *_R1_* and *_R2_*

#### If you use Visual Studio Code, you can download each .html file straight from your folder and then visualize it in https://view.qiime2.org/

### 3 - Actiave QIIME 2 environment and importa data into it

conda activate qiime2-amplicon-2024.5

#### obs 1 - if you are in a different folder, make sure to give the whole path to manifest.tsv

qiime tools import\
   --type 'SampleData[PairedEndSequencesWithQuality]'\
   --input-path manifest.tsv\
   --output-path demux\
   --input-format PairedEndFastqManifestPhred33V2

### 4 - Quality control - summarize

#### Check quality of reads to determine how much to truncate before the paired ends read are joined
#### It will always depend on the length of amplicon and quality of reads

qiime demux summarize\
   --i-data /home/users/ibunholi/risor_microbio/risorgut_sponge/dada2_250_220/demux.qza\
   --o-visualization /home/users/ibunholi/risor_microbio/risorgut_sponge/dada2_250_220/demux 

mkdir visuals

mv demux.qzv visuals/ 

#### If you use Visual Studio Code, you can download your demux.qzv file from the explorer side tab and visualize it at https://view.qiime2.org/

#### 5 - Cutadapt - to remove primers

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/demux-skin.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r CCGYCAATTYMTTTRAGTTT \
  --p-discard-untrimmed \
  --p-match-adapter-wildcards \
  --p-match-read-wildcards \
  --p-cores 25 \
  --o-trimmed-sequences /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/demux-skin-cutadapt.qza \
  --verbose

### 6 - Denoising data - DADA2

#### It will fix base errors, trim primers, remove low quality regions, and filter out chimeras
#### Make sure to run dada2 separately for each sequencing run
#### Before choosing trunc, take a look at demux.qzv (interactive quality plot)
#f 250 and r 220

qiime dada2 denoise-paired\
   --i-demultiplexed-seqs /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/demux-skin-cutadapt.qza \
   --p-trunc-len-f 250 \
   --p-trunc-len-r 220 \
   --p-trim-left-f 0 \
   --p-trim-left-r 0 \
   --p-n-threads 16 \
   --o-denoising-stats /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/dns-skin \
   --o-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/table-skin \
   --o-representative-sequences /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/rep-seqs-skin

#### Visualizing DADA2 results

##### metadata on denoising
qiime metadata tabulate \
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/dns-skin.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-new/skin_P1/dns-skin

qiime metadata tabulate \
   --m-input-file /home/users/ibunholi/risor_microbio/elacatinus_p/dada2_220_150/table.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/elacatinus_p/dada2_220_150/table2   

##### unique sequences across all samples
qiime feature-table tabulate-seqs \
   --i-data /home/users/ibunholi/risor_microbio/elacatinus_p/dada2_220_150/rep-seqs.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/elacatinus_p/dada2_220_150/rep-seqs

#Denoise done separately for each sample type

##### dns.qzv will show the number of reads that passed filtration and how many were merged. Based on that you can evaluate if you need to change trunc options.
##### If relatively few are failing to merge, it suggests that the selected truncation length was good (check what is the best accepted percentage)

#### Merging denoised data of each sample type ####

qiime feature-table merge \
    --i-tables /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/skin_P1/table-skin.qza \
    --i-tables /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/gut_sponge_P2/table-gut-sponge.qza \
    --i-tables /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/prochilos/table-prochilos.qza \
    --i-tables /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/swater/table-seawater.qza \
    --o-merged-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all.qza

qiime feature-table merge-seqs \
    --i-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/skin_P1/rep-seqs-skin.qza \
    --i-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/gut_sponge_P2/rep-seqs-gut-sponge.qza \
    --i-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/prochilos/rep-seqs-prochilos.qza \
    --i-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/swater/rep-seqs-seawater.qza \
    --o-merged-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-rep-seqs-all.qza

qiime metadata tabulate \
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all

### 7 - Taxonomic assignment

#### vsearch uses a fast heuristic based on strings shared by the query and target sequences in order to quickly identify similar sequences

#16S
qiime feature-classifier classify-consensus-vsearch\
   --i-query /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-rep-seqs-all.qza \
   --i-reference-reads /home/share/databases/qiime2-databases/SILVA_databases/silva-138-99-seqs.qza \
   --i-reference-taxonomy /home/share/databases/qiime2-databases/SILVA_databases/silva-138-99-tax.qza \
   --p-maxaccepts 5 \
   --p-top-hits-only \
   --p-query-cov 0.4 \
   --p-perc-identity 0.7 \
   --o-classification taxonomy-all \
   --o-search-results tophits-all \
   --p-threads 72

#tophits output will be used to filter out ASVs that don't match our threshold for each taxonomic level (later in R)

qiime metadata tabulate \
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy-all.qza \
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-rep-seqs-all.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-feature-rep-tax.qzv

#### Visualising vsearch results

##### taxonomy.qza is the main file 
qiime metadata tabulate\
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy-all.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy.qzv

qiime metadata tabulate\
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/tophits-all.qza\
   --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/tophits.qzv

##### taxa-barplot can be visualize at https://view.qiime2.org/
qiime taxa barplot \
   --i-table table.qza \
   --i-taxonomy taxonomy.qza\
   --o-visualization taxa-barplot\
   --m-metadata-file /home/users/ibunholi/risor_microbio/skin_sponge/mdat_skin_sponge.tsv

### OPTIONAL 8 - Exporting/converting .qza files to .tsv 
#this can be done OR .qza files can be open directly in R using package qiime2R (read_qza)

#### table.qza
qiime tools export \
    --input-path /home/users/ibunholi/risor_microbio/skin_sponge/dada2_250_220/table.qza \
    --output-path ./ #converting .qza to .biom file

biom convert -i /home/users/ibunholi/risor_microbio/elacatinus/feature-table.biom -o filtered-table-elacatinus.tsv --to-tsv #converting .biom to .tsv file

#### taxonomy .qza
qiime tools export \
    --input-path /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/tophits-all.qza\
    --output-path ./ #converting .qza to .biom file - it seems to be converting straight to .tsv file in the new qiime version

#biom convert -i feature-taxonomy.biom -o taxonomy_250_220_gut_skin_sponge.tsv --to-tsv #converting .biom to .tsv file 

####ADDITIONAL STEPS ####
####Important: these following steps can also be done in R directly using the taxonomy.tsv and table.tsv files

#### 9 - Diversity analysis

#### Filter out unassigned, bacteria only, archeae, and eukaryotes

# D_0__Bacteria and unassigned
qiime taxa filter-table\
   --i-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all.qza\
   --i-taxonomy /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy-all.qza\
   --p-exclude "Unassigned"\
   --p-mode exact\
   --o-filtered-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/no-unassigned-table

# Eukaryota
qiime taxa filter-table\
   --i-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all.qza \
   --i-taxonomy /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy-all.qza \
   --p-exclude "Eukaryota" \
   --o-filtered-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/no-eukaryota-table

# Archaea
qiime taxa filter-table\
   --i-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/no-eukaryota-table.qza \
   --i-taxonomy /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy-all.qza \
   --p-exclude "Archaea"\
   --o-filtered-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/no-archaea-table

#### Filtering out singletons 

qiime feature-table filter-features \
  --i-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all.qza \
  --p-min-frequency 2 \
  --o-filtered-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-table-all-nosingle.qza

#### Check how much was filtered

#qiime feature-table summarize\
 #  --i-table filtered-table.qza\
 #  --m-sample-metadata-file /home/users/ibunholi/risor_microbio/skin_sponge/mdat_skin_sponge.tsv\
 #  --o-visualization filtered-table

#qiime taxa barplot \
 #  --i-table filtered-table.qza\
 #  --i-taxonomy taxonomy.qza\
 #  --o-visualization filtered-taxa-barplot\
 #  --m-metadata-file /home/users/ibunholi/risor_microbio/gut_sponge/mdat_gut_sponge.tsv

#### Filter the sequences to reflect the new table and generate feature-seq-tax file

qiime feature-table filter-seqs\
   --i-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/filtered-table-all.qza \
   --i-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-rep-seqs-all.qza \
   --o-filtered-data /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/filtered-rep-seqs-all

qiime metadata tabulate \
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/taxonomy-all.qza \
   --m-input-file /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/filtered-rep-seqs-all.qza \
   --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/filtered-feature-rep-tax.qzv

#qiime feature-table tabulate-seqs\
#   --i-data filtered-rep-seqs.qza\
#   --o-visualization filtered-rep-seqs

#qiime feature-table summarize \
#  --i-table filtered-table.qza \
#  --o-visualization table.qzv \
#  --m-sample-metadata-file /home/users/ibunholi/risor_microbio/16S_sponge/dada2_250_220/mdat_sponge.tsv

#### Creating phylogenetic tree

qiime phylogeny align-to-tree-mafft-fasttree\
  --i-sequences /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/merged-rep-seqs-all.qza\
  --o-alignment /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/raw-aligned-rep-seqs.qza\
  --o-masked-alignment /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/raw-masked-aligned-rep-seqs.qza\
   --o-tree /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/raw-unrooted-tree.qza\
   --o-rooted-tree /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/raw-rooted-tree.qza\
 --p-n-threads 18

qiime tools export \
  --input-path /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/raw-rooted-tree.qza \
  --output-path /home/users/ibunholi/risor_microbio/raw-skin-gut-sponge-prochilos-sw-new

qiime feature-table summarize \
  --i-table /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/filtered-table-all.qza \
  --o-visualization /home/users/ibunholi/risor_microbio/skin-gut-sponge-prochilos-sw-new/table_summary.qzv


#### Creating rarefaction curve
#alpha-rarefaction - generate interactive alpha rarefaction curves by computing rarefactions between min_depth and max_depth
#This is important to make samples compared to each other (rarefy samples to the same sampling depth)

#qiime diversity alpha-rarefaction\
#   --i-table filtered-table.qza\
#   --i-phylogeny rooted-tree.qza   \
#   --p-max-depth 5000   \
#   --p-steps 100 \
#   --m-metadata-file mdat_gut_skin_sponge.tsv   \
#   --o-visualization alpha-rarefaction.qzv

