# Risor_Sponge_Microbiome

This readme file was generated on 2025-09-02 by Ingrid Bunholi

This repository contains the code and data to reproduce all tables and figures presented in Bunholi et al. *"Sponge host shapes the skin and gut microbiome of a resident goby"*

**GENERAL INFORMATION**

Title of Dataset: Risor_Sponge_Microbiome

Name: Ingrid Bunholi ORCID: 0000-0001-5489-276X Email: ingridbunholi@gmail.com

Institution: University of Texas at Austin Marine Science Institute (UTMSI); Address: 750 Channel View Dr, Port Aransas, TX 78373

**DETAILS**

The elements of this project include:

1.sequence data: raw sequences can be found under the XXXXX NCBI accession number

2.metadata: this folder includes all the metadata necessary to run the scripts

3.scripts: this folder includes all scripts to generate the cleaned metadatas, figures and summary statistics

4.output: this folder includes all outputs from the R-script - include data generation, figures and summary statistics


**DATA & FILE OVERVIEW**

*Main metadata generation*

- filtered-table-all.qza: qiime2 object with raw ASVs for all sample types
- taxonomy-all.qza: qiime2 objecte with raw taxonomy assingment after filtration in qiime2
- mdat_all.tsv: Raw sample metadata with additional information for each sample (e.g., depth of collection, speciment total length, lat/long)
- tophits-all.tsv: file with vsearch metrics generated in qiime2
- tree.nwk: phylogenetic tree with raw ASVs generated in qiime2
- asv.all.filtered.csv: Raw ASV table after filtration in qiime2. This file contains all sample types
- taxon.all.filtered.csv: Raw taxonomy assingment after filtration in qiime2.
- genus.filtered.50ab.original.csv: final ASV dataset after all filtration steps grouped by genus for all sample types.
- metadata_filtered_all_2000.csv: metadata only for samples with acceptable sequencing depth. 

*Scripts*

- microbiome_script.sh: includes qiime2 QC and taxonomy assingment
- RisorMicrobiome_Cleaning_Filtering.R: includes taxonomy curation and additional filtration (e.g., negative controls, sequencing depth)
- RisorMicrobiome_Diversity_Sharedness.R: includes alpha/beta diversity and sharedness among samples (Figure 1 and 2, and 4)
- RisorMicrobiome_Null_ModelR: includes functions and null model analysis (Figure 3)
- RisorMicrobiome_Matrix_and_model_ontogeny.R: includes pairwise analysis and bayesian model of relative abundance and body size (Figure 5)

**METHODOLOGICAL INFORMATION**

*Description of methods used for collection/generation of data:* 

- Where/when: Four patch reefs near Carrie Bow Cay, Belize (July 2023; <10 m).

- Hosts: Ircinia campana, I. felix, I. strobilina with resident Risor ruber.

- Field workflow: Photograph sponge with ruler, record depth, tag to avoid resampling (all tags removed post-study). Anesthetize residents by spraying 5:1 ethanol:clove oil around excurrent pores; collect via SCUBA.

- Fish sampling: Skin microbiome—sterile swab across one body side (3 passes). Gut microbiome—dissect entire GI tract (fish ~1.6 cm).

- Sponge sampling: Cut cross-sections for microbiome.

- Storage: Skin swabs, GI tracts, and sponge sections into sterile 2 mL tubes with preservation buffer.

- DNA & amplicons: QIAGEN DNeasy PowerSoil Pro extraction (minor mods). Amplified 16S rRNA V4 region with 515F/806R primers following Earth Microbiome Project protocols

*Methods for processing the data:* 

- Pipeline: QIIME2 v2024.5 for QC and ASV taxonomy.

- Comparators: (i) Elacatinus prochilos skin/gut from 23 fish (Xavier et al. 2019); (ii) six seawater samples from Carrie Bow Cay (Gantt et al. 2019). Analyzed alongside our R. ruber–Ircinia samples to test whether microbial overlap differs from another goby and from the environment.

- Taxonomy thresholds: ≥97% for genus/species; ≥85% for family–class; ≥70% for phylum. ASVs were annotated to the highest rank meeting these cutoffs.

For full details, see Materials & Methods in the manuscript.

*Instrument- or software-specific information needed to interpret the data:* we used the statistical software R version 4.5 for all data analysis and visualization. Necessary packages to run the scripts: tidyverse, ggplot2, purrr, taxize, vegan, stringr, qiime2R, phyloseq, reshape, ape, UpSetR, picante, MASS, circlize, patchwork, glmmTB, brms, and tidybayes.

**Please contact me for any issue or question (ingridbunholi@gmail.com)**

