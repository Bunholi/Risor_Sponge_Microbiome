############################################################
################# Risor microbiome project #################
####### Cleaning and filtering outputs from QIIME2 #########
###########################################################
#Updated by: Ingrid Bunholi
#Date: 07/15/26

#Required packages
library(tidyverse)
library(dplyr)
library(readr)
library(purrr)
library(taxize)
library(ggplot2)
library(vegan)
library(stringr)
library(qiime2R)
library(phyloseq)
library(reshape)

#### 1 - Environment set up and data loading ####
setwd("~/Desktop/Risor_microbiome_paper/Risor_Sponge_Microbiome/metadata") #your working directory

# Metadata
meta.all <- read_tsv("mdat_all.tsv")

# ASV table after QIIME removal of Eukaryota and Archaea
asv.all <- read_qza("no-archaea-table.qza")

asv.all.df <- as.data.frame(asv.all$data) %>%
  rownames_to_column("Feature.ID") %>%
  mutate(across(-Feature.ID, as.numeric))

#Representative sequences
seq.all <- read_qza("merged-rep-seqs-all.qza")

#Taxonomy
taxon.all <- read_qza("taxonomy-all.qza")

taxon.all.df <- as.data.frame(taxon.all$data) %>%
  separate(
    Taxon,
    into = c(
      "Domain", "Phylum", "Class", "Order",
      "Family", "Genus", "Species"
    ),
    sep = ";\\s*",
    fill = "right"
  ) %>%
  mutate(
    Domain = str_remove(Domain, "^d__"),
    Phylum = str_remove(Phylum, "^p__"),
    Class = str_remove(Class, "^c__"),
    Order = str_remove(Order, "^o__"),
    Family = str_remove(Family, "^f__"),
    Genus = str_remove(Genus, "^g__"),
    Species = str_remove(Species, "^s__")
  )

# Top-hit information
top <- read_tsv("tophits-all.tsv") %>%
  filter(pident >= 70) %>%
  group_by(Feature.ID) %>%
  slice_max(pident, n = 1, with_ties = FALSE) %>%
  ungroup()

#### 2 - Remove singleton and contaminant ASVs ####
contaminant_asvs <- c(
  "4ab56f9b012674eb6981ba9d3ef36170",
  "7ac9047b8f42a496b1640c20ea18156a",
  "94cb34fb244cbbfaeaafc9599a8f986b",
  "ac164649e33237f799c1ffee4f344d5d")

blank_samples <- c(
  "IBELBlank2",
  "SBELBlank",
  "RBELBlank1")

asv.clean <- asv.all.df %>%
  filter(!Feature.ID %in% contaminant_asvs) %>% #Remove ASVs detected as blank contaminants
  filter(rowSums(across(-Feature.ID), na.rm = TRUE) > 1) %>% #Remove singleton ASVs across the entire dataset
  dplyr::select(-any_of(blank_samples)) #Remove blank samples from the abundance table

#### 3 - Remove low-depth samples ####
#Calculate sequencing depth for each sample
sample_read_sums <- colSums(
  asv.clean %>% dplyr::select(-Feature.ID)
)

#Retain samples with at least 5,000 reads
samples_to_keep <- names(
  sample_read_sums[sample_read_sums >= 5000]
)

asv.clean <- asv.clean %>%
  dplyr::select(Feature.ID, all_of(samples_to_keep))

meta.clean <- meta.all %>%
  filter(sample_id %in% samples_to_keep)

#### 4 - Taxonomy standardization threshold ####

tax.clean <- taxon.all.df %>%
  filter(Feature.ID %in% asv.clean$Feature.ID) %>%
  left_join(
    top %>% dplyr::select(Feature.ID, pident),
    by = "Feature.ID"
  ) %>%
  mutate(
    Species = if_else(
      !is.na(Species) & pident >= 97,
      Species,
      NA_character_
    ),
    
    Genus = if_else(
      !is.na(Genus) & pident >= 97,
      Genus,
      NA_character_
    ),
    
    Family = if_else(
      !is.na(Family) & pident >= 85,
      Family,
      NA_character_
    ),
    
    Order = if_else(
      !is.na(Order) & pident >= 85,
      Order,
      NA_character_
    ),
    
    Class = if_else(
      !is.na(Class) & pident >= 85,
      Class,
      NA_character_
    ),
    
    Phylum = if_else(
      !is.na(Phylum) & pident >= 70,
      Phylum,
      NA_character_
    )
  )

#### 5 - Removing unwanted taxa ####
taxonomy_text <- tax.clean %>%
  unite(
    "taxonomy_string",
    Domain:Species,
    sep = ";",
    remove = FALSE,
    na.rm = TRUE
  )

taxa_to_remove <- taxonomy_text %>%
  filter(
    str_detect(
      taxonomy_string,
      regex("chloroplast|mitochondria", ignore_case = TRUE)
    )
  ) %>%
  pull(Feature.ID)

asv.clean <- asv.clean %>%
  filter(!Feature.ID %in% taxa_to_remove)

tax.clean <- tax.clean %>%
  filter(Feature.ID %in% asv.clean$Feature.ID)

#### 6 - Build standardized phyloseq objects ####
#Abundance matrix: taxa as rows, samples as columns
ASV.all <- asv.clean %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

storage.mode(ASV.all) <- "numeric"

# Taxonomy matrix
TAX.all <- tax.clean %>%
  distinct(Feature.ID, .keep_all = TRUE) %>%
  dplyr::select(
    Feature.ID,
    Domain, Phylum, Class, Order,
    Family, Genus, Species
  ) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# Metadata
META.all <- meta.clean %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  column_to_rownames("sample_id")

#checking if names match
common_taxa <- intersect(
  rownames(ASV.all),
  rownames(TAX.all))

common_samples <- intersect(
  colnames(ASV.all),
  rownames(META.all))

ASV.all <- ASV.all[
  common_taxa,
  common_samples,
  drop = FALSE]

TAX.all <- TAX.all[
  common_taxa,
  ,
  drop = FALSE]

META.all <- META.all[
  common_samples,
  ,
  drop = FALSE]

#creating main object: master ASV object, all samples, including Domain-unassigned ASVs
physeq.asv.all <- phyloseq(
  otu_table(ASV.all, taxa_are_rows = TRUE),
  tax_table(TAX.all),
  sample_data(META.all)) 

#adding phylogenetic tree
qiime_tree <- read.tree("tree.nwk")

tree_taxa <- intersect(
  taxa_names(physeq.asv.all),
  qiime_tree$tip.label)

physeq.asv.all <- prune_taxa(
  tree_taxa,
  physeq.asv.all)

qiime_tree <- keep.tip(
  qiime_tree,
  tree_taxa)

phy_tree(physeq.asv.all) <- qiime_tree

physeq.asv.all #73153 ASVs and 183 samples -master ASV object, all samples, including Domain-unassigned ASVs
setdiff(
  taxa_names(physeq.asv.all),
  phy_tree(physeq.asv.all)$tip.label) #names match tip label

#overall metrics
data.frame(
  total_samples = nsamples(physeq.asv.all),
  total_ASVs = ntaxa(physeq.asv.all),
  total_reads = sum(sample_sums(physeq.asv.all)),
  mean_reads_per_sample = mean(sample_sums(physeq.asv.all)),
  median_reads_per_sample = median(sample_sums(physeq.asv.all)),
  bacterial_ASVs = sum(
    tax_table(physeq.asv.all)[, "Domain"] == "Bacteria",
    na.rm = TRUE
  ),
  unassigned_ASVs = sum(
    tax_table(physeq.asv.all)[, "Domain"] == "Unassigned",
    na.rm = TRUE
  )
)

#Bacterial-assigned ASV object
tax_df <- as.data.frame(
  tax_table(physeq.asv.all))

bacterial_taxa <- rownames(tax_df)[
  !is.na(tax_df$Domain) &
    tax_df$Domain == "Bacteria"]

physeq.asv.bacteria <- prune_taxa(
  bacterial_taxa,
  physeq.asv.all)
physeq.asv.bacteria #13780 ASVs and 183 samples

#Risor-sponge ASV object
rs_samples <- grepl(
  "^(SBEL|RBEL|IBEL)",
  sample_names(physeq.asv.bacteria))

physeq.asv.bacteria.rs <- prune_samples(
  rs_samples,
  physeq.asv.bacteria)

physeq.asv.bacteria.rs <- prune_taxa(
  taxa_sums(physeq.asv.bacteria.rs) > 0,
  physeq.asv.bacteria.rs)

physeq.asv.bacteria.rs #11145 ASVs and 131 samples

#### 7 - Saving phyloseq objects ####
#Create folder for processed objects
dir.create("processed_objects", showWarnings = FALSE)

physeq.asv.all #73153 ASVs and 183 samples
saveRDS(
  physeq.asv.all,
  "processed_objects/physeq.asv.all.rds")

physeq.asv.bacteria #13780 ASVs and 183 sample
saveRDS(
  physeq.asv.bacteria,
  "processed_objects/physeq.asv.bacteria.rds")

physeq.asv.bacteria.rs #11145 and 131 samples
saveRDS(
  physeq.asv.bacteria.rs,
  "processed_objects/physeq.asv.bacteria.rs.rds")
