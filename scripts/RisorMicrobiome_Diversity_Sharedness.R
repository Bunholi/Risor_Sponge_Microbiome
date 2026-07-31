####################################################
######### Risor microbiome project #################
####### Diveristy analysis and sharednesss #########
####################################################
#Updated by: Ingrid Bunholi
#Date: 07/15/26

#Required packages
library(tidyverse)
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(vegan)
library(stringr)
library(qiime2R)
library(phyloseq)
library(reshape)
library(ape)
library(data.table)
library(BiocManager)
library(tibble)
library(UpSetR)
library(picante)
library(MASS)
library(eulerr)

#### Environment set up and data loading ####
setwd("~/Desktop/Risor_microbiome_paper/Risor_Sponge_Microbiome/metadata") #your working directory

#Loading processed phyloseq objects

#All samples and all ASVs
physeq.asv.all<- readRDS("processed_objects/physeq.asv.all.rds")
physeq.asv.all

# All samples and bacterial ASVs
physeq.asv.bacteria<- readRDS("processed_objects/physeq.asv.bacteria.rds")
physeq.asv.bacteria
table(sample_data(physeq.asv.bacteria)$type)

# Risor + sponge, bacterial ASVs
physeq.asv.bacteria.rs <- readRDS(
  "processed_objects/physeq.asv.bacteria.rs.rds")

#### Figure 1 C and D - Diversity metrics ####

physeq_fig1 <- physeq.asv.bacteria

physeq_fig1_qc <- prune_samples(
  sample_sums(physeq_fig1) >= 1000,
  physeq_fig1)

physeq_fig1_qc <- prune_taxa(
  taxa_sums(physeq_fig1_qc) > 0,
  physeq_fig1_qc)

physeq_fig1_qc
table(sample_data(physeq_fig1_qc)$type) 

#Calculate alpha diversity metrics for Shannon and Observed
alpha_metrics <- estimate_richness(physeq_fig1_qc, measures = c("Shannon", "Observed"))

#Combine alpha metrics with sample metadata
alpha_data <- cbind(alpha_metrics, as(sample_data(physeq_fig1_qc), "data.frame"))

alpha_data2<- alpha_data %>%
  tibble::rownames_to_column(var = "SampleID")

#Phylogenetic diversity - ASV level
otu <- as(otu_table(physeq_fig1_qc), "matrix") #isolating asv table
if(taxa_are_rows(physeq_fig1_qc)) { otu <- t(otu) }

tree <- phy_tree(physeq_fig1_qc) #isolating tree

pd_result <- pd(otu, tree, include.root = TRUE)

pd_result<- pd_result %>%
  tibble::rownames_to_column(var = "SampleID")

alpha_merged<- left_join(alpha_data2, pd_result, by = "SampleID") #merging PD with the other alpha diversity metrics

#Testing the effect of SAMPLE TYPE on ALPHA DIVERSITY using Linear Models

#Shannon
hist(alpha_merged$Shannon) #checking data distribution

lm_shannon<- lm(Shannon ~ type, data = alpha_merged) 
hist(resid(lm_shannon))
qqnorm(resid(lm_shannon)); qqline(resid(lm_shannon))

summary(lm_shannon)

#Observed
hist(alpha_merged$Observed) #checking data distribution

glm_nb <- glm.nb(Observed ~ type, data = alpha_merged) #glm with negative binomial distribution
summary(glm_nb)

#PD
hist(alpha_merged$PD)

glm_pd <- glm(PD ~ type, data = alpha_merged, family = Gamma(link = "log")) #glm with gamma distribution
summary(glm_pd)

##Plotting alpha diversity for ASV - sample type - Figure 1C##
alpha_long <- pivot_longer(alpha_merged, 
                           cols = c("Observed", "Shannon", "PD"), 
                           names_to = "metric", 
                           values_to = "value")

alpha_long <- alpha_long %>%
  mutate(type = factor(type,
      levels = c(
        "Risor gut",
        "Risor skin",
        "E.prochilos gut",
        "E.prochilos skin",
        "sponge",
        "seawater")))


alpha_plot <- ggplot(alpha_long, aes(x = type, y = value, fill = type, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) + 
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +  
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  labs(x = "Sample Type", y = "Alpha Diversity") +
  scale_fill_manual(values = c("Risor skin" = "#DDCC77", 
                               "sponge" = "#1D2D55", 
                               "Risor gut" = "#44AA99",
                               "E.prochilos skin" = "#B8887A",
                               "E.prochilos gut" = "darkslategrey", 
                               "seawater" = "royalblue"
                               )) +
  scale_color_manual(values = c("Risor skin" = "#DDCC77", 
                                "sponge" = "#1D2D55", 
                                "Risor gut" = "#44AA99",
                                "E.prochilos skin" = "#B8887A",
                                "E.prochilos gut" = "darkslategrey", 
                                "seawater" = "royalblue"))

alpha_plot +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Observed" ~ scale_y_continuous(limits = c(0, 1000)),
      metric == "PD"       ~ scale_y_continuous(limits = c(0, 45)),
      metric == "Shannon"  ~ scale_y_continuous()
    )
  )
ggsave("alpha_boxplot_ASV_Bac_Risor_Sponge_Prochilos_Seawater.pdf", alpha_plot , width = 3.5, height = 9, useDingbats = FALSE)

# Beta diversity - Figure 1D 

pcoa.asv<- ordinate(physeq_fig1_qc, "PCoA", "unifrac", weighted=FALSE) #Unifrac distance
pcoa.plot.asv<- plot_ordination(physeq_fig1_qc, pcoa.asv, color="type")
pcoa.plot.asv<- pcoa.plot.asv + geom_point(size=2)+
  theme_bw()+
  scale_colour_manual(values = c("Risor skin" = "#DDCC77", 
                                 "sponge" = "#1D2D55", 
                                 "Risor gut" = "#44AA99",
                                 "E.prochilos skin" = "#B8887A",
                                 "E.prochilos gut" = "darkslategrey", 
                                 "seawater" = "royalblue"))+
  labs(color="Sample Type")+
  stat_ellipse(
    aes(color = type),
    type = "norm",
    level = 0.95,
    linetype = 2,
    linewidth = 0.2,
    show.legend = FALSE
  )+
  theme(legend.position = "none")
pcoa.plot.asv

ggsave("PCoA_ASV_Bac_Risor_Sponge_Prochilos_Seawater.pdf", pcoa.plot.asv , width = 5.5, height = 5, useDingbats = FALSE)


#Checking whether unifrac differences are statistically different - PERMANOVA
unifrac_dist <- phyloseq::distance(
  physeq_fig1_qc,
  method = "unifrac",
  weighted = FALSE
)

metadata_rs <- data.frame(
  SampleID = sample_names(physeq_fig1_qc),
  type = sample_data(physeq_fig1_qc)$type
)

rownames(metadata_rs) <- metadata_rs$SampleID

set.seed(123)

adonis_res<-adonis2(
  unifrac_dist ~ type,
  data = metadata_rs,
  permutations = 999
)

print(adonis_res)


##Beta diversity at the genus level

physeq_fig1qc_genus <- physeq_fig1_qc %>%
  tax_glom(taxrank = "Genus", NArm = FALSE)

physeq_fig1qc_genus <- subset_taxa(
  physeq_fig1qc_genus,
  !is.na(Genus) &
    Genus != "" &
    Genus != "uncultured" &
    Genus != "Unassigned"
)

pcoa.gen<- ordinate(physeq_fig1qc_genus, "PCoA", "unifrac", weighted=FALSE) #Unifrac distance
pcoa.plot.gen<- plot_ordination(physeq_fig1qc_genus, pcoa.gen, color="type")
pcoa.plot.gen<- pcoa.plot.gen + geom_point(size=2)+
  theme_bw()+
  scale_colour_manual(values = c("Risor skin" = "#DDCC77", 
                                 "sponge" = "#1D2D55", 
                                 "Risor gut" = "#44AA99",
                                 "E.prochilos skin" = "#B8887A",
                                 "E.prochilos gut" = "darkslategrey", 
                                 "seawater" = "royalblue"))+
  labs(color="Sample Type")+
  stat_ellipse(
    aes(color = type),
    type = "norm",
    level = 0.68,
    linetype = 2,
    linewidth = 0.5,
    show.legend = FALSE
  )+
  theme(legend.position = "top")
pcoa.plot.gen

ggsave("PCoA_Genus_Bac_Risor_Sponge_Prochilos_Seawater.pdf", pcoa.plot.gen , width = 5.5, height = 5, useDingbats = FALSE)


#Checking whether unifrac differences are statistically different - PERMANOVA
unifrac_dist <- phyloseq::distance(
  physeq_fig1_qc,
  method = "unifrac",
  weighted = FALSE
)

metadata_rs <- data.frame(
  SampleID = sample_names(physeq_fig1_qc),
  type = sample_data(physeq_fig1_qc)$type
)

rownames(metadata_rs) <- metadata_rs$SampleID

set.seed(123)

adonis_res<-adonis2(
  unifrac_dist ~ type,
  data = metadata_rs,
  permutations = 999
)

print(adonis_res)

#Checking whether differences are due to within group dispersion
dispersion <- betadisper(
  unifrac_dist,
  group = metadata_rs$type
)

permutest(dispersion, permutations = 999)
anova(dispersion)


##Beta diversity sponge - supplemental material (Figure S2)
custom_colors <- c(
  "Ircinia campana" = "lightsalmon3",
  "Ircinia felix" = "brown",
  "Ircinia ruetzleri" = "green4",
  "Ircinia strobilina" = "gray20"
)

# Remove samples where sponge_spp is NA
physeq_fig1_qc.sponge <- subset_samples(physeq_fig1_qc, !is.na(sponge_spp))

physeq_fig1_qc.sponge@sam_data$sponge_spp <- factor(physeq_fig1_qc.sponge@sam_data$sponge_spp, 
                                                    levels = c("Ircinia campana", "Ircinia felix", "Ircinia ruetzleri", 
                                                               "Ircinia strobilina"))

pcoa.gen.sponge<- ordinate(physeq_fig1_qc.sponge, "PCoA", "unifrac", weighted=FALSE)

pcoa.gen.sponge<- plot_ordination(physeq_fig1_qc.sponge, pcoa.gen.sponge, color="sponge_spp")

pcoa.gen.plot.sponge<- pcoa.gen.sponge + geom_point(size=2)+
  theme_bw()+
  scale_colour_manual(values = custom_colors)+
  labs(color="Sponge species")+
  stat_ellipse(aes(group = sponge_spp), level = 0.95, linetype = 2)
pcoa.gen.plot.sponge

ggsave("pcoa_spongespp_ASV_new.pdf", pcoa.gen.plot.sponge, width = 8, height = 6, useDingbats = FALSE)


#Checking whether unifrac differences are statistically different
unifrac_dist <- phyloseq::distance(physeq_fig1_qc.sponge, method = "unifrac", weighted = FALSE)

sample_meta <- as(sample_data(physeq_fig1_qc.sponge), "data.frame")

#PERMANOVA using adonis
set.seed(123)  # For reproducibility of the permutation test
adonis_res_sponge <- adonis2(unifrac_dist ~ sponge_spp, data = sample_meta, permutations = 999)

print(adonis_res_sponge)

#### Figure 2 A and B - Shared quantity ####

physeq_fig2_asv <- physeq_fig1_qc

#Remove Cyanobacteria

tax_df <- as.data.frame(tax_table(physeq_fig2_asv))

is_cyanobacteria <- !is.na(tax_df$Phylum) &
  grepl("Cyanobacter", tax_df$Phylum, ignore.case = TRUE)

physeq_fig2_asv <- prune_taxa(
  taxa_names(physeq_fig2_asv)[!is_cyanobacteria],
  physeq_fig2_asv)

physeq_fig2_asv <- prune_taxa(
  taxa_sums(physeq_fig2_asv) > 0,
  physeq_fig2_asv)

physeq_fig2_asv #10692 ASVs and 121 samples (only 395 ASVs removed)

#checking if any samples fall below 1000 reads after filtering - no they did not
data.frame(
  type = sample_data(physeq_fig2_asv)$type,
  reads = sample_sums(physeq_fig2_asv)
) %>%
  group_by(type) %>%
  summarise(
    n = n(),
    min_reads = min(reads),
    median_reads = median(reads),
    n_below_1000 = sum(reads < 1000),
    .groups = "drop")

#Figure 2A - ASV level
#create a sample-by-sample binary matrix for venndiagram
asv_counts <- as(otu_table(physeq_fig2_asv), "matrix")

if (taxa_are_rows(physeq_fig2_asv)) {
  asv_counts <- t(asv_counts)}

metadata_fig2 <- data.frame(sample_data(physeq_fig2_asv))

metadata_fig2 <- metadata_fig2[
  rownames(asv_counts),
  ,
  drop = FALSE]

asv_binary <- asv_counts > 0 #rows are samples and columns are ASVs

#Separate by sample type 
groups <- c("Risor gut", "Risor skin", "sponge")

asv_sets <- setNames(
  lapply(groups, function(group_name) {
    
    group_samples <- metadata_fig2$type == group_name
    
    colnames(asv_binary)[
      colSums(asv_binary[group_samples, , drop = FALSE]) > 0
    ]
  }),
  groups
)

sapply(asv_sets, length)

#plot the Euler diagram
euler_asv <- eulerr::euler(asv_sets)

pdf(file = "Euler_2A.pdf", width = 4, height = 3)
plot(
  euler_asv,
  quantities = TRUE,
  labels = TRUE,
  fills = list(
    fill = c(
      "Risor gut" = "#44AA99",
      "Risor skin" = "#DDCC77",
      "sponge" = "#1D2D55"), alpha = 0.8))
dev.off()

#Euler - genus level
physeq_fig2_genus=tax_glom(physeq_fig2_asv, taxrank="Genus") #931 genera
physeq_fig2_genus<- subset_taxa(physeq_fig2_genus, Genus != "uncultured") #825 genera

genus_counts <- as(otu_table(physeq_fig2_genus), "matrix")
if (taxa_are_rows(physeq_fig2_genus)) { genus_counts <- t(genus_counts)}

groups <- c("Risor gut", "Risor skin", "sponge")

genus_sets <- psmelt(physeq_fig2_genus) %>%
  filter(
    type %in% groups,
    Abundance > 0
  ) %>%
  distinct(type, Genus) %>%
  group_by(type) %>%
  summarise(
    genera = list(Genus),
    .groups = "drop"
  ) %>%
  deframe()

#Plot Euler for genus level
euler_genus <- eulerr::euler(genus_sets)

#cairo_pdf(file = "euler_2B.pdf", width = 4, height = 3)
#plot(
#  euler_genus,
#  quantities = TRUE,
#  labels = TRUE,
#  fills = list(
#    fill = c(
#      "Risor gut" = "#44AA99",
#      "Risor skin" = "#DDCC77",
#      "sponge" = "#1D2D55"), alpha = 0.8))
#dev.off()

#### Figure 2 B - Shared proportion ####

# ASV level

#identify all ASVs detected in at least one sponge sample
sponge_samples <- rownames(metadata_fig2)[
  metadata_fig2$type == "sponge"]

sponge_asvs <- colnames(asv_counts)[
  colSums(
    asv_counts[sponge_samples, , drop = FALSE] > 0) > 0]

length(sponge_asvs)

#calculate the proportion separately for every fish sample
fish_samples <- rownames(metadata_fig2)[
  metadata_fig2$type %in% c("Risor gut", "Risor skin")]

asv_shared_reads <- data.frame(
  sample_id = fish_samples,
  type = metadata_fig2[fish_samples, "type"],
  total_reads = rowSums(
    asv_counts[fish_samples, , drop = FALSE]),
  shared_reads = rowSums(
    asv_counts[
      fish_samples,
      sponge_asvs,
      drop = FALSE])) %>%
  mutate(shared_percent = 100 * shared_reads / total_reads,
    type = factor(type, levels = c("Risor gut", "Risor skin")))

asv_shared_reads %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean_percent = mean(shared_percent),
    sd_percent = sd(shared_percent),
    median_percent = median(shared_percent),
    IQR_percent = IQR(shared_percent),
    min_percent = min(shared_percent),
    max_percent = max(shared_percent),
    .groups = "drop"
  )

#visualizing in the boxplot
plot_shared_asv <- ggplot(
  asv_shared_reads,
  aes(x = type,
    y = shared_percent,
    fill = type,
    color = type)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.6, width = 0.6) +
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.8) +
  scale_fill_manual(
    values = c("Risor gut" = "#44AA99","Risor skin" = "#DDCC77")) +
  scale_color_manual(values = c("Risor gut" = "#44AA99", "Risor skin" = "#DDCC77")) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.03))) +
  labs(
    x = NULL,
    y = paste(
      "Reads belonging to ASVs",
      "also detected in sponges (%)")) +
  theme_bw() +
  theme(legend.position = "none")

plot_shared_asv

ggsave("Shared_ASV_boxplot.pdf", plot_shared_asv, width = 4, height = 3, useDingbats = FALSE)

#Genus level
sponge_samples <- rownames(metadata_fig2)[
  metadata_fig2$type == "sponge"]

sponge_genus <- colnames(genus_counts)[
  colSums(
    genus_counts[sponge_samples, , drop = FALSE] > 0) > 0]

length(sponge_genus)

#calculate the proportion separately for every fish sample
fish_samples <- rownames(metadata_fig2)[
  metadata_fig2$type %in% c("Risor gut", "Risor skin")]

genus_shared_reads <- data.frame(
  sample_id = fish_samples,
  type = metadata_fig2[fish_samples, "type"],
  total_reads = rowSums(
    genus_counts[fish_samples, , drop = FALSE]),
  shared_reads = rowSums(
    genus_counts[
      fish_samples,
      sponge_genus,
      drop = FALSE])) %>%
  mutate(shared_percent = 100 * shared_reads / total_reads,
         type = factor(type, levels = c("Risor gut", "Risor skin")))

genus_shared_reads %>%
  group_by(type) %>%
  summarise(
    n = n(),
    mean_percent = mean(shared_percent),
    sd_percent = sd(shared_percent),
    median_percent = median(shared_percent),
    IQR_percent = IQR(shared_percent),
    min_percent = min(shared_percent),
    max_percent = max(shared_percent),
    .groups = "drop"
  )

#visualizing in the boxplot
plot_shared_genus <- ggplot(
  genus_shared_reads,
  aes(x = type,
      y = shared_percent,
      fill = type,
      color = type)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.6, width = 0.6) +
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.8) +
  scale_fill_manual(
    values = c("Risor gut" = "#44AA99","Risor skin" = "#DDCC77")) +
  scale_color_manual(values = c("Risor gut" = "#44AA99", "Risor skin" = "#DDCC77")) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.03))) +
  labs(
    x = NULL,
    y = paste(
      "Genus-assigned reads belonging to genera",
      "also detected in sponges (%)")) +
  theme_bw() +
  theme(legend.position = "none")

plot_shared_genus

#ggsave("Shared_Genus_boxplot.pdf", plot_shared_genus, width = 4, height = 3, useDingbats = FALSE)

#overall summary - Genus level

shared_summary <- bind_rows(
  asv_shared_reads %>%
    mutate(resolution = "ASV"),
  genus_shared_reads %>%
    mutate(resolution = "Genus")) %>%
  group_by(resolution, type) %>%
  summarise(
    n = n(),
    mean_percent = mean(shared_percent),
    sd_percent = sd(shared_percent),
    median_percent = median(shared_percent),
    IQR_percent = IQR(shared_percent),
    min_percent = min(shared_percent),
    max_percent = max(shared_percent),
    .groups = "drop")

shared_summary

#### Sensitivity analysis - Sharedness ####

# Threshold combinations
prevalence_cutoffs <- c(0.05, 0.10, 0.20)
read_cutoffs <- c(10, 50, 100)

genus_summary <- genus_counts %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(
    -sample_id,
    names_to = "Genus",
    values_to = "reads") %>%
  left_join(
    metadata_fig2 %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      dplyr::select(sample_id, type),
    by = "sample_id") %>%
  mutate(
    type = recode(
      type,
      "Risor gut" = "gut",
      "Risor skin" = "skin",
      "sponge" = "sponge")) %>%
  group_by(type, Genus) %>%
  summarise(
    pooled_reads = sum(reads),
    prevalence = mean(reads > 0),
    .groups = "drop") %>%
  pivot_wider(
    names_from = type,
    values_from = c(pooled_reads, prevalence),
    values_fill = 0)

#Function to calculate sensitivity results for one tissue
make_sensitivity_table <- function(data, tissue, tissue_label) {
  
  tissue_reads <- paste0("pooled_reads_", tissue)
  
  # No threshold, prevalence-only, and read-only criteria
  individual_thresholds <- tibble(
    `Filtering criterion` = c(
      "No additional threshold",
      "Sponge prevalence >=5%",
      "Sponge prevalence >=10%",
      "Sponge prevalence >=20%",
      ">10 reads in both groups",
      ">50 reads in both groups",
      ">100 reads in both groups"
    ),
    
    `Shared genera` = c(
      sum(
        data[[tissue_reads]] > 0 &
          data$pooled_reads_sponge > 0
      ),
      
      sum(
        data[[tissue_reads]] > 0 &
          data$prevalence_sponge >= 0.05
      ),
      
      sum(
        data[[tissue_reads]] > 0 &
          data$prevalence_sponge >= 0.10
      ),
      
      sum(
        data[[tissue_reads]] > 0 &
          data$prevalence_sponge >= 0.20
      ),
      
      sum(
        data[[tissue_reads]] > 10 &
          data$pooled_reads_sponge > 10
      ),
      
      sum(
        data[[tissue_reads]] > 50 &
          data$pooled_reads_sponge > 50
      ),
      
      sum(
        data[[tissue_reads]] > 100 &
          data$pooled_reads_sponge > 100
      )
    )
  )
  
  # All prevalence × read-count combinations
  combined_thresholds <- crossing(
    prevalence_cutoff = prevalence_cutoffs,
    read_cutoff = read_cutoffs
  ) %>%
    rowwise() %>%
    mutate(
      `Filtering criterion` = paste0(
        "Sponge prevalence >=",
        prevalence_cutoff * 100,
        "% and >",
        read_cutoff,
        " reads in both groups"
      ),
      
      `Shared genera` = sum(
        data$prevalence_sponge >= prevalence_cutoff &
          data[[tissue_reads]] > read_cutoff &
          data$pooled_reads_sponge > read_cutoff
      )
    ) %>%
    ungroup() %>%
    dplyr::select(
      `Filtering criterion`,
      `Shared genera`
    )
  
  bind_rows(
    individual_thresholds,
    combined_thresholds
  ) %>%
    mutate(
      Tissue = tissue_label,
      .before = 1
    )
}

#then apply for each sample type
gut_sensitivity <- make_sensitivity_table(
  data = genus_summary,
  tissue = "gut",
  tissue_label = "Risor gut"
)

skin_sensitivity <- make_sensitivity_table(
  data = genus_summary,
  tissue = "skin",
  tissue_label = "Risor skin"
)

sensitivity_table <- bind_rows(
  gut_sensitivity,
  skin_sensitivity
)

sensitivity_table
#write_tsv(sensitivity_table, "sensitivity_table_complete.tsv")

#gut sensitivity
gut_sensitivity <- genus_summary %>%
  summarise(
    `No additional threshold` =
      sum(pooled_reads_gut > 0 & pooled_reads_sponge > 0),
    `Sponge prevalence >=5%` =
      sum(pooled_reads_gut > 0 & prevalence_sponge >= 0.05),
    `Sponge prevalence >=10%` =
      sum(pooled_reads_gut > 0 & prevalence_sponge >= 0.10),
    `Sponge prevalence >=20%` =
      sum(pooled_reads_gut > 0 & prevalence_sponge >= 0.20),
    `>10 reads in both groups` =
      sum(pooled_reads_gut > 10 & pooled_reads_sponge > 10),
    `>50 reads in both groups` =
      sum(pooled_reads_gut > 50 & pooled_reads_sponge > 50),
    `>100 reads in both groups` =
      sum(pooled_reads_gut > 100 & pooled_reads_sponge > 100),
    `>=10% prevalence and >50 reads` =
      sum(
        prevalence_sponge >= 0.10 &
          pooled_reads_gut > 50 &
          pooled_reads_sponge > 50)) %>%
  pivot_longer(
    everything(),
    names_to = "Filtering criterion",
    values_to = "Shared genera") %>%
  mutate(Tissue = "Risor gut")

#skin sensitivity
skin_sensitivity <- genus_summary %>%
  summarise(
    `No additional threshold` =
      sum(pooled_reads_skin > 0 & pooled_reads_sponge > 0),
    `Sponge prevalence >=5%` =
      sum(pooled_reads_skin > 0 & prevalence_sponge >= 0.05),
    `Sponge prevalence >=10%` =
      sum(pooled_reads_skin > 0 & prevalence_sponge >= 0.10),
    `Sponge prevalence >=20%` =
      sum(pooled_reads_skin > 0 & prevalence_sponge >= 0.20),
    `>10 reads in both groups` =
      sum(pooled_reads_skin > 10 & pooled_reads_sponge > 10),
    `>50 reads in both groups` =
      sum(pooled_reads_skin > 50 & pooled_reads_sponge > 50),
    `>100 reads in both groups` =
      sum(pooled_reads_skin > 100 & pooled_reads_sponge > 100),
    `>=10% prevalence and >50 reads` =
      sum(
        prevalence_sponge >= 0.10 &
          pooled_reads_skin > 50 &
          pooled_reads_sponge > 50)) %>%
  pivot_longer(
    everything(),
    names_to = "Filtering criterion",
    values_to = "Shared genera") %>%
  mutate(Tissue = "Risor skin")

#combine
sensitivity_table <- bind_rows(
  gut_sensitivity,
  skin_sensitivity) %>%
  dplyr::select(Tissue, `Filtering criterion`, `Shared genera`)

sensitivity_table

write_tsv(sensitivity_table, "sensitivity_analysis.tsv")

#### Figure 3 Overall sharedness ####

physeq.fig3.g<- physeq.asv.bacteria #13780 ASVs
physeq.fig3.g=tax_glom(physeq.fig3.g, taxrank="Genus")

physeq.fig3.g <- subset_taxa(physeq.fig3.g, Genus != "NA") #removing unwanted taxa
physeq.fig3.g <- subset_taxa(physeq.fig3.g, Genus != "uncultured") 
physeq.fig3.g <- subset_taxa(physeq.fig3.g, Phylum != "Cyanobacteria") #914

#check thresholds: 5% of 43 sponges represent 3 sponges
n_sponges <- sum(sample_data(physeq.fig3.g)$type == "sponge")
n_sponges
minimum_sponges<- ceiling(0.05 * n_sponges)

sponge_prevalence <- psmelt(physeq.fig3.g) %>%
  filter(
    type == "sponge",
    Abundance > 0
  ) %>%
  distinct(Sample, Genus) %>%
  count(Genus, name = "n_sponges_present") %>%
  mutate(
    prevalence = n_sponges_present / n_sponges
  )


genera_keep <- sponge_prevalence %>%
  filter(n_sponges_present >= minimum_sponges) %>%
  pull(Genus)

length(genera_keep)

#retain genera present in at least 5% of sponges
genera_keep <- sponge_prevalence %>%
  filter(prevalence >= 0.05) %>%
  pull(Genus)

physeq.fig3.g.5 <- subset_taxa(physeq.fig3.g, Genus %in% genera_keep)

physeq.fig3.g.5 <- prune_taxa(taxa_sums(physeq.fig3.g.5) > 0, physeq.fig3.g.5) #98 genera

sort(unique(sample_data(physeq.fig3.g.5)$type))

fig3_groups <- c(
  "sponge",
  "Risor gut",
  "Risor skin",
  "E.prochilos gut",
  "E.prochilos skin",
  "seawater"
)

#create genus sets
genus_sets_fig3 <- psmelt(physeq.fig3.g.5) %>%
  filter(
    type %in% fig3_groups,
    Abundance > 0
  ) %>%
  distinct(type, Genus) %>%
  group_by(type) %>%
  summarise(
    genera = list(Genus),
    .groups = "drop"
  ) %>%
  deframe()

genus_sets_fig3 <- genus_sets_fig3[fig3_groups]

names(genus_sets_fig3) <- c(
  "Sponge",
  "R. ruber gut",
  "R. ruber skin",
  "E. prochilos gut",
  "E. prochilos skin",
  "Seawater"
)

#checking the number of genera detected in each group
sapply(genus_sets_fig3, length)

#convert to upset format
upset_data_fig3 <- fromList(genus_sets_fig3)


upset_data_shared <- upset_data_fig3[
  upset_data_fig3$Sponge == 1 &
    rowSums(upset_data_fig3[, set_order]) >= 2,
  ,
  drop = FALSE
]

#plotting and saving upset plot
cairo_pdf(file = "upset_plot_fig3.pdf", width = 9, height = 4.5)

upset(
  upset_data_shared,
  
  sets = c(
    "E. prochilos gut",
    "E. prochilos skin",
    "Seawater",
    "R. ruber gut",
    "R. ruber skin",
    "Sponge"
  ),
  
  keep.order = TRUE,
  nsets = 6,
  nintersects = NA,
  order.by = "freq",
  
  mainbar.y.label = "Number of genera",
  sets.x.label = "Genera detected",
  
  main.bar.color = "#1D2D55",
  sets.bar.color = "grey40",
 # matrix.color = "#1D2D55",
  
  point.size = 3.5,
  line.size = 1,
  mb.ratio = c(0.65, 0.35),
  
  text.scale = c(
    1.3,  # intersection-size title
    1.1,  # intersection-size numbers
    1.0,  # set-size title
    1.0,  # set-size numbers
    1.1,  # set names
    1.0   # matrix labels
  )
)

dev.off()

#### Figure 4 Key Shared Genera ####

fig4_long<- psmelt(physeq.fig3.g)

#calculate sponge prevalence - now staying with the 10% threhsold
n_sponges <- fig4_long %>%
  filter(type == "sponge") %>%
  distinct(Sample) %>%
  nrow()

min_sponges_10 <- ceiling(0.10 * n_sponges)

sponge_prevalence <- fig4_long %>%
  filter(
    type == "sponge",
    Abundance > 0
  ) %>%
  distinct(Sample, OTU, Genus, Phylum, Family) %>%
  count(
    OTU,
    Genus,
    Phylum,
    Family,
    name = "n_sponges_present"
  ) %>%
  mutate(
    sponge_prevalence = n_sponges_present / n_sponges)

#calculate pooled reads by sample type
genus_counts <- fig4_long %>%
  group_by(
    OTU,
    Genus,
    Phylum,
    Family,
    type
  ) %>%
  summarise(
    reads = sum(Abundance),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = reads,
    values_fill = 0
  ) %>%
  left_join(
    sponge_prevalence,
    by = c("OTU", "Genus", "Phylum", "Family")
  )

#key gut and sponge
key_gut <- genus_counts %>%
  filter(
    n_sponges_present >= min_sponges_10,
    `Risor gut` > 50,
    sponge > 50
  )

#key skin and sponge
key_skin <- genus_counts %>%
  filter(
    n_sponges_present >= min_sponges_10,
    `Risor skin` > 50,
    sponge > 50
  )

key_gut %>%
  dplyr::select(Genus, n_sponges_present, `Risor gut`, sponge)

key_skin %>%
  dplyr::select(Genus, n_sponges_present, `Risor skin`, sponge)

#restricting to only Risor and sponge interactions
key_gut_filtered <- key_gut %>%
  filter(
    seawater == 0,
    `E.prochilos gut` == 0,
    `E.prochilos skin` == 0,
    sponge > `Risor gut`)%>%
  mutate(
    Risor_gut_reads = `Risor gut`,
    Sponge_reads = sponge,
    Risor_gut_RA = Risor_gut_reads / sum(Risor_gut_reads),
    Sponge_RA = Sponge_reads / sum(Sponge_reads)
  )%>%
  filter(!(Genus == "EF100-94H03"))

key_skin_filtered <- key_skin %>%
  filter(
    seawater == 0,
    `E.prochilos gut` == 0,
    `E.prochilos skin` == 0,
    sponge > `Risor skin`)%>%
  mutate(
    Risor_skin_reads = `Risor skin`,
    Sponge_reads = sponge,
    Risor_skin_RA = Risor_skin_reads / sum(Risor_skin_reads),
    Sponge_RA = Sponge_reads / sum(Sponge_reads)
  )%>%
  filter(!(Genus == "Microbulbifer"))%>%
  filter(!(Genus == "Fulvivirga"))

key_shared_all <- bind_rows(
  key_gut_filtered,
  key_skin_filtered)

write_tsv(key_shared_all, "key_shared_all.tsv")

#calculate relative abundance

plot_gut <- key_gut_filtered %>%
  mutate(
    Risor_gut_reads = `Risor gut`,
    Sponge_reads = sponge,
    Risor_gut_RA = Risor_gut_reads / sum(Risor_gut_reads),
    Sponge_RA = Sponge_reads / sum(Sponge_reads)
  )%>%
  filter(!(Genus == "EF100-94H03"))

plot_skin <- key_skin_filtered %>%
  mutate(
    Risor_skin_reads = `Risor skin`,
    Sponge_reads = sponge,
    Risor_skin_RA = Risor_skin_reads / sum(Risor_skin_reads),
    Sponge_RA = Sponge_reads / sum(Sponge_reads)
  )%>%
  filter(!(Genus == "Microbulbifer"))%>%
  filter(!(Genus == "Fulvivirga"))




selected_reads <- key_shared_all %>%
  summarise(
    `Risor gut` = sum(`Risor gut`, na.rm = TRUE),
    `Risor skin` = sum(`Risor skin`, na.rm = TRUE),
    sponge = sum(sponge, na.rm = TRUE)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "type",
    values_to = "selected_reads"
  )

total_reads <- fig4_long %>%
  filter(type %in% c("Risor gut", "Risor skin", "sponge")) %>%
  group_by(type) %>%
  summarise(
    total_reads = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

read_representation <- selected_reads %>%
  left_join(total_reads, by = "type") %>%
  mutate(
    percent_reads = 100 * selected_reads / total_reads
  )

read_representation

read_representation %>%
  mutate(percent_reads = round(percent_reads, 2))

#Adjust data for plotting

#Skin
long_data_mirrored.skin <- plot_skin %>%
  mutate(
    Family = case_when(
      Family == "PAUC43f_marine_benthic_group" ~ "PAUC43f",
      Family == "BD2-11_terrestrial_group" ~ "BD2-11",
      TRUE ~ Family
    ),
    
    Genus = case_when(
      Genus == "PAUC43f_marine_benthic_group" ~ "PAUC43f_marine",
      Genus == "BD2-11_terrestrial_group" ~ "BD2-11_terrestrial",
      TRUE ~ Genus
    ),
    
    Phylum_Genus = paste(Phylum, Family, Genus, sep = ":"),
    total_ab = Risor_skin_reads + Sponge_reads
  ) %>%
  arrange(total_ab) %>%
  mutate(
    Phylum_Genus = factor(
      Phylum_Genus,
      levels = unique(Phylum_Genus)
    )
  ) %>%
  dplyr::select(
    Phylum_Genus,
    Risor_skin_RA,
    Sponge_RA
  ) %>%
  pivot_longer(
    cols = c(Risor_skin_RA, Sponge_RA),
    names_to = "Sample_Type",
    values_to = "Occurrences"
  ) %>%
  mutate(
    Sample_Type = recode(
      Sample_Type,
      Risor_skin_RA = "Risor skin",
      Sponge_RA = "Sponge"
    ),
    
    Occurrences = if_else(
      Sample_Type == "Risor skin",
      -Occurrences,
      Occurrences
    )
  )

# Plotting the mirrored barplot - Figure 4b
shared_barplot.skin <- ggplot(long_data_mirrored.skin, aes(x = Occurrences, y = Phylum_Genus, fill = Sample_Type)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = abs, expand = c(0, 0)) +
  scale_fill_manual(values = c("#DDCC77","#1D2D55")) + 
  labs(
    x = "Relative abundance",
    fill = "Sample Type"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("Shared.genus.skin.sponge.filtered.pdf", shared_barplot.skin , width = 8, height = 6, useDingbats = FALSE)

#Adjust data for plotting - gut

long_data_mirrored.gut <- plot_gut %>%
  mutate(
    Family = case_when(
      Family == "PAUC43f_marine_benthic_group" ~ "PAUC43f",
      Family == "BD2-11_terrestrial_group" ~ "BD2-11",
      TRUE ~ Family
    ),
    
    Genus = case_when(
      Genus == "PAUC43f_marine_benthic_group" ~ "PAUC43f_marine",
      Genus == "BD2-11_terrestrial_group" ~ "BD2-11_terrestrial",
      TRUE ~ Genus
    ),
    
    Phylum_Genus = paste(Phylum, Family, Genus, sep = ":"),
    total_ab = Risor_gut_reads + Sponge_reads
  ) %>%
  arrange(total_ab) %>%
  mutate(
    Phylum_Genus = factor(
      Phylum_Genus,
      levels = unique(Phylum_Genus)
    )
  ) %>%
  dplyr::select(
    Phylum_Genus,
    Risor_gut_RA,
    Sponge_RA
  ) %>%
  pivot_longer(
    cols = c(Risor_gut_RA, Sponge_RA),
    names_to = "Sample_Type",
    values_to = "Occurrences"
  ) %>%
  mutate(
    Sample_Type = recode(
      Sample_Type,
      Risor_gut_RA = "Risor gut",
      Sponge_RA = "Sponge"
    ),
    
    Occurrences = if_else(
      Sample_Type == "Risor gut",
      -Occurrences,
      Occurrences
    )
  )

# Plotting the mirrored barplot 
shared_barplot.gut <- ggplot(long_data_mirrored.gut, aes(x = Occurrences, y = Phylum_Genus, fill = Sample_Type)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = abs, expand = c(0, 0)) +
  scale_fill_manual(values = c("#44AA99","#1D2D55")) + 
  labs(
    x = "Relative abundance",
    fill = "Sample Type"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave("Shared.genus.gut.sponge.filtered.pdf", shared_barplot.gut , width = 8, height = 6, useDingbats = FALSE)




