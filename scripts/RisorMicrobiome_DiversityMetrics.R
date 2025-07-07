######################################################################
################# Risor microbiome project ###########################
########### DATA Analyses - Diversity Metrics (Figure 2) #############
#####################################################################

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

#### 1- Environment set up and data ####
setwd("User/workingdirectory") #your working directory

#ASV table - filtered
asv.all<- read.csv("asv.all.filtered.csv")

asv.all<- asv.all%>%
  column_to_rownames("ASV")
asv_m.all<- as.matrix(asv.all)

#Taxon table - filtered
tax.all<- read.csv("taxon.all.filtered.csv")

tax.all<- tax.all%>%
  column_to_rownames("ASV")
tax_m.all<- as.matrix(tax.all)

#Metadata table - filtered
meta.all<- read.csv("metadata_filtered_all_2000.csv")

meta.all<- meta.all%>%
  column_to_rownames("sample_id")
meta.all<-sample_data(meta.all)

#Generating phyloseq object with ASV and Taxonomy 
ASV.all<- otu_table(asv_m.all, taxa_are_rows = TRUE)
TAX.all<- tax_table(tax_m.all)

physeq.all<- phyloseq(ASV.all, TAX.all)

#Importing tree from qiime to add to phyloseq object
qiime_tree <- read.tree("tree.nwk")

head(qiime_tree$tip.label) #check tip labels

qiime_tree_filtered <- prune_taxa(taxa_names(physeq.all), phy_tree(phyloseq(qiime_tree))) #Filter the tree to only keep taxa present in your phyloseq object

#ASV LEVEL - Phyloseq object
physeq.all<- merge_phyloseq(physeq.all, meta.all, qiime_tree_filtered) 
physeq.all #all the data #14499 ASVs and 165 samples (all)

tree_clean <- ape::keep.tip(phy_tree(physeq1.all), taxa_names(physeq1.all)) #Get tree and prune it to match physeq

physeq.all <- subset_taxa(physeq.all, Genus != "Chloroplast") #removing unwanted taxa
physeq.all <- subset_taxa(physeq.all, Genus != "Mitochondria") 
physeq.all <- subset_taxa(physeq.all, Phylum != "Cyanobacteria") 

tree_clean <- ape::keep.tip(phy_tree(physeq.all), taxa_names(physeq.all)) #obtaining cleaned tree

physeq1.all<- merge_phyloseq(physeq.all, tree_clean) #generating final phyloseq object

#GENUS LEVEL - Phyloseq object
physeq1.all.g <- subset_taxa(physeq1.all, Genus != "NA") #removing unwanted taxa
physeq1.all.g <- subset_taxa(physeq1.all.g, Genus != "uncultured") 
physeq1.all.g <- subset_taxa(physeq1.all.g, Genus != "Chloroplast") 
physeq1.all.g <- subset_taxa(physeq1.all.g, Genus != "Mitochondria") 
physeq1.all.g <- subset_taxa(physeq1.all.g, Phylum != "Cyanobacteria") 

physeq1.all.g=tax_glom(physeq1.all.g, taxrank="Genus") #914 genus

gen_df <- as.data.frame(otu_table(physeq1.all.g))

#Filtering out genera with less than 50 reads 
taxa_to_keep_all <- gen_df %>%
  rownames_to_column(var = "Taxa") %>%  # Keep track of taxa names
  pivot_longer(-Taxa, names_to = "Sample", values_to = "Abundance") %>% 
  filter(Abundance > 0) %>%  # Consider only present taxa
  group_by(Taxa) %>% 
  summarise(MaxAbundance = max(Abundance)) %>%  # Find max abundance per taxon
  filter(MaxAbundance > 50) %>%  # Keep only those where max is > 50
  pull(Taxa)  # Extract taxa names as a vector

physeq1.all.g.f <- prune_taxa(taxa_to_keep_all, physeq1.all.g) #594 genera
gen_df_f <- as.data.frame(otu_table(physeq1.all.g.f))
#write.csv(gen_df_f, "genus.filtered.50ab.original.csv")


#### 2 - General Abundance metrics ####

#Obtain most abundant phylum for each sample type
physeq_rskin <- subset_samples(physeq1.all, type == "Risor skin")
#physeq_rgut <- subset_samples(physeq1.all, type == "Risor gut")
#physeq_sponge <- subset_samples(physeq1.all, type == "sponge")
#physeq_eskin <- subset_samples(physeq1.all, type == "E.prochilos skin")
#physeq_egut <- subset_samples(physeq1.all, type == "E.prochilos gut")
#physeq_seawater <- subset_samples(physeq1.all, type == "seawater")

#Extract the taxonomic table and the sample counts - for each of sample types
tax_table_rskin <- tax_table(physeq_rskin)
otu_table_rskin <- otu_table(physeq_rskin)

tax_df_rskin <- as.data.frame(tax_table_rskin)

#Summarize ASV counts by phylum
phyla_abundance_rskin <- otu_table_rskin %>%
  as.data.frame() %>%
  rownames_to_column("ASV") %>%
  gather(key = "Sample", value = "Abundance", -ASV) %>%
  left_join(tax_df %>% rownames_to_column("ASV"), by = "ASV") %>%
  filter(!is.na(Phylum)) %>% 
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup() %>%
  arrange(desc(Abundance))

#Calculate relative abundance
phyla_abundance_rskin <- phyla_abundance_rskin %>%
  group_by(Sample) %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>%
  ungroup()

#Get the top 10 phyla by relative abundance 
top_10_phyla <- phyla_abundance %>%
  group_by(Phylum) %>%
  summarise(TotalRelativeAbundance = sum(RelativeAbundance)) %>%
  arrange(desc(TotalRelativeAbundance)) %>%
  head(10)
#write.csv(top_10_phyla2, "top_10_phyla.csv")


#### 3 - Basic stats - Alpha diversity  ####

##Alpha diversity - ASV level##

#Calculate alpha diversity metrics for Shannon and Observed
alpha_metrics <- estimate_richness(physeq1.all, measures = c("Shannon", "Observed"))

#Combine alpha metrics with sample metadata
alpha_data <- cbind(alpha_metrics, as(sample_data(physeq1.all), "data.frame"))

alpha_data2<- alpha_data %>%
  tibble::rownames_to_column(var = "SampleID")

#Phylogenetic diversity - ASV level
otu <- as(otu_table(physeq1.all), "matrix") #isolating asv table
if(taxa_are_rows(physeq1.all)) { otu <- t(otu) }

tree <- phy_tree(physeq1.all) #isolating tree

pd_result <- pd(otu, tree, include.root = TRUE)

pd_result<- pd_result %>%
  tibble::rownames_to_column(var = "SampleID")

alpha_merged<- left_join(alpha_data2, pd_result, by = "SampleID") #merging PD with the other alpha diversity metrics

##Testing the effect of SAMPLE TYPE on ALPHA DIVERSITY using Linear Models##
#The same functions were used to plot alpha diversity with raw ASVs - supp material

#Shannon
hist(alpha_merged$Shannon) #checking data distribution

lm_shannon<- lm(Shannon ~ type, data = alpha_merged) #trying linear model first
hist(resid(lm_shannon))
qqnorm(resid(lm_shannon)); qqline(resid(lm_shannon))

glmShannon <- glm(Shannon ~ sponge_spp, family = Gamma(link = "log"), data = alpha_sponge) #trying glm with gamma distribution

AIC(lm_shannon, glmShannon) #smaller AIC for lm, so sticking with it
summary(lm_shannon)

#Observed
hist(alpha_sponge$Observed) #checking data distribution

glm_nb <- glm.nb(Observed ~ type, data = alpha_merged) #glm with negative binomial distribution
summary(glm_nb)

#PD
hist(alpha_sponge$PD)

glm_pd <- glm(PD ~ sponge_spp, data = alpha_sponge, family = Gamma(link = "log")) #glm with gamma distribution
summary(glm_pd)


##Plotting alpha diversity - sample type - Figure 2A##
alpha_long <- pivot_longer(alpha_merged, 
                           cols = c("Observed", "Shannon", "PD"), 
                           names_to = "metric", 
                           values_to = "value")

alpha_plot <- ggplot(alpha_long, aes(x = type, y = value, fill = type, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + 
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) + 
  facet_wrap(~ metric, scales = "free_y") +           
  theme_bw() +
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  labs(x = "Sample Type", y = "Alpha Diversity") +
  scale_fill_manual(values = c("Risor skin" = "#DDCC77", 
                               "sponge" = "#661100", 
                               "Risor gut" = "#44AA99",
                               "E.prochilos skin" = "#B8887A", 
                               "E.prochilos gut" = "darkslategrey", 
                               "seawater" = "darkblue")) +
  scale_color_manual(values = c("Risor skin" = "#DDCC77", 
                                "sponge" = "#661100", 
                                "Risor gut" = "#44AA99", 
                                "E.prochilos skin" = "#B8887A", 
                                "E.prochilos gut" = "darkslategrey", 
                                "seawater" = "darkblue"))

#ggsave("alpha_boxplot_genus_all_metrics_ASV_new.pdf", alpha_plot , width = 10, height = 4, useDingbats = FALSE)

##Testing the effect of SPONGE SPECIES on ALPHA DIVERSITY - Supp figure##

alpha_sponge<- alpha_merged%>%
  filter(type == "sponge")

#Shannon
hist(alpha_sponge$Shannon)

lm_shannon<- lm(Shannon ~ sponge_spp, data = alpha_sponge)
hist(resid(lm_shannon))
qqnorm(resid(lm_shannon)); qqline(resid(lm_shannon))

glmShannon <- glm(Shannon ~ sponge_spp, family = Gamma(link = "log"), data = alpha_sponge)
AIC(lm_shannon, glmShannon) #smaller AIC for lm, so sticking with it
summary(lm_shannon)

#Observed
hist(alpha_sponge$Observed)

glm_nb <- glm.nb(Observed ~ sponge_spp, data = alpha_sponge)
summary(glm_nb)

#PD
hist(alpha_sponge$PD)
glm_pd <- glm(PD ~ sponge_spp, data = alpha_sponge, family = Gamma(link = "log"))
summary(glm_pd)

##Plotting alpha diversity - sponge species##
custom_colors <- c(
  "I. campana" = "lightsalmon3",
  "I. felix" = "saddlebrown",
  "I. ruetzleri" = "peachpuff3",
  "I. strobilina" = "gray38"
)

alpha_long_sp <- alpha_long %>%
  filter(!is.na(sponge_spp))%>%
  mutate(sponge_spp = fct_recode(sponge_spp,
                                 "I. campana"    = "Ircinia campana",
                                 "I. felix"      = "Ircinia felix",
                                 "I. ruetzleri"  = "Ircinia ruetzleri",
                                 "I. strobilina" = "Ircinia strobilina"
  ))


alpha_plot_sp <- ggplot(alpha_long_sp, aes(x = sponge_spp, y = value, fill = sponge_spp, color = sponge_spp)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +   
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) + 
  facet_wrap(~ metric, scales = "free_y") +          
  theme(legend.position = "none", 
        axis.title.x = element_blank()) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Sample Type", y = "Alpha Diversity") 
alpha_plot_sp

#ggsave("alpha_boxplot_asv_sponge_species.pdf", alpha_plot_sp , width = 10, height = 4, useDingbats = FALSE)


#### 4 - Basic stats - Beta diversity  ####

##ASV level and sample type - supp material

physeq1.all@sam_data$type <- factor(physeq1.all@sam_data$type, 
                                    levels = c("Risor gut", "Risor skin", "sponge", 
                                               "E.prochilos gut", "E.prochilos skin", "seawater"))

pcoa.asv.p.w<- ordinate(physeq1.all, "PCoA", "unifrac", weighted=FALSE) #Unifrac distance
pcoa.plot.asv.p.w<- plot_ordination(physeq1.all, pcoa.asv.p.w, color="type")
pcoa.plot.asv.p.w<- pcoa.plot.asv.p.w + geom_point(size=2)+
  theme_bw()+
  scale_colour_manual(values = c("Risor skin" = "#DDCC77", 
                                 "sponge" = "#661100", 
                                 "Risor gut" = "#44AA99",
                                 "E.prochilos skin" = "#B8887A", 
                                 "E.prochilos gut" = "darkslategrey", 
                                 "seawater" = "darkblue"))+
  labs(color="Sample Type")+
  stat_ellipse(aes(group = type), level = 0.95, linetype = 2)
pcoa.plot.asv.p.w

#ggsave("pcoa_sampletype_all_ASV_new.pdf", pcoa.plot.asv.p.w , width = 8, height = 6, useDingbats = FALSE)

alpha_pca<- alpha_plot / pcoa.plot.asv.p.w

##Genus level and sample type - Figure 2B
physeq1.all.g@sam_data$type <- factor(physeq1.all.g@sam_data$type, 
                                      levels = c("Risor gut", "Risor skin", "sponge", 
                                                 "E.prochilos gut", "E.prochilos skin", "seawater"))

pcoa.gen<- ordinate(physeq1.all.g, "PCoA", "unifrac", weighted=FALSE)
pcoa.gen<- plot_ordination(physeq1.all.g, pcoa.gen, color="type")
pcoa.gen.plot<- pcoa.gen + geom_point(size=2)+
  theme_bw()+
  scale_colour_manual(values = c("Risor skin" = "#DDCC77", 
                                 "sponge" = "#661100", 
                                 "Risor gut" = "#44AA99",
                                 "E.prochilos skin" = "#B8887A", 
                                 "E.prochilos gut" = "darkslategrey", 
                                 "seawater" = "darkblue"))+
  labs(color="Sample Type")+
  stat_ellipse(aes(group = type), level = 0.95, linetype = 2)
pcoa.gen.plot
ggsave("pcoa_sampletype_all_Genus_new.pdf", pcoa.gen.plot , width = 8, height = 6, useDingbats = FALSE)

#Checking whether unifrac differences are statistically different
unifrac_dist <- phyloseq::distance(physeq1.all.g, method = "unifrac", weighted = FALSE)

sample_meta <- as(sample_data(physeq1.all.g), "data.frame")

#PERMANOVA using adonis
set.seed(123)  # For reproducibility of the permutation test
adonis_res <- adonis2(unifrac_dist ~ type, data = sample_meta, permutations = 999)

print(adonis_res)

#Checking whether differences are due to within group dispersion
disp <- betadisper(unifrac_dist, sample_meta$type)
anova(disp)
permutest(disp, pairwise = TRUE, permutations = 999) #pairwise
plot(disp)
boxplot(disp) #assumption of homogeneous dispersion is violated by all variables except seawater and sponge


##Genus level and Sponge species - supp material

custom_colors <- c(
  "Ircinia campana" = "lightsalmon3",
  "Ircinia felix" = "saddlebrown",
  "Ircinia ruetzleri" = "peachpuff3",
  "Ircinia strobilina" = "gray38"
)

# Remove samples where sponge_spp is NA
physeq1.g.filtered <- subset_samples(physeq1.all.g, !is.na(sponge_spp))

physeq1.g.filtered@sam_data$sponge_spp <- factor(physeq1.g.filtered@sam_data$sponge_spp, 
                                                 levels = c("Ircinia campana", "Ircinia felix", "Ircinia ruetzleri", 
                                                            "Ircinia strobilina"))

pcoa.gen.sponge<- ordinate(physeq1.g.filtered, "PCoA", "unifrac", weighted=FALSE)

pcoa.gen.sponge<- plot_ordination(physeq1.g.filtered, pcoa.gen.sponge, color="sponge_spp")

pcoa.gen.plot.sponge<- pcoa.gen.sponge + geom_point(size=2)+
  theme_bw()+
  scale_colour_manual(values = custom_colors)+
  labs(color="Sponge species")+
  stat_ellipse(aes(group = sponge_spp), level = 0.95, linetype = 2)
pcoa.gen.plot.sponge

#ggsave("pcoa_spongespp_Genus_new.pdf", pcoa.gen.plot.sponge, width = 8, height = 6, useDingbats = FALSE)


#Checking whether unifrac differences are statistically different
unifrac_dist <- phyloseq::distance(physeq1.g.filtered, method = "unifrac", weighted = FALSE)

sample_meta <- as(sample_data(physeq1.g.filtered), "data.frame")

#PERMANOVA using adonis
set.seed(123)  # For reproducibility of the permutation test
adonis_res_sponge <- adonis2(unifrac_dist ~ sponge_spp, data = sample_meta, permutations = 999)

print(adonis_res_sponge)