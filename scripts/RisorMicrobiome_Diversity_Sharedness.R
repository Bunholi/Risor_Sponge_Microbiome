
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

######################################################################
################# Risor microbiome project ###########################
########### DATA Analyses - Diversity Metrics (Figure 1) #############
#####################################################################

#### Environment set up and data wrangling ####
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
#saveRDS(ps, file = "phyloseq_genus_50ab_all.rds")
gen_df_f <- as.data.frame(otu_table(physeq1.all.g.f))
#write.csv(gen_df_f, "genus.filtered.50ab.original.csv")


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


#### Alpha diversity - Figure 1C ####

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


##Plotting alpha diversity - sample type - Figure 1C##
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


#### Beta diversity - Figure 1D ####

##ASV level and sample type - supp material Figure S7

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

##Genus level and sample type - Figure 1D
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


##Genus level and Sponge species - supp material Figure S5

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



##############################################################################
######################### Risor microbiome project ###########################
############## DATA Analyses - Shared Microbiome (Figure 2 and 4) ############
##############################################################################

#### Shared microbiome sample types (UpsetR plot) - Figure 2A ####
tax_table_gen.all #594 genus - Genus with more than 50 reads overall
tax_table_gen.all <- tax_table_gen.all%>%
  dplyr::select(-c("Species", "pident"))

# Define the categorical columns to exclude from transformation
categorical_columns <- c("Genus", "Phylum", "Class", "Family", "Order", "ASV")

# Convert ASV counts to presence-absence (1 for presence, 0 for absence)
df_pa <- tax_table_gen.all %>%
  mutate(across(-all_of(categorical_columns), ~ ifelse(. > 0, 1, 0)))%>%
  dplyr::select(-c(ASV, Phylum, Class, Order, Family, Genus))

df_numeric <- df_pa %>%
  dplyr::select(where(is.numeric))

df_filtered <- df_numeric %>%
  rowwise() %>%
  mutate(total_presence = sum(c_across(everything()), na.rm = TRUE)) %>%  # Count in how many sets each taxon appears
  filter(total_presence > 1)%>%  # Remove taxa that appear in only one set
  dplyr::select(-total_presence)  # Drop helper column

df_filtered <- as.data.frame(df_filtered)


# Create the UpSet plot with custom styling
pdf("upset_plot_all_v4.pdf", width = 9, height = 4.5, useDingbats = FALSE)

upset(
  df_filtered, 
  sets = c("seawater", "E.prochilos skin", "E.prochilos gut", "sponge", "Risor skin", "Risor gut"), 
  sets.bar.color = "#CC6677",   # Custom bar color
  order.by = "freq",           # Order intersections by frequency
  empty.intersections = "off",  # Show empty sets
  mainbar.y.label = "Intersection Size",  # Custom Y-axis label for main bar plot
  sets.x.label = "Genus Count",           # Custom X-axis label for set size plot
  text.scale = c(1.5, 1.5, 1.2, 1.2, 1.2, 1.0),  # Adjust text sizes (titles, axis, labels)
  matrix.color = "#888888",  # Change the matrix dots to a custom color
  point.size = 3.5,  # Increase the matrix point size
  line.size = 1,  # Make intersection lines thicker
  keep.order = TRUE,  # Keep the order of sets as defined
)

dev.off()

#### Proportion shared (donut plot) - Figure 2B and 2C ####

#Donut plot - SKIN

#Selecting sponge IDs
skin <- meta.all %>%
  filter(grepl("^SBEL", sample_id))%>%
  pull(sample_id)

#Subset physeq object specific for sponges
physeq_skin<- prune_samples(skin, physeq1.all.g.f)

gen_skin<- otu_table(physeq_skin)
gen_skin[gen_skin > 0] <- 1 #convert counts to presence/absence

#Calculate the proportion of sponge samples each taxon is found in
gen_prev.skin<- rowSums(gen_skin) / length(skin)

gen_to_keep.01.skin<- taxa_names(physeq1.all.g)[gen_prev.skin >= 0.1]

physeq1.all.g.01.skin<- prune_taxa(gen_to_keep.01.skin, physeq1.all.g)

physeq_skin_type = merge_samples(physeq1.all.g.01.skin, "type")
sample_data(physeq_skin_type)
table_physeq.s <- t(otu_table(physeq_skin_type, ))

table_physeq.s<- as.data.frame(table_physeq.s)
table_physeq.s<- rownames_to_column(table_physeq.s)
table_physeq.s<- table_physeq.s %>%
  dplyr::rename(ASV = rowname)

physeq_skin_tax<- tax_table(physeq1.all.g.01.skin)
physeq_skin_tax<- as.data.frame(physeq_skin_tax)
physeq_skin_tax<- rownames_to_column(physeq_skin_tax)
physeq_skin_tax<- physeq_skin_tax %>%
  dplyr::rename(ASV = rowname)

tax_table_gen.s<- physeq_skin_tax %>%
  inner_join(table_physeq.s, by = "ASV")%>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, pident, `Risor skin`, `Risor gut`, `sponge`)

filtered_tax_table_gen.skin.only<- tax_table_gen.s %>%
  filter(`Risor skin` > 50) #209 genus
write.csv(filtered_tax_table_gen.skin.only, "skin_only_genus.csv")

tax_table_gen.all.f.skin #shared between risor skin and sponge >50 reads - 17

count_skin.only<- nrow(filtered_tax_table_gen.skin.only)
count_shared.s.s<- nrow(tax_table_gen.all.f.skin)

# Create a data frame for visualization
donut_data.skin <- data.frame(
  Category = c("Risor skin & Sponge", "Risor skin"),
  count = c(count_shared.s.s, count_skin.only - count_shared.s.s)
)

# Calculate percentages and cumulative values
donut_data.skin$fraction <- donut_data.skin$count / sum(donut_data.skin$count)
donut_data.skin$ymax <- cumsum(donut_data.skin$fraction)
donut_data.skin$ymin <- c(0, head(donut_data.skin$ymax, n=-1))
donut_data.skin$label_position <- (donut_data.skin$ymax + donut_data.skin$ymin) / 2  # Midpoint for text labels
donut_data.skin$label_text <- paste0(round(donut_data.skin$fraction * 100, 1), "%")  # Percentage text

# Create a donut plot - figure 2B
donut_plot.skin <- ggplot(donut_data.skin, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
  geom_rect() +
  geom_text(aes(x = 3.5, y = label_position, label = label_text), size = 5, color = "white") +  # Add percentage labels
  coord_polar(theta = "y") +  # Convert to donut
  xlim(c(2, 5)) +  # Adjust for donut shape
  theme_void() +
  scale_fill_manual(values = c("Risor skin & Sponge" = "black", "Risor skin" = "#DDCC77")) +
  theme(legend.position = "none")  # Hide legend if not needed

ggsave("Donut_plot_skin_filtered.pdf", donut_plot.skin , width = 8, height = 6, useDingbats = FALSE)

#Donut plot - GUT

#Selecting sponge IDs
gut <- meta.all %>%
  filter(grepl("^RBEL", sample_id))%>%
  pull(sample_id)

#Subset physeq object specific for sponges
physeq_gut<- prune_samples(gut, physeq1.all.g.f)

gen_gut<- otu_table(physeq_gut)
gen_gut[gen_gut > 0] <- 1 #convert counts to presence/absence

#Calculate the proportion of gut samples each taxon is found in
gen_prev.gut<- rowSums(gen_gut) / length(gut)

gen_to_keep.01.gut<- taxa_names(physeq1.all.g.f)[gen_prev.gut >= 0.1] #at least 10% of risor skin

physeq1.all.g.01.gut<- prune_taxa(gen_to_keep.01.gut, physeq1.all.g.f) #Only 106 genera

physeq_gut_type = merge_samples(physeq1.all.g.01.gut, "type")
sample_data(physeq_gut_type)
table_physeq.g <- t(otu_table(physeq_gut_type, ))

table_physeq.g<- as.data.frame(table_physeq.g)
table_physeq.g<- rownames_to_column(table_physeq.g)
table_physeq.g<- table_physeq.g %>%
  dplyr::rename(ASV = rowname)

physeq_gut_tax<- tax_table(physeq1.all.g.01.gut)
physeq_gut_tax<- as.data.frame(physeq_gut_tax)
physeq_gut_tax<- rownames_to_column(physeq_gut_tax)
physeq_gut_tax<- physeq_gut_tax %>%
  dplyr::rename(ASV = rowname)

tax_table_gen.g<- physeq_gut_tax %>%
  inner_join(table_physeq.g, by = "ASV")%>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, pident, `Risor gut`)

filtered_tax_table_gen.gut.only<- tax_table_gen.g %>%
  filter(`Risor gut` > 50) #103 genus
#write.csv(filtered_tax_table_gen.gut.only, "skin_only_genus.csv")

tax_table_gen.all.gut.f  #shared between risor gut and sponge - 11

count_gut.only<- nrow(filtered_tax_table_gen.gut.only)
count_shared.g.s<- nrow(tax_table_gen.all.gut.f)

# Create a data frame for visualization
donut_data.gut <- data.frame(
  Category = c("Risor gut & Sponge", "Risor gut"),
  count = c(count_shared.g.s, count_gut.only - count_shared.g.s)
)

# Calculate percentages and cumulative values
donut_data.gut$fraction <- donut_data.gut$count / sum(donut_data.gut$count)
donut_data.gut$ymax <- cumsum(donut_data.gut$fraction)
donut_data.gut$ymin <- c(0, head(donut_data.gut$ymax, n=-1))
donut_data.gut$label_position <- (donut_data.gut$ymax + donut_data.gut$ymin) / 2  # Midpoint for text labels
donut_data.gut$label_text <- paste0(round(donut_data.gut$fraction * 100, 1), "%")  # Percentage text

# Create a donut plot - Figure 2c
donut_plot.gut <- ggplot(donut_data.gut, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
  geom_rect() +
  geom_text(aes(x = 3.5, y = label_position, label = label_text), size = 5, color = "white") +  # Add percentage labels
  coord_polar(theta = "y") +  # Convert to donut
  xlim(c(2, 5)) +  # Adjust for donut shape
  theme_void() +
  scale_fill_manual(values = c("Risor gut & Sponge" = "black", "Risor gut" = "#44AA99")) +
  theme(legend.position = "none")  # Hide legend if not needed

ggsave("Donut_plot_gut_filtered.pdf", donut_plot.gut , width = 8, height = 6, useDingbats = FALSE)


#### Key shared microbial genera (mirrored barplot) - Figure 4 ####

##Read filtered phyloseq object - 594 genera
physeq1.all.g.f

##Filtering by sponge abundance before comparison

#Selecting sponge IDs
sponge <- meta.all %>%
  filter(grepl("^IBEL", sample_id))%>%
  pull(sample_id)

#Subset physeq object specific for sponges
physeq_sponge<- prune_samples(sponge, physeq1.all.g.f)

#Compute presence/absence matrix for sponge samples
gen_sponge<- otu_table(physeq_sponge)
gen_sponge[gen_sponge > 0] <- 1 #convert counts to presence/absence

#Calculate the proportion of sponge samples each taxon is found in
gen_prev<- rowSums(gen_sponge) / length(sponge)

#Select genera present from main physeq object in which at least x% of sponge samples 

#For 10% of sponges
gen_to_keep.01<- taxa_names(physeq1.all.g.f)[gen_prev >= 0.1]

physeq1.all.g.01<- prune_taxa(gen_to_keep.01, physeq1.all.g.f) #Only 94 genera (97 when allowing for >10 reads)

###Agglomerate by sample type ###

#All
physeq1.all.type = merge_samples(physeq1.all.g.f, "type")
sample_data(physeq1.all.type)
table_physeq.all <- t(otu_table(physeq1.all.type, ))

#Only taxa present in at least 10% of sponges
physeq1.all.type.01 = merge_samples(physeq1.all.g.01, "type")
sample_data(physeq1.all.type.01)
table_physeq.all.01 <- t(otu_table(physeq1.all.type.01, ))

###saving df related to agglomerate###

#All
table_physeq1.all<- as.data.frame(table_physeq.all)
table_physeq1.all<- rownames_to_column(table_physeq1.all)
table_physeq1.all<- table_physeq1.all %>%
  dplyr::rename(ASV = rowname)

#Only taxa present in at least 10% of sponges
table_physeq1.all.01<- as.data.frame(table_physeq.all.01)
table_physeq1.all.01<- rownames_to_column(table_physeq1.all.01)
table_physeq1.all.01<- table_physeq1.all.01 %>%
  dplyr::rename(ASV = rowname)

#All
physeq1.tax.all<- tax_table(physeq1.all.g.f)
physeq1.tax.all<- as.data.frame(physeq1.tax.all)
physeq1.tax.all<- rownames_to_column(physeq1.tax.all)
physeq1.tax.all<- physeq1.tax.all %>%
  dplyr::rename(ASV = rowname)

#Only taxa present in at least 10% of sponges
physeq1.tax.all.01<- tax_table(physeq1.all.g.01)
physeq1.tax.all.01<- as.data.frame(physeq1.tax.all.01)
physeq1.tax.all.01<- rownames_to_column(physeq1.tax.all.01)
physeq1.tax.all.01<- physeq1.tax.all.01 %>%
  dplyr::rename(ASV = rowname)

#All
tax_table_gen.all<- physeq1.tax.all %>%
  inner_join(table_physeq1.all, by = "ASV")%>%
  dplyr::select(ASV, Phylum, Class, Order, Family, Genus, Species, pident, `Risor gut`, `Risor skin`, sponge, seawater, `E.prochilos gut`, `E.prochilos skin`)

tax_table_gen.all.ra<- tax_table_gen.all%>%
  mutate(
    Risor_gut_ra = `Risor gut` / sum(`Risor gut`),
    Sponge_ra = sponge / sum(sponge),
    Risor_skin_ra = `Risor skin` / sum(`Risor skin`)) 

#write.csv(tax_table_gen.all.ra, "all_Genera_table_ra.csv")

#Only taxa present in at least 10% of sponges
tax_table_gen.all.01<- physeq1.tax.all.01 %>%
  inner_join(table_physeq1.all.01, by = "ASV")%>%
  dplyr::select(ASV, Phylum, Class, Order, Family, Genus, Species, pident, `Risor gut`, `Risor skin`, sponge, seawater, `E.prochilos gut`, `E.prochilos skin`)

#Filtering out genus present in others than risor and sponge###

tax_table_gen.all.f <- tax_table_gen.all.01 %>%
  filter(seawater <= 0, `E.prochilos gut` <= 0, `E.prochilos skin` <= 0) #24 genus total

#### RISOR GUT AND SPONGE

#Filter taxa that is shared between risor gut and sponge
tax_table_gen.all.gut <- tax_table_gen.all.f %>%
  filter(`Risor gut` > 0 & sponge > 0) #14

#Filter taxa that is shared meaningful shared between risor gut and sponge (>50 reads)
tax_table_gen.all.gut.f <- tax_table_gen.all.gut %>%
  filter(`Risor gut` > 50 & sponge > 50) #11 genus

#write.csv(tax_table_gen.all.gut.f, "key_shared_gut.csv")

#Calculating relative abundance
aggregated_tax_table_gen.gut <- tax_table_gen.all.gut.f  %>%
  group_by(Genus, Phylum, Family) %>%  # Include Phylum in the grouping
  summarize(
    gut = sum(`Risor gut`, na.rm = TRUE),
    sponge = sum(sponge, na.rm = TRUE)
  ) %>%
  mutate(total_ab = gut + sponge) %>%
  ungroup()

aggregated_tax_table_gen.gut.ra <- aggregated_tax_table_gen.gut %>%
  mutate(
    `Risor gut` = gut / sum(gut),
    Sponge = sponge / sum(sponge)) %>%
  filter(!(Genus == "Cm1-21")) %>%
  filter(!(Genus == "Cohaesibacter"))

#write.csv(aggregated_tax_table_gen.gut.ra, "key_shared_gut2.csv")

#Plotting mirrored barplot

long_data_mirrored.gut <- aggregated_tax_table_gen.gut.ra %>%
  mutate(Family = ifelse(Family == "PAUC43f_marine_benthic_group", "PAUC43f", Family))%>%
  mutate(Family = ifelse(Family == "BD2-11_terrestrial_group", "BD2-11", Family))%>%
  mutate(Genus = ifelse(Genus == "PAUC43f_marine_benthic_group", "PAUC43f_marine", Genus))%>%
  mutate(Genus = ifelse(Genus == "BD2-11_terrestrial_group", "BD2-11_terrestrial", Genus))%>%
  mutate(Phylum_Genus = paste(Phylum, Family, Genus, sep = ":")) %>% # Combine Phylum and Genus
  arrange(total_ab) %>% # Sort by ra_diff
  mutate(Phylum_Genus = factor(Phylum_Genus, levels = unique(Phylum_Genus))) %>% # Factorize with sorted levels
  dplyr::select(Phylum_Genus, `Risor gut`, Sponge) %>%
  pivot_longer(cols = c(`Risor gut`, Sponge), names_to = "Sample_Type", values_to = "Occurrences") %>%
  mutate(Occurrences = if_else(Sample_Type == "Risor gut", -Occurrences, Occurrences)) # Mirror Risor skin

# Plotting the mirrored barplot - Figure 4b
shared_barplot.gut <- ggplot(long_data_mirrored.gut, aes(x = Occurrences, y = Phylum_Genus, fill = Sample_Type)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = abs, expand = c(0, 0)) +
  scale_fill_manual(values = c("#44AA99","#661100")) + 
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

### RISOR SKIN AND SPONGE

#Filter taxa that is shared between risor skin and sponge
tax_table_gen.all.f.skin <- tax_table_gen.all.f %>%
  filter(`Risor skin` > 0 & sponge > 0) #19 genus (20 when filtering for 10 reads)

#Filter taxa that is shared meaningful shared between risor skin and sponge (>50 reads)
tax_table_gen.all.f.skin<- tax_table_gen.all.f.skin %>%
  filter(`Risor skin` > 50 & sponge > 50) #17 genus

#write.csv(tax_table_gen.all.f.skin, "key_shared_skin.csv")

#Filter taxa that is bigger in the sponge than risor skin
tax_table_gen.all.f.skin <- tax_table_gen.all.f.skin %>%
  filter(sponge > `Risor skin`) #11 genus

#Calculating relative abundance

aggregated_tax_table_gen.skin <- tax_table_gen.all.f.skin %>%
  #select(-c(ASV, Class, Order, Species)) %>%
  group_by(Genus, Phylum, Family) %>%  #Include Phylum in the grouping
  summarize(
    skin = sum(`Risor skin`, na.rm = TRUE),
    sponge = sum(sponge, na.rm = TRUE),
  ) %>%
  mutate(total_ab = skin + sponge) %>%
  ungroup()

aggregated_tax_table_gen.skin.ra <- aggregated_tax_table_gen.skin %>%
  mutate(
    `Risor skin` = skin / sum(skin),
    Sponge = sponge / sum(sponge))%>%
  filter(!(Genus == "Microbulbifer"))%>%
  filter(!(Genus == "Fulvivirga"))

#write.csv(aggregated_tax_table_gen.skin.ra, "key_shared_skin2.csv")

#Plotting mirrored barplot - Figure 4

long_data_mirrored.skin <- aggregated_tax_table_gen.skin.ra %>%
  mutate(Phylum_Genus = paste(Phylum, Family, Genus, sep = ":")) %>% # Combine Phylum and Genus
  arrange(total_ab) %>% # Sort by ra_diff
  mutate(Phylum_Genus = factor(Phylum_Genus, levels = unique(Phylum_Genus))) %>% # Factorize with sorted levels
  select(Phylum_Genus, `Risor skin`, Sponge) %>%
  pivot_longer(cols = c(`Risor skin`, Sponge), names_to = "Sample_Type", values_to = "Occurrences") %>%
  mutate(Occurrences = if_else(Sample_Type == "Risor skin", -Occurrences, Occurrences)) # Mirror Risor skin

# Plotting the mirrored barplot - Figure 4a
shared_barplot.skin <- ggplot(long_data_mirrored.skin, aes(x = Occurrences, y = Phylum_Genus, fill = Sample_Type)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = abs, expand = c(0, 0)) +
  scale_fill_manual(values = c("#DDCC77","#661100")) + 
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
