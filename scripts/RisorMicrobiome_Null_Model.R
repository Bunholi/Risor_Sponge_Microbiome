library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(ggtext)
library(ggnewscale)
library(bipartite)
library(patchwork)

######################################################################
################# Risor microbiome project ###########################
########### DATA Analyses - Null Model (Figure 3) ###################
#####################################################################

setwd("User/workingdirectory") #your working directory

# Load data
asv_table <- read.table("genus.filtered.50ab.original.csv", header = TRUE, row.names = 1, sep = ",")
metadata <- read.table("metadata_filtered_all_5000.csv", header = TRUE, sep = ",")

#asv_table <- read.table("genus.filtered.50ab.5000.csv", header = TRUE, row.names = 1, sep = ",")
#metadata <- read.table("metadata_filtered_all_5000_nm.csv", header = TRUE, sep = ",")

# Transpose ASV table
asv_table <- t(asv_table)

# Check that rownames of ASV table match metadata
all(rownames(asv_table) %in% metadata$sample_id)  # Should be TRUE

# Binarize ASV data
asv_bin <- (asv_table > 0) * 1

# Compute shared taxa matrix (samples x samples)
shared_taxa_matrix <- asv_bin %*% t(asv_bin)

####Custom fully random null model####

#Fixed probability - minimal microbial sharing
ifallsame <- function(y){
  if (length(unique(y)) / ncol(y) == 1){
    if (y[1,1] == 1){
      y[sample(1:length(y), 1)] <- 0
    } else {
      y[sample(1:length(y), 1)] <- 1
    }
  }
  return (y)
}

#null - gives each cell a 50/50 chance of being a 1
nullmod <- function(y){
  for(i in 1:nrow(y)){ 
    for(j in 1:ncol(y)){ 
      y[i,j] <- sample(0:1,1)  # default p = 0.5
    }
  }
  ifallsame(y)
  return (y) 
}

nullmod <- function(y, prob = 0.5) {
  for (i in 1:nrow(y)) { 
    for (j in 1:ncol(y)) { 
      y[i, j] <- sample(0:1, 1, prob = c(1 - prob, prob))
    }
  }
  ifallsame(y)
  return(y)
}

# Generate 1000 null matrices
make.null <- function(x){
  l <- list()
  for (i in 1:1000){ 
    mat <- x
    l[[i]] <- nullmod(mat)
  }
  return (l)
}

make.null <- function(x, n = 1000, prob = 0.5) {
  l <- vector("list", n)
  for (i in 1:n) { 
    l[[i]] <- nullmod(x, prob = prob)
  }
  return(l)
}



#### Risor skin - Figure 3A ####

#Prepare metadata
skin_ids <- metadata %>% filter(type == "Risor skin") %>% pull(sample_id)
sponge_ids_skin <- metadata %>% filter(type == "sponge") %>% pull(sample_id)

# All combinations of fish and sponge pairs
fish_sponge_pairs_skin <- tidyr::crossing(
  Fish_skin = skin_ids,
  Sponge = sponge_ids_skin
)

# Use rowwise() to compute shared taxa for each pair
fish_sponge_pairs_skin <- fish_sponge_pairs_skin %>%
  rowwise() %>%
  mutate(SharedTaxa = shared_taxa_matrix[Fish_skin, Sponge]) %>%
  ungroup()

# Generate all unpaired combinations
fish_sponge_pairs_skin <- fish_sponge_pairs_skin %>%
  mutate(Link = ifelse(SharedTaxa >= 5, 1, 0))

#Build binary interaction matrix (Sponge as rows, Fish as columns)
interaction_matrix_skin <- fish_sponge_pairs_skin %>%
  group_by(Sponge, Fish_skin) %>%
  summarise(Link = max(Link), .groups = "drop") %>%
  pivot_wider(names_from = Fish_skin, values_from = Link, values_fill = list(Link = 0)) %>%
  column_to_rownames(var = "Sponge") %>%
  as.matrix()

#Compute observed metrics
obs_total_links_skin <- sum(interaction_matrix_skin)

#Wang's function
null_list_p10_skin <- make.null(interaction_matrix_skin, n = 1000, prob = 0.1)
null_list_p20_skin <- make.null(interaction_matrix_skin, n = 1000, prob = 0.2)
null_list_p30_skin <- make.null(interaction_matrix_skin, n = 1000, prob = 0.3)
null_list_p50_skin <- make.null(interaction_matrix_skin, n = 1000, prob = 0.5)

null_total_p10_skin <- sapply(null_list_p10_skin, sum)
p_val_p10 <- mean(null_total_p10_skin >= obs_total_links_skin)

null_total_p20_skin <- sapply(null_list_p20_skin, sum)
p_val_p20 <- mean(null_total_p20_skin >= obs_total_links_skin)

null_total_p30_skin <- sapply(null_list_p30_skin, sum)
p_val_p30 <- mean(null_total_p30_skin >= obs_total_links_skin)

null_total_p50_skin <- sapply(null_list_p50_skin, sum)
p_val_p50 <- mean(null_total_p50_skin >= obs_total_links_skin)

mean(null_total_links_skin)

#Plot with Wang's function - Figure 3A
plot_links_skin05 <- ggplot(data.frame(Null = null_total_p50_skin), aes(x = Null)) +
  geom_density(fill = "#DDCC77", alpha = 0.5) +
  geom_vline(xintercept = obs_total_links_skin, linetype = "dashed", color = "red") +
  labs(title = "R. ruber skin and sponge", x = "Number of Links (p = 0.5)", y = "Density")+
  annotate("text", x = obs_total_links_skin, y = 0,
           label = paste0("Observed = ", obs_total_links_skin,
                          "\nP = ", round(p_val_p50, 3)),
           hjust = 0, vjust = -1.5, fontface = "bold", size = 4, color = "black")+
  theme_bw()

plot_links_skin05 #Wang's function


#### Risor gut - Figure 3B ####

gut_ids <- metadata %>% filter(type == "Risor gut") %>% pull(sample_id)
sponge_ids <- metadata %>% filter(type == "sponge") %>% pull(sample_id)

# All combinations of fish and sponge pairs
fish_sponge_pairs_gut <- tidyr::crossing(
  Fish_gut = gut_ids,
  Sponge = sponge_ids
)

# Use rowwise() to compute shared taxa for each pair
fish_sponge_pairs_gut <- fish_sponge_pairs_gut %>%
  rowwise() %>%
  mutate(SharedTaxa = shared_taxa_matrix[Fish_gut, Sponge]) %>%
  ungroup()

# Generate all unpaired combinations
fish_sponge_pairs_gut <- fish_sponge_pairs_gut %>%
  mutate(Link = ifelse(SharedTaxa >= 5, 1, 0))

#Build binary interaction matrix (Sponge as rows, Fish as columns)
interaction_matrix_gut <- fish_sponge_pairs_gut %>%
  group_by(Sponge, Fish_gut) %>%
  summarise(Link = max(Link), .groups = "drop") %>%
  pivot_wider(names_from = Fish_gut, values_from = Link, values_fill = list(Link = 0)) %>%
  column_to_rownames(var = "Sponge") %>%
  as.matrix()

#Compute observed metrics
obs_total_links_gut <- sum(interaction_matrix_gut)
#Wang's function
#null_list_gut<- make.null(interaction_matrix_gut)
#null_total_links_gut <- sapply(null_list_gut, sum)
#p_val_links_gut <- mean(null_total_links_gut >= sum(interaction_matrix_gut))
#mean(null_total_links_gut)

null_list_p10_gut <- make.null(interaction_matrix_gut, n = 1000, prob = 0.1)
null_list_p20_gut <- make.null(interaction_matrix_gut, n = 1000, prob = 0.2)
null_list_p30_gut <- make.null(interaction_matrix_gut, n = 1000, prob = 0.3)
null_list_p50_gut <- make.null(interaction_matrix_gut, n = 1000, prob = 0.5)

null_total_p10_gut <- sapply(null_list_p10_gut, sum)
p_val_p10 <- mean(null_total_p10_gut >= obs_total_links_gut)

null_total_p20_gut <- sapply(null_list_p20_gut, sum)
p_val_p20 <- mean(null_total_p20_gut >= obs_total_links_gut)

null_total_p30_gut <- sapply(null_list_p30_gut, sum)
p_val_p30 <- mean(null_total_p30_gut >= obs_total_links_gut)

null_total_p50_gut <- sapply(null_list_p50_gut, sum)
p_val_p50 <- mean(null_total_p50_gut >= obs_total_links_gut)


#Plot with Wang's function
plot_links_gut05 <- ggplot(data.frame(Null = null_total_p30_gut), aes(x = Null)) +
  geom_density(fill = "#44AA99", alpha = 0.5) +
  geom_vline(xintercept = obs_total_links_gut, linetype = "dashed", color = "red") +
  labs(title = "R. ruber gut and sponge", x = "Number of Links (p = 0.5)", y = "Density")+
  annotate("text", x = obs_total_links_gut, y = 0,
           label = paste0("Observed = ", obs_total_links_gut,
                          "\nP = ", round(p_val_p50, 3)),
           hjust = 0, vjust = -1.5, fontface = "bold", size = 4, color = "black")+
  theme_bw()

plot_links_gut05 #Wang's function


#### E. prochilos skin - Figure 3C ####
# Get sample IDs
e.skin_ids <- metadata %>% filter(type == "E.prochilos skin") %>% pull(sample_id)
sponge_ids <- metadata %>% filter(type == "sponge") %>% pull(sample_id)

# All combinations of fish and sponge pairs
fish_sponge_pairs_e.skin <- tidyr::crossing(
  Fish_skin = e.skin_ids,
  Sponge = sponge_ids
)

# Use rowwise() to compute shared taxa for each pair
fish_sponge_pairs_e.skin <- fish_sponge_pairs_e.skin %>%
  rowwise() %>%
  mutate(SharedTaxa = shared_taxa_matrix[Fish_skin, Sponge]) %>%
  ungroup()

# Generate all unpaired combinations
fish_sponge_pairs_e.skin <- fish_sponge_pairs_e.skin %>%
  mutate(Link = ifelse(SharedTaxa >= 5, 1, 0))

#Build binary interaction matrix (Sponge as rows, Fish as columns)
interaction_matrix_e.skin <- fish_sponge_pairs_e.skin %>%
  group_by(Sponge, Fish_skin) %>%
  summarise(Link = max(Link), .groups = "drop") %>%
  pivot_wider(names_from = Fish_skin, values_from = Link, values_fill = list(Link = 0)) %>%
  column_to_rownames(var = "Sponge") %>%
  as.matrix()

#Compute observed metrics
obs_total_links_e.skin <- sum(interaction_matrix_e.skin)

#Wang's function
#null_list_e.skin<- make.null(interaction_matrix_e.skin)
#null_total_links_e.skin <- sapply(null_list_e.skin, sum)
#p_val_links_e.skin <- mean(null_total_links_e.skin >= sum(interaction_matrix_e.skin))
#mean(null_total_links_e.skin)

null_list_p10_e.skin <- make.null(interaction_matrix_e.skin, n = 1000, prob = 0.1)
null_list_p20_e.skin <- make.null(interaction_matrix_e.skin, n = 1000, prob = 0.2)
null_list_p30_e.skin <- make.null(interaction_matrix_e.skin, n = 1000, prob = 0.3)
null_list_p50_e.skin <- make.null(interaction_matrix_e.skin, n = 1000, prob = 0.5)

null_total_p10_e.skin <- sapply(null_list_p10_e.skin, sum)
p_val_p10 <- mean(null_total_p10_e.skin >= obs_total_links_e.skin)

null_total_p20_e.skin <- sapply(null_list_p20_e.skin, sum)
p_val_p20 <- mean(null_total_p20_e.skin >= obs_total_links_e.skin)

null_total_p30_e.skin<- sapply(null_list_p30_e.skin, sum)
p_val_p30 <- mean(null_total_p30_e.skin >= obs_total_links_e.skin)

null_total_p50_e.skin<- sapply(null_list_p50_e.skin, sum)
p_val_p50 <- mean(null_total_p50_e.skin >= obs_total_links_e.skin)


#Plot with Wang's function
plot_links_e.skin <- ggplot(data.frame(Null = null_total_p10_e.skin), aes(x = Null)) +
  geom_density(fill = "#B8887A", alpha = 0.5) +
  geom_vline(xintercept = obs_total_links_e.skin, linetype = "dashed", color = "red") +
  labs(title = "E. prochilos skin and sponge", x = "Number of Links (p = 0.1)", y = "Density")+
  annotate("text", x = obs_total_links_e.skin, y = 0,
           label = paste0("Observed = ", obs_total_links_e.skin,
                          "\nP = ", round(p_val_p10, 3)),
           hjust = 0, vjust = -1.5, fontface = "bold", size = 4, color = "black")+
  theme_bw()

plot_links_e.skin #Wang's function

ggsave("nullmodel_p01_E.skin.pdf", plot_links_e.skin, width = 10, height = 8)


#### E. prochilos gut - Figure 3D ####
# Get sample IDs
e.gut_ids <- metadata %>% filter(type == "E.prochilos gut") %>% pull(sample_id)
sponge_ids <- metadata %>% filter(type == "sponge") %>% pull(sample_id)

# All combinations of fish and sponge pairs
fish_sponge_pairs_e.gut <- tidyr::crossing(
  Fish_gut = e.gut_ids,
  Sponge = sponge_ids
)

# Use rowwise() to compute shared taxa for each pair
fish_sponge_pairs_e.gut <- fish_sponge_pairs_e.gut %>%
  rowwise() %>%
  mutate(SharedTaxa = shared_taxa_matrix[Fish_gut, Sponge]) %>%
  ungroup()

# Generate all unpaired combinations
fish_sponge_pairs_e.gut <- fish_sponge_pairs_e.gut %>%
  mutate(Link = ifelse(SharedTaxa >= 5, 1, 0))

#Build binary interaction matrix (Sponge as rows, Fish as columns)
interaction_matrix_e.gut <- fish_sponge_pairs_e.gut %>%
  group_by(Sponge, Fish_gut) %>%
  summarise(Link = max(Link), .groups = "drop") %>%
  pivot_wider(names_from = Fish_gut, values_from = Link, values_fill = list(Link = 0)) %>%
  column_to_rownames(var = "Sponge") %>%
  as.matrix()

#Compute observed metrics
obs_total_links_e.gut <- sum(interaction_matrix_e.gut)
#Wang's function

null_list_p10_e.gut <- make.null(interaction_matrix_e.gut, n = 1000, prob = 0.1)
null_list_p20_e.gut <- make.null(interaction_matrix_e.gut, n = 1000, prob = 0.2)
null_list_p30_e.gut <- make.null(interaction_matrix_e.gut, n = 1000, prob = 0.3)
null_list_p50_e.gut <- make.null(interaction_matrix_e.gut, n = 1000, prob = 0.5)

null_total_p10_e.gut <- sapply(null_list_p10_e.gut, sum)
p_val_p10 <- mean(null_total_p10_e.gut >= obs_total_links_e.gut)

null_total_p20_e.gut <- sapply(null_list_p20_e.skin, sum)
p_val_p20 <- mean(null_total_p20_e.gut >= obs_total_links_e.gut)

null_total_p30_e.gut<- sapply(null_list_p30_e.skin, sum)
p_val_p30 <- mean(null_total_p30_e.gut >= obs_total_links_e.gut)

null_total_p50_e.gut<- sapply(null_list_p50_e.skin, sum)
p_val_p50 <- mean(null_total_p50_e.gut >= obs_total_links_e.gut)

mean(null_total_links_e.gut)

#Plot with Wang's function
plot_links_e.gut <- ggplot(data.frame(Null = null_total_p10_e.gut), aes(x = Null)) +
  geom_density(fill = "darkslategrey", alpha = 0.5) +
  geom_vline(xintercept = obs_total_links_e.gut, linetype = "dashed", color = "red") +
  labs(title = "E. prochilos gut and sponge", x = "Number of Links (p = 0.1)", y = "Density")+
  annotate("text", x = obs_total_links_e.gut, y = 0,
           label = paste0("Observed = ", obs_total_links_e.gut,
                          "\nP = ", round(p_val_p10, 3)),
           hjust = 0, vjust = -1.5, fontface = "bold", size = 4, color = "black")+
  theme_bw()

plot_links_e.gut #Wang's function


#Joining plots - Risor only
combined_genus<- plot_links_skin / plot_links_gut

ggsave("nullmodel_Risor_all_genus.pdf", combined_genus, width = 10, height = 8)

risor_alternative_p_skin<- (plot_links_skin01 + plot_links_skin02) / (plot_links_skin03 + plot_links_skin05)

ggsave("nullmodel_Risor_alternative_p_skin.pdf", risor_alternative_p_skin, width = 10, height = 8)

risor_alternative_p_gut<- (plot_links_gut01 + plot_links_gut02) / (plot_links_gut03 + plot_links_gut05)

ggsave("nullmodel_Risor_alternative_p_gut.pdf", risor_alternative_p_gut, width = 10, height = 8)

#Joining plots - Elacatinus only

combined_alternative_p<- plot_links_e.skin / plot_links_e.gut
ggsave("nullmodel_prochilos_p01.pdf", combined_alternative_p, width = 10, height = 8)

#Joining plots - ALL
combined_genus_all<- (plot_links_skin + plot_links_e.skin) / (plot_links_gut + plot_links_e.gut)

ggsave("nullmodel_all_fish_genus.pdf", combined_genus_all, width = 10, height = 8)


