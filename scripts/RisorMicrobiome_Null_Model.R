######################################################################
################# Risor microbiome project ###########################
########### DATA Analyses - Null Model (Figure 3) ###################
#####################################################################

#Updated by: Ingrid Bunholi
#Date: 07/15/26

library(tidyverse)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(ggtext)
library(ggnewscale)
library(bipartite)
library(patchwork)

#### Environment set up and data loading ####
setwd("~/Desktop/Risor_microbiome_paper/Risor_Sponge_Microbiome/metadata") #your working directory

# Risor + sponge, bacterial ASVs
physeq.asv.bacteria.rs <- readRDS(
  "processed_objects/physeq.asv.bacteria.rs.rds")

physeq_null<- physeq.asv.bacteria.rs

rare_depth<- 5000

physeq_null_rare <- prune_samples(
  sample_sums(physeq_null) >= rare_depth,
  physeq_null
)

physeq_null_rare <- prune_taxa(taxa_sums(physeq_null_rare) > 0, physeq_null_rare)
#95 samples

sample_sums(physeq_null_rare)
set.seed(123)

physeq_null_rare <- rarefy_even_depth(
  physeq_null_rare,
  sample.size = rare_depth,
  rngseed = 123,
  replace = FALSE,
  trimOTUs = TRUE,
  verbose = TRUE
)

# Extract metadata
metadata <- data.frame(sample_data(physeq_null_rare)) %>%
  rownames_to_column("sample_id")

# Extract ASV table
asv_table <- as(otu_table(physeq_null_rare), "matrix")

# Make sure samples are rows and ASVs are columns
if (taxa_are_rows(physeq_null_rare)) {
  asv_table <- t(asv_table)
}

# Check that rownames match metadata sample IDs
all(rownames(asv_table) %in% metadata$sample_id)

# Binarize ASV data
asv_bin <- (asv_table > 0) * 1

# Compute shared ASV matrix: sample x sample
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

# Null model: each cell has probability = prob of being assigned 1
nullmod <- function(y, prob = 0.5) {
  for (i in 1:nrow(y)) {
    for (j in 1:ncol(y)) {
      y[i, j] <- sample(
        0:1,
        1,
        prob = c(1 - prob, prob)
      )
    }
  }
  ifallsame(y)
  return(y)
}

# Generate null matrices
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
p_val_p50_skin <- (sum(null_total_p50_skin >= obs_total_links_skin) + 1) /
  (length(null_total_p50_skin) + 1) #p = 0.00099


#Plot with Wang's function - Figure 3A
plot_links_skin05 <- ggplot(data.frame(Null = null_total_p50_skin), aes(x = Null)) +
  geom_density(fill = "#DDCC77", alpha = 0.5) +
  geom_vline(xintercept = obs_total_links_skin, linetype = "dashed", color = "red") +
  labs(x = "Number of Links (p = 0.5)", y = "Density")+
  annotate("text", x = obs_total_links_skin, y = 0,
           label = paste0("Observed = ", obs_total_links_skin,
                          "\nP = ", round(p_val_p50_skin, 3)),
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

p_val_p50_gut <- (sum(null_total_p50_gut >= obs_total_links_gut) + 1) /
  (length(null_total_p50_gut) + 1) #p = 0.00099


#Plot with Wang's function
plot_links_gut05 <- ggplot(data.frame(Null = null_total_p50_gut), aes(x = Null)) +
  geom_density(fill = "#44AA99", alpha = 0.5) +
  geom_vline(xintercept = obs_total_links_gut, linetype = "dashed", color = "red") +
  labs(x = "Number of Links (p = 0.5)", y = "Density")+
  annotate("text", x = obs_total_links_gut, y = 0,
           label = paste0("Observed = ", obs_total_links_gut,
                          "\nP = ", round(p_val_p50, 3)),
           hjust = 0, vjust = -1.5, fontface = "bold", size = 4, color = "black")+
  theme_bw()

plot_links_gut05 #Wang's function


#Checking metrics

range(interaction_matrix_gut)
table(as.vector(interaction_matrix_gut))

range(interaction_matrix_skin)
table(as.vector(interaction_matrix_skin))

nrow(interaction_matrix_skin) * ncol(interaction_matrix_skin)
sum(interaction_matrix_skin)

nrow(interaction_matrix_gut) * ncol(interaction_matrix_gut)
sum(interaction_matrix_gut)

summary(fish_sponge_pairs_skin$SharedTaxa)
summary(fish_sponge_pairs_gut$SharedTaxa)

quantile(fish_sponge_pairs_skin$SharedTaxa, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
quantile(fish_sponge_pairs_gut$SharedTaxa, probs = c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))

threshold_summary <- function(pair_data) {
  thresholds <- c(1, 2, 3, 5, 10)
  
  data.frame(
    threshold = thresholds,
    links = sapply(thresholds, function(th) sum(pair_data$SharedTaxa >= th)),
    total_pairs = nrow(pair_data)
  ) %>%
    mutate(proportion_links = links / total_pairs)
}

threshold_summary(fish_sponge_pairs_skin)
threshold_summary(fish_sponge_pairs_gut)

#Joining plots - Risor only
combined<- plot_links_skin05 / plot_links_gut05

ggsave("nullmodel_Risor_ASV.pdf", combined, width = 8, height = 3)

#Run for different probabilities

run_null_probs <- function(interaction_matrix,
                           tissue_name,
                           probs = c(0.1, 0.2, 0.3, 0.5),
                           n = 1000,
                           threshold = 5) {
  
  obs_links <- sum(interaction_matrix)
  
  results <- lapply(probs, function(p) {
    
    null_list <- make.null(interaction_matrix, n = n, prob = p)
    null_total <- sapply(null_list, sum)
    
    p_val <- (sum(null_total >= obs_links) + 1) /
      (length(null_total) + 1)
    
    data.frame(
      tissue = tissue_name,
      link_threshold = paste0(">=", threshold, " ASVs"),
      null_probability = p,
      observed_links = obs_links,
      mean_null_links = mean(null_total),
      sd_null_links = sd(null_total),
      min_null_links = min(null_total),
      max_null_links = max(null_total),
      empirical_p = p_val
    )
  })
  
  bind_rows(results)
}

skin_null_results <- run_null_probs(
  interaction_matrix_skin,
  tissue_name = "R. ruber skin",
  probs = c(0.1, 0.2, 0.3, 0.5),
  n = 1000,
  threshold = 5
)

gut_null_results <- run_null_probs(
  interaction_matrix_gut,
  tissue_name = "R. ruber gut",
  probs = c(0.1, 0.2, 0.3, 0.5),
  n = 1000,
  threshold = 5
)

null_sensitivity_results <- bind_rows(
  skin_null_results,
  gut_null_results
)

null_sensitivity_results

write.csv(null_sensitivity_results, "nullmodel_sensitivity.csv")
