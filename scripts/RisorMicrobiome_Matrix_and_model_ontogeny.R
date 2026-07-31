######################################################################
################# Risor microbiome project ###########################
########### DATA Analyses - Mechanisms and Ontogeny (Figure 5) #######
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
library(glmmTMB)
library(brms)
library(ggplot2)
library(tidybayes)

#### Environment set up and data loading ####
setwd("~/Desktop/Risor_microbiome_paper/Risor_Sponge_Microbiome/metadata") #your working directory

#Loading processed phyloseq object
physeq.asv.bacteria.rs <- readRDS(
  "processed_objects/physeq.asv.bacteria.rs.rds")

physeq_rs<- physeq.asv.bacteria.rs

#sequencing depth 
physeq_rs <- prune_samples(
  sample_sums(physeq_rs) >= 1000,
  physeq_rs)

physeq_rs <- prune_taxa(
  taxa_sums(physeq_rs) > 0,
  physeq_rs)

#Remove Cyanobacteria

tax_df <- as.data.frame(tax_table(physeq_rs))

is_cyanobacteria <- !is.na(tax_df$Phylum) &
  grepl("Cyanobacter", tax_df$Phylum, ignore.case = TRUE)

physeq_rs<- prune_taxa(
  taxa_names(physeq_rs)[!is_cyanobacteria],
  physeq_rs)

physeq_rs <- prune_taxa(
  taxa_sums(physeq_rs) > 0,
  physeq_rs)

physeq_rs #10692 ASVs and 121 samples (only 395 ASVs removed)

#Extract metadata directly from phyloseq object
df_meta <- data.frame(sample_data(physeq_rs)) %>%
  rownames_to_column("sample_name")

#Create sponge-skin pairs:
sponge_pairs <- df_meta %>%
  filter(type == "sponge") %>%
  transmute(
    sample_id = sample_name,        
    ref_skin = as.character(ref_skin), #corresponding SBEL skin
    risor_id,
    sponge_id,
    risor_tl)

#Create sponge-gut pairs:
gut_pairs <- df_meta %>%
  filter(type == "Risor gut") %>%
  transmute(
    ref_skin = as.character(ref_skin), #corresponding SBEL skin
    risor_id,
    sponge_id,
    ref_gut = sample_name)

#join them using the shared skin ID and individual identifier:
pairs <- sponge_pairs %>%
  left_join(
    gut_pairs,
    by = c("ref_skin", "risor_id", "sponge_id"))

#keeping only complete triplets
pairs <- pairs %>%
  filter(
    !is.na(ref_gut),
    sample_id %in% sample_names(physeq_rs),
    ref_skin %in% sample_names(physeq_rs),
    ref_gut %in% sample_names(physeq_rs))

#organize by body size
pairs <- pairs %>%
  mutate(
    sample_id = as.character(sample_id),
    ref_skin = as.character(ref_skin),
    ref_gut = as.character(ref_gut),
    risor_tl = as.numeric(risor_tl)
  ) %>%
  distinct(sample_id, ref_skin, ref_gut, .keep_all = TRUE) %>%
  arrange(risor_tl)

#write_tsv(pairs, "triplet_pairs_microbiome.tsv")

#then restrict the phyloseq object
samples_to_keep <- unique(c(
  pairs$sample_id,
  pairs$ref_skin,
  pairs$ref_gut))

physeq_pairs <- prune_samples(
  sample_names(physeq_rs) %in% samples_to_keep,
  physeq_rs)

physeq_pairs <- prune_taxa(
  taxa_sums(physeq_pairs) > 0,
  physeq_pairs) #8153 ASVs and 96 samples


# Keep only R. ruber gut samples
gut_ps <- subset_samples(
  physeq_pairs,
  type == "Risor gut"
)

# Keep Epulopiscium
epulo_ps <- prune_taxa(
  as.character(tax_table(gut_ps)[, "Genus"]) == "Epulopiscium",
  gut_ps
)

# Presence in each gut sample
epulo_presence <- sample_sums(epulo_ps) > 0

# Summary
data.frame(
  gut_samples = length(epulo_presence),
  samples_with_Epulopiscium = sum(epulo_presence),
  prevalence_percent = mean(epulo_presence) * 100
)

#checking read depth
read_counts.pairs <- tibble(
  sample = sample_names(physeq_pairs),
  reads = sample_sums(physeq_pairs)
)

pairs_reads <- pairs %>%
  left_join(
    read_counts.pairs,
    by = c("sample_id" = "sample")
  ) %>%
  dplyr::rename(sponge_reads = reads) %>%
  
  left_join(
    read_counts.pairs,
    by = c("ref_skin" = "sample")
  ) %>%
  dplyr::rename(skin_reads = reads) %>%
  
  left_join(
    read_counts.pairs,
    by = c("ref_gut" = "sample")
  ) %>%
  dplyr::rename(gut_reads = reads)

#physeq_pairs <- tax_glom(
 # physeq_pairs ,
  #taxrank = "Genus",
  #NArm = FALSE) #1212 genus and 96 samples

#convert counts to relative abundance within phyloseq

physeq_asv_rel <- transform_sample_counts(
  physeq_pairs,
  function(x) x / sum(x))

#extract ASV table
asv_matrix <- as(otu_table(physeq_asv_rel), "matrix")

if (!taxa_are_rows(physeq_asv_rel)) {
  asv_matrix <- t(asv_matrix)}

df_asv_rel <- as.data.frame(asv_matrix) %>%
  rownames_to_column("ASV")

####Helper function####

#Now use a function to loop through each pair and check if both columns from pairs are present in the asv table
#if present, pulls out those columns into lists and then combine all side by side (IBELX, SBELX, IBELY, SBELY,..)
extract_asv_by_pairs <- function(df_rel, pairs) {
  
  df_list <- list()
  valid_pairs <- list()
  
  for (i in seq_len(nrow(pairs))) {
    
    cols <- c(
      pairs$sample_id[i],
      pairs$ref_skin[i],
      pairs$ref_gut[i]
    )
    
    if (all(cols %in% colnames(df_rel))) {
      
      df_list[[length(df_list) + 1]] <- df_rel %>%
        dplyr::select(all_of(cols))
      
      valid_pairs[[length(valid_pairs) + 1]] <- pairs[i, ]
      
    } else {
      
      message(
        "Missing sample(s) for pair: ",
        paste(cols[!cols %in% colnames(df_rel)], collapse = ", ")
      )
    }
  }
  
  df_combined <- if (length(df_list) > 0) {
    bind_cols(df_list)
  } else {
    tibble()
  }
  
  df_pairs_filtered <- if (length(valid_pairs) > 0) {
    bind_rows(valid_pairs)
  } else {
    tibble()
  }
  
  list(
    asv_data = df_combined,
    pairs_data = df_pairs_filtered
  )
}

####Matrix pairing relative abundance shared - Figure S3####

result_asv <- extract_asv_by_pairs(
  df_rel = df_asv_rel,
  pairs = pairs)

df_asv_filtered <- result_asv$asv_data
filtered_pairs <- result_asv$pairs_data %>%
  arrange(risor_tl)

sponge_order <- filtered_pairs$sample_id
skin_order   <- filtered_pairs$ref_skin
gut_order    <- filtered_pairs$ref_gut

# Count co-occurences between sample and skin columns 
#counts how many times both sponge and skin column have values >0 across rows and build a df from it
new_cols<- colnames(df_asv_filtered)
# Extract sample and skin columns by prefix
sample_new_list <- new_cols[startsWith(new_cols, "IBEL")]
skin_new_list <- new_cols[startsWith(new_cols, "SBEL")]

# Initialize output tibble
data_skin <- tibble(sample_id = character(),
               ref_skin = character(),
               count = integer())

# Count co-occurrence of > 0 values between each sample/skin pair
for (sample in sample_new_list) {
  for (skin in skin_new_list) {
    sample_gen <- df_asv_filtered[[sample]]
    skin_gen <- df_asv_filtered[[skin]]
    
    count <- sum(pmin(sample_gen, skin_gen), na.rm = TRUE)
    
    data_skin <- add_row(data_skin,
                    sample_id = sample,
                    ref_skin = skin,
                    count = count)
  }
}

#Generating matrix

#Pivot to matrix format
df_matrix_skin <- data_skin %>%
  pivot_wider(names_from = ref_skin, values_from = count, values_fill = 0)

#Prepare long-format version of confusion matrix
df_long_skin <- df_matrix_skin %>%
  pivot_longer(-sample_id, names_to = "ref_skin", values_to = "count")

#Then mark diagonal AFTER ref_skin is reordered
df_long_skin <- df_long_skin %>%
  left_join(pairs %>% mutate(is_diagonal = TRUE), by = c("sample_id", "ref_skin")) %>%
  mutate(is_diagonal = ifelse(is.na(is_diagonal), FALSE, is_diagonal))

#Create the ordered list of ref_skin based on risor_tl
ordered_ref_skin <- df_long_skin %>%
  filter(is_diagonal) %>%
  arrange(risor_tl) %>%
  pull(ref_skin)

#Convert ref_skin to a factor with that custom order
df_long_skin$ref_skin <- factor(df_long_skin$ref_skin, levels = unique(ordered_ref_skin))

df_long_skin <- df_long_skin %>%
  mutate(
    ref_skin = factor(
      as.character(ref_skin),
      levels = skin_order
    ),
    sample_id = factor(
      as.character(sample_id),
      levels = rev(sponge_order)
    )
  )


#Plotting
matrix_skin_asv <- ggplot(df_long_skin, aes(x = ref_skin, y = sample_id, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(count, 3)), size = 1.5) +
  scale_fill_gradient(low = "white", high = "blue") +
  coord_fixed() +
  labs(
    x = "Corresponding Risor skin ID (ordered by TL)",
    y = "Sponge ID"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

matrix_skin_asv

ggsave("confusion_matrix_skin_ASV.pdf", width = 8, height = 6, units = "in")

#GUT and SPONGE

# Count co-occurences between sample and skin columns 
#counts how many times both sponge and skin column have values >0 across rows and build a df from it
new_cols<- colnames(df_asv_filtered)
# Extract sample and skin columns by prefix
sample_new_list <- new_cols[startsWith(new_cols, "IBEL")]
gut_new_list <- new_cols[startsWith(new_cols, "RBEL")]

# Initialize output tibble
data_gut <- tibble(sample_id = character(),
               ref_gut = character(),
               count = integer())

# Count co-occurrence of > 0 values between each sample/skin pair
for (sample in sample_new_list) {
  for (gut in gut_new_list) {
    sample_gen <- df_asv_filtered[[sample]]
    gut_gen <- df_asv_filtered[[gut]]
    
    count <- sum(pmin(sample_gen, gut_gen), na.rm = TRUE)
    
    data_gut <- add_row(data_gut,
                    sample_id = sample,
                    ref_gut = gut,
                    count = count)
  }
}

#Confusion matrix 

#Pivot to confusion matrix format
df_matrix_gut <- data_gut %>%
  pivot_wider(names_from = ref_gut, values_from = count, values_fill = 0)

#Prepare long-format version of confusion matrix
df_long_gut <- df_matrix_gut %>%
  pivot_longer(-sample_id, names_to = "ref_gut", values_to = "count")

#Then mark diagonal AFTER ref_skin is reordered
df_long_gut <- df_long_gut %>%
  left_join(pairs %>% mutate(is_diagonal = TRUE), by = c("sample_id", "ref_gut"))%>%
  mutate(is_diagonal = ifelse(is.na(is_diagonal), FALSE, is_diagonal))

#Create the ordered list of ref_skin based on risor_tl
ordered_ref_gut <- df_long_gut %>%
  filter(is_diagonal) %>%
  arrange(risor_tl) %>%
  pull(ref_gut)

#Convert ref_skin to a factor with that custom order
df_long_gut$ref_gut <- factor(df_long_gut$ref_gut, levels = unique(ordered_ref_gut))

df_long_gut <- df_long_gut %>%
  mutate(
    ref_gut = factor(
      as.character(ref_gut),
      levels = gut_order
    ),
    sample_id = factor(
      as.character(sample_id),
      levels = rev(sponge_order)
    )
  )

#Plot
matrix_gut_asv <- ggplot(df_long_gut, aes(x = ref_gut, y = sample_id, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(count, 3)), size = 1.5) +
  scale_fill_gradient(low = "white", high = "blue") +
  coord_fixed() +
  labs(
    x = "Corresponding Risor gut ID (ordered by TL)",
    y = "Sponge ID"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

matrix_gut_asv

ggsave("confusion_matrix_gut_ASV.pdf", width = 8, height = 6, units = "in")


#### Bayesian beta regression - relative abundance vs body size - Figure 5####

### SKIN MICROBIOME

df_skin_model <- df_long_skin %>%
  mutate(
    ref_skin = factor(ref_skin),
    sample_id = factor(sample_id)
  ) %>%
  filter(
    !is.na(count),
    !is.na(risor_tl),
    count > 0,
    count < 1
  )

range(df_skin_model$count)

sum(df_skin_model$count == 0)
sum(df_skin_model$count == 1)


skin.tl <- brm(
  count ~ risor_tl,
  family = Beta(link = "logit"),
  data = df_skin_model,
  chains = 5,
  iter = 10000,
  warmup = 5000,
  cores = 5,
  seed = 123,
  control = list(adapt_delta = 0.99)
)

summary(skin.tl)

newdata.skin <- tibble(
  risor_tl = seq(
    min(df_skin_model$risor_tl),
    max(df_skin_model$risor_tl),
    length.out = 100
  )
)

posterior_lines.skin <- add_epred_draws(
  newdata.skin,
  object = skin.tl,
  re_formula = NA,
  ndraws = 500
)

#Plot spaghetti + raw data
#plotting model's predicted mean (.epred) - from the posterior 
#plotting observed data - black dots - Figure 5A

skin_ont<-ggplot(posterior_lines.skin, aes(x = risor_tl, y = .epred, group = .draw)) +
  geom_line(alpha = 0.1, color = "#DDCC77") +  # spaghetti
  geom_point(data = df_long_skin, aes(x = risor_tl, y = count), inherit.aes = FALSE, alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Risor ruber - Total length (mm)",
    y = "Relative abundance shared with sponge"
  )

ggsave("RA_skin_TL_model_ASV.pdf", width = 6, height = 4, units = "in")

### GUT MICROBIOME

df_gut_model <- df_long_gut %>%
  mutate(
    ref_gut = factor(ref_gut),
    sample_id = factor(sample_id)
  ) %>%
  filter(
    !is.na(count),
    !is.na(risor_tl),
    count > 0,
    count < 1
  )

range(df_gut_model$count)

sum(df_skin_model$count == 0)
sum(df_skin_model$count == 1)


gut.tl <- brm(
  count ~ risor_tl,
  family = Beta(link = "logit"),
  data = df_gut_model,
  chains = 5,
  iter = 10000,
  warmup = 5000,
  cores = 5,
  seed = 123,
  control = list(adapt_delta = 0.99)
)

summary(gut.tl)

newdata.gut <- tibble(
  risor_tl = seq(
    min(df_gut_model$risor_tl),
    max(df_gut_model$risor_tl),
    length.out = 100
  )
)

posterior_lines.gut <- add_epred_draws(
  newdata.gut,
  object = gut.tl,
  re_formula = NA,
  ndraws = 500
)

#Plot spaghetti + raw data

#plotting model's predicted mean (.epred) - from the posterior 
#plotting observed data - black dots - Figure 5B

gut_ont<-ggplot(posterior_lines.gut, aes(x = risor_tl, y = .epred, group = .draw)) +
  geom_line(alpha = 0.1, color = "#44AA99") +  # spaghetti
  geom_point(data = df_long_gut, aes(x = risor_tl, y = count), inherit.aes = FALSE, alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Risor ruber - Total length (mm)",
    y = "Relative abundance shared with sponge"
  )

ggsave("RA_gut_TL_model_ASV.pdf", width = 6, height = 4, units = "in")


#Sequencing depth adjusted model - alternative

gut_depth_check <- df_gut_model %>%
  left_join(
    pairs_reads %>%
      dplyr::select(
        sample_id,
        ref_gut,
        gut_reads,
        sponge_reads
      ),
    by = c("sample_id", "ref_gut")
  )

gut.tl.depth <- brm(
  count ~ risor_tl + log10(gut_reads),
  family = Beta(link = "logit"),
  data = gut_depth_check,
  chains = 5,
  iter = 10000,
  warmup = 5000,
  cores = 5,
  seed = 123,
  control = list(adapt_delta = 0.99)
)

fixef(gut.tl)["risor_tl", ]

fixef(gut.tl.depth)["risor_tl", ]

fixef(gut.tl.depth)

#body-size effect is attenuated after accounting for gut sequencing depth, but it remains clearly negative


