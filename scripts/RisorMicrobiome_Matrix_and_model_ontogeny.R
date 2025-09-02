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

######################################################################
################# Risor microbiome project ###########################
########### DATA Analyses - Mechanisms and Ontogeny (Figure 5) #######
#####################################################################

setwd("User/workingdirectory") #your working directory

####Heatmap matrix - Figure S6 ####

# Read metadata table
df_meta <- read_csv("metadata_filtered_all_5000.csv")

# Read ASV table
df_gen <- read_csv("genus.filtered.50ab.original.csv")

#Normalize each column by its column sum (i.e., relative abundance)
df_gen_rel <- df_gen %>%
  column_to_rownames(var = colnames(df_gen)[1]) %>% # Temporarily set ASV IDs as rownames
  mutate(across(everything(), ~ . / sum(.))) %>%    # Normalize columns
  rownames_to_column(var = colnames(df_gen)[1])     # Restore ASV ID column

#if want counts instead of RA, use df_gen object instead

# Filter metadata by sample_id starting with 'IBEL'
df_meta <- df_meta %>% filter(str_starts(sample_id, "IBEL"))

#Extrac sample_id (sponge) and ref_skin or ref_gut
pairs <- df_meta %>% select(sample_id, ref_gut, ref_skin, risor_tl)

#SKIN and SPONGE
#Extract ASV values in pair

#Now use a function to loop through each pair and check if both columns from pairs are present in the asv table
#if present, pulls out those columns into lists and then combine all side by side (IBELX, SBELX, IBELY, SBELY,..)

#for both types
extract_asv_by_pairs <- function(df_gen_rel, pairs) {
  df_list <- list()
  
  for (i in 1:nrow(pairs)) {
    pair <- pairs[i, ]
    cols <- c(pair$sample_id, pair$ref_skin, pair$ref_gut)
    
    if (all(cols %in% colnames(df_gen_rel))) {
      df_instance <- df_gen_rel %>% select(all_of(cols))
      df_list[[length(df_list) + 1]] <- df_instance
    } else {
      cat("Pair not found in df_gen:", pair$sample_id, ",", pair$ref_skin, ",", pair$ref_gut, "\n")
    }
  }
  
  if (length(df_list) > 0) {
    return(bind_cols(df_list))
  } else {
    warning("No matching pairs found — result will be empty.")
    return(tibble())
  }
}

#for both types - more complete 
extract_asv_by_pairs <- function(df_gen_rel, pairs) {
  df_list <- list()
  valid_pairs <- list()
  
  for (i in 1:nrow(pairs)) {
    pair <- pairs[i, ]
    cols <- c(pair$sample_id, pair$ref_skin, pair$ref_gut)
    
    if (all(cols %in% colnames(df_gen_rel))) {
      df_instance <- df_gen_rel %>% select(all_of(cols))
      df_list[[length(df_list) + 1]] <- df_instance
      valid_pairs[[length(valid_pairs) + 1]] <- pair  # Save the valid pair
    } else {
      cat("Pair not found in df_gen:", pair$sample_id, ",", pair$ref_skin, ",", pair$ref_gut, "\n")
    }
  }
  
  df_combined <- if (length(df_list) > 0) bind_cols(df_list) else tibble()
  df_pairs_filtered <- if (length(valid_pairs) > 0) bind_rows(valid_pairs) else tibble()
  
  return(list(asv_data = df_combined, pairs_data = df_pairs_filtered))
}


result <- extract_asv_by_pairs(df_gen_rel, pairs)

asv_data <- result$asv_data        # The ASV table with matched columns
filtered_pairs <- result$pairs_data  # Only pairs where all 3 sample IDs were found

write.csv(filtered_pairs, "filtered_pairs.csv")

#fow just one type
extract_asv_by_pairs <- function(df_gen_rel, pairs) {
  df_list <- list()
  
  for (i in 1:nrow(pairs)) {
    pair <- pairs[i, ]
    cols <- c(pair$sample_id, pair$ref_gut)
    
    if (all(cols %in% colnames(df_gen_rel))) {
      df_instance <- df_gen_rel %>% select(all_of(cols))
      df_list[[length(df_list) + 1]] <- df_instance
    } else {
      cat("Pair not found in df_gen:", pair$sample_id, ",", pair$ref_gut, "\n")
    }
  }
  
  if (length(df_list) > 0) {
    return(bind_cols(df_list))
  } else {
    warning("No matching pairs found — result will be empty.")
    return(tibble())
  }
}


#Apply the function
df_gen_filtered <- extract_asv_by_pairs(df_gen_rel, pairs) 

# Count co-occurences between sample and skin columns 
#counts how many times both sponge and skin column have values >0 across rows and build a df from it
new_cols<- colnames(df_gen_filtered)
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
    sample_gen <- df_gen_filtered[[sample]]
    skin_gen <- df_gen_filtered[[skin]]
    
    count <- sum(pmin(sample_gen, skin_gen), na.rm = TRUE)
    
    data_skin <- add_row(data_skin,
                    sample_id = sample,
                    ref_skin = skin,
                    count = count)
  }
}

#Generating matrix

# Pivot to matrix format
df_matrix_skin <- data_skin %>%
  pivot_wider(names_from = ref_skin, values_from = count, values_fill = 0)

# Step 1: Prepare long-format version of confusion matrix
df_long_skin <- df_matrix_skin %>%
  pivot_longer(-sample_id, names_to = "ref_skin", values_to = "count")

# Step 3: Then mark diagonal AFTER ref_skin is reordered
df_long_skin <- df_long_skin %>%
  left_join(pairs %>% mutate(is_diagonal = TRUE), by = c("sample_id", "ref_skin")) %>%
  mutate(is_diagonal = ifelse(is.na(is_diagonal), FALSE, is_diagonal))

# Step 1: Create the ordered list of ref_skin based on risor_tl
ordered_ref_skin <- df_long_skin %>%
  filter(is_diagonal) %>%
  arrange(risor_tl) %>%
  pull(ref_skin)

# Step 2: Convert ref_skin to a factor with that custom order
df_long_skin$ref_skin <- factor(df_long_skin$ref_skin, levels = unique(ordered_ref_skin))

# Step 3: Replot using your ggplot code
matrix_skin_asv <- ggplot(df_long_skin, aes(x = ref_skin, y = sample_id, fill = count)) +
  geom_tile(color = "white") +
  geom_tile(data = df_long %>% filter(is_diagonal),
            aes(x = ref_skin, y = sample_id),
            fill = NA, color = "#DDCC77", linewidth = 1.2) +
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

ggsave("confusion_matrix_gen_skin_complete_pair_5000_RA_TLorganize.pdf", width = 8, height = 6, units = "in")

sum(apply(df_confusion, 1, which.max) == seq_len(nrow(df_confusion)))

#GUT and SPONGE

# Count co-occurences between sample and skin columns 
#counts how many times both sponge and skin column have values >0 across rows and build a df from it
new_cols<- colnames(df_gen_filtered)
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
    sample_gen <- df_gen_filtered[[sample]]
    gut_gen <- df_gen_filtered[[gut]]
    
    count <- sum(pmin(sample_gen, gut_gen), na.rm = TRUE)
    
    data_gut <- add_row(data_gut,
                    sample_id = sample,
                    ref_gut = gut,
                    count = count)
  }
}

#Confusion matrix 

# Pivot to confusion matrix format
df_matrix_gut <- data_gut %>%
  pivot_wider(names_from = ref_gut, values_from = count, values_fill = 0)

# Step 1: Prepare long-format version of confusion matrix
df_long_gut <- df_matrix_gut %>%
  pivot_longer(-sample_id, names_to = "ref_gut", values_to = "count")

# Step 3: Then mark diagonal AFTER ref_skin is reordered
df_long_gut <- df_long_gut %>%
  left_join(pairs %>% mutate(is_diagonal = TRUE), by = c("sample_id", "ref_gut"))%>%
  mutate(is_diagonal = ifelse(is.na(is_diagonal), FALSE, is_diagonal))

# Step 1: Create the ordered list of ref_skin based on risor_tl
ordered_ref_gut <- df_long_gut %>%
  filter(is_diagonal) %>%
  arrange(risor_tl) %>%
  pull(ref_gut)

# Step 2: Convert ref_skin to a factor with that custom order
df_long_gut$ref_gut <- factor(df_long_gut$ref_gut, levels = unique(ordered_ref_gut))

# Step 3: Replot using your ggplot code
matrix_gut_asv <- ggplot(df_long_gut, aes(x = ref_gut, y = sample_id, fill = count)) +
  geom_tile(color = "white") +
  geom_tile(data = df_long_gut %>% filter(is_diagonal),
            aes(x = ref_gut, y = sample_id),
            fill = NA, color = "#44AA99", linewidth = 1.2) +
  geom_text(aes(label = round(count, 3)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "blue") +
  coord_fixed() +
  labs(
    x = "Corresponding Risor gut ID (ordered by TL)",
    y = "Sponge ID"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

matrix_gut_asv

ggsave("confusion_matrix_gen_gut_complete_pair_5000_RA_TLorganize.pdf", width = 8, height = 6, units = "in")


#### Bayesian beta regression - relative abundance vs body size - Figure 5####

### SKIN MICROBIOME
df_long_skin <- df_long_skin %>%
  group_by(ref_skin) %>%
  mutate(risor_tl = unique(na.omit(risor_tl))[1]) %>%
  ungroup()

hist(df_long_skin$count) #check distribution - beta regression (0-1)

df_long_skin <- df_long_skin %>%
  mutate(ref_skin = as.factor(ref_skin))%>%
  filter(count > 0) #remove zero values

#Run bayesian model with beta distribution
skin.tl <- brm(count ~ risor_tl + (1 | ref_skin), 
              family = Beta(),
              data = df_long_skin,
              seed = 123)

summary(skin.tl)

#Create a new dataset across the range of risor_tl
newdata.skin <- expand_grid(
  risor_tl = seq(min(df_long_skin$risor_tl), max(df_long_skin$risor_tl), length.out = 100),
  ref_skin = unique(df_long_skin$ref_skin)
)

#Sample expected predictions from the posterior (500 draws)
posterior_lines.skin <- add_epred_draws(
  newdata.skin,
  object = skin.tl,
  re_formula = NA,
  ndraws = 500
)

summary(posterior_lines.skin$.epred)

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

ggsave("RA_skin_TL_model.pdf", width = 8, height = 6, units = "in")

### GUT MICROBIOME

df_long_gut <- df_long_gut %>%
  group_by(ref_gut) %>%
  mutate(risor_tl = unique(na.omit(risor_tl))[1]) %>%
  ungroup()

hist(df_long_gut$count) #check distribution - beta regression (0-1)

df_long_gut <- df_long_gut %>%
  mutate(ref_gut = as.factor(ref_gut))

#Run bayesian model with beta distribution
gut.tl <- brm(count ~ risor_tl + (1 | ref_gut), 
              family = Beta(),
              data = df_long_gut,
              seed = 123)

summary(gut.tl)

gut.tl <- brm(count ~ risor_tl + (1|ref_gut),
              df_long_gut,
             iter = 10000,
             chains = 5,
             family = Beta(),
             seed = 123)

summary(gut.tl)

plot(gut.tl)
conditional_effects(gut.tl, spaghetti = T, ndraws = 500)

#Create a new dataset across the range of risor_tl
newdata <- expand_grid(
  risor_tl = seq(min(df_long_gut$risor_tl), max(df_long_gut$risor_tl), length.out = 100),
  ref_gut = unique(df_long_gut$ref_gut)
)

#Sample expected predictions from the posterior (500 draws)
posterior_lines <- add_epred_draws(
  newdata,
  object = gut.tl,
  re_formula = NA,
  ndraws = 500
)

summary(posterior_lines$.epred)

#Plot spaghetti + raw data

#plotting model's predicted mean (.epred) - from the posterior 
#plotting observed data - black dots -  Figure 5B

gut_ont<-ggplot(posterior_lines, aes(x = risor_tl, y = .epred, group = .draw)) +
  geom_line(alpha = 0.1, color = "#44AA99") +  # spaghetti
  geom_point(data = df_long_gut, aes(x = risor_tl, y = count), inherit.aes = FALSE, alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Risor ruber - Total length (mm)",
    y = "Relative abundance shared with sponge"
  )

ggsave("RA_gut_TL_model.pdf", width = 8, height = 6, units = "in")
