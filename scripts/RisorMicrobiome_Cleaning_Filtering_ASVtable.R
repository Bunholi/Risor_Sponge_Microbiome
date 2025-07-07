############################################################
################# Risor microbiome project #################
####### Cleaning and filtering outputs from QIIME2 #########
###########################################################

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
setwd("User/workingdirectory") #your working directory

#Reading metadata
meta.all<- read_tsv("mdat_all.tsv") #includes all sample types

#Reading ASV table
asv.all<-read_qza("filtered-table-all.qza") #includes all sample types
asv.all.df<- as.data.frame(asv.all$data) #transforming qza into df

asv.all.df<- asv.all.df%>%
  rownames_to_column(var = "Feature.ID")
seq.all<- read_tsv("seq-all.tsv") #sequences from each ASV

#Reading taxonomy
taxon.all<- read_qza("taxonomy-all.qza") 

taxon.all.df<- as.data.frame(taxon.all$data)

taxon.all.df <- taxon.all.df %>% #organizing dataframe output
  separate(Taxon, into =c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species"), sep= ";\\s*")%>%
  mutate(Domain = str_replace(Domain, pattern="d__", replacement="")) %>%
  mutate(Phylum = str_replace(Phylum, pattern="p__", replacement="")) %>%
  mutate(Class = str_replace(Class, pattern="c__", replacement="")) %>%
  mutate(Order = str_replace(Order, pattern="o__", replacement="")) %>%
  mutate(Family = str_replace(Family, pattern="f__", replacement="")) %>%
  mutate(Genus = str_replace(Genus, pattern="g__", replacement="")) %>%
  mutate(Species = str_replace(Species, pattern="s__", replacement=""))%>%
  filter(!(Domain == "Unassigned")) %>% #filtering out Unassigned ASVs
  filter(!(Feature.ID == "4ab56f9b012674eb6981ba9d3ef36170"))%>% #filtering out ASVs amplified in in blanks
  filter(!(Feature.ID == "7ac9047b8f42a496b1640c20ea18156a"))%>%
  filter(!(Feature.ID == "94cb34fb244cbbfaeaafc9599a8f986b"))%>%
  filter(!(Feature.ID == "ac164649e33237f799c1ffee4f344d5d")) 

taxon.all.df.f<- taxon.all.df %>%
  semi_join(asv.all.df, by = "Feature.ID") #keeping only ASV in taxonomy file that are also present in ASV table

#Loading top hits file from qiime2 to obtain taxonomic assingment information  
top<- read_tsv("tophits-all.tsv")
top<- top %>%
  filter(pident >=70) #filtering only stuff with pident >70

top.f <- top %>% #Keeping informating from the ASV with the highest pident for each Feature.ID
  group_by(Feature.ID) %>%
  filter(row_number(desc(pident)) == 1) %>%
  ungroup()


#### 2 - Filtration process ####

#Filtering out ASVs present in blank samples
blank_check <- asv.all.df %>% summarise(across(everything(), ~ any(. > 0)))
blank_check_sponge<- asv.all.df %>% filter(IBELBlank2 !=0)
blank_check_skin<- asv.all.df %>% filter(SBELBlank !=0)
blank_check_gut<- asv.all.df %>% filter(RBELBlank1 !=0)

asv.all.df.f <- asv.all.df %>% #for each blank column, subtract from all other columns row-wise
  mutate(across(-c(Feature.ID, IBELBlank2), ~ pmax(0, . - IBELBlank2)))

asv.all.df.f <- asv.all.df.f %>%
  mutate(across(-c(Feature.ID, SBELBlank), ~ pmax(0, . - SBELBlank)))

asv.all.df.f <- asv.all.df.f %>%
  mutate(across(-c(Feature.ID, RBELBlank1), ~ pmax(0, . - RBELBlank1)))

asv.all.df.f<- asv.all.df %>%
  filter(!(Feature.ID == "4ab56f9b012674eb6981ba9d3ef36170"))%>% #Removing from previous ASV table
  filter(!(Feature.ID == "7ac9047b8f42a496b1640c20ea18156a"))%>%
  filter(!(Feature.ID == "94cb34fb244cbbfaeaafc9599a8f986b"))%>%
  filter(!(Feature.ID == "ac164649e33237f799c1ffee4f344d5d"))%>%
  filter(Feature.ID %in% taxon.all.df.f$Feature.ID)

##Cleaning low coverage ASVs - Summing up the occurrences in each individual/column to remove the ones with less than 2k

asv.all.df.n = asv.all.df.f %>% #185 individuals
  mutate_at(vars(2:186), as.numeric)
asv.all.sums <- colSums(asv.all.df.n[2:186])

indv_to_rem <- asv.all.sums < 2000

indv_to_rem_df <- data.frame(
  ID = names(asv.all.sums),
  Total_Sum = asv.all.sums,
  row.names = NULL,
  Remove = indv_to_rem ) %>%
  filter(Remove == TRUE)
#write.csv(indv_to_rem_df, "indv_less_2000_all_final.csv")

indv_read_sum <- data.frame(
  ID = names(asv.all.sums),
  Total_Sum = asv.all.sums,
  row.names = NULL)
#write.csv(indv_read_sum, "indv_read_sum_all.csv")

asv.all.df.f<- asv.all.df.n %>%
  select(-c(indv_to_rem_df$ID)) #Filtering out removed individuals from the main table
#write.csv(asv.all.df.f, "ASV_df_all_2000_final.csv") #165 total individuals

asv.all.df<- asv.all.df.f%>% 
  column_to_rownames(var = "Feature.ID")

asv.all.df.t<- as.data.frame(t(asv.all.df))
asv.all.df.t<- asv.all.df.t%>%
  rownames_to_column(var = "sample_id")

meta.all.f<- meta.all %>% #Filtering out removed individuals from the metadata
  filter(sample_id %in% asv.all.df.t$sample_id)
#write.csv(meta.all.f, "metadata_filtered_all_5000.csv") #saving new metadata

rem_indv<- meta.all %>%
  filter(!(sample_id %in% asv.all.df.t$sample_id)) #all samples removed
#write.csv(rem_indv, "metadata_remv_2000.csv") #saving samples removed

##Merging ASV, taxonomy, and top hits tables after cleaning
asv.taxon.all.f<- asv.all.df.f %>%
  inner_join(taxon.all.df.f, by = "Feature.ID")

asv.taxon.all.f.2<- asv.taxon.all.f %>%
  inner_join(top.f, by = "Feature.ID") %>%
  dplyr::rename(ASV = Feature.ID)
#write.csv(asv.taxon.g.f.2, "asv.taxon.g.f.2.csv")

##### 3 - Filtering final taxonomy based on percent identity #####
asv.taxon.all.match<- asv.taxon.all.f.2 %>%
  mutate(asv.species = case_when(!is.na(Species) & pident >= 97 ~ Species,
                                 TRUE ~ "NA"),
         asv.genus = case_when(!is.na(Genus) & pident >= 97 ~ Genus,
                               TRUE ~ "NA"),
         asv.family = case_when(!is.na(Family) & pident >= 85 ~ Family,
                                TRUE ~ "NA"),
         asv.order = case_when(!is.na(Order) & pident >= 85 ~ Order,
                               TRUE ~ "NA"),
         asv.class = case_when(!is.na(Class) & pident >= 85 ~ Class,
                               TRUE ~ "NA"),
         asv.phylum = case_when(!is.na(Phylum) & pident >= 70 ~ Phylum,
                                TRUE ~ "NA"))%>%
  filter(!if_all(c(asv.species, asv.genus, asv.family, asv.order, asv.class, asv.phylum), ~ . == "NA"))


asv.taxon.all.final<- asv.taxon.all.match %>%
  left_join(seq.all, by="ASV")
#write.csv(asv.taxon.all.final, "asv.taxon.all.match.filtered.5000.csv") #14499 ASVs

##### 4 - Calculating and plotting rarefaction curve #####

#Checking the number of ASVs per individual ID
total.s.asv <- asv.taxon.s.f %>%
  select(-c(Feature.ID, Domain, Phylum, Class, Order, Family, Genus, Species, Consensus)) %>%
  summarize_all(funs(sum)) %>%
  mutate(sum = rowSums(.))

#Transform to long format
total.s.asv.long <- total.s.asv %>%
  pivot_longer(cols = -sum, names_to = "sample_id", values_to = "ASV_Count")

#write.csv(total.g.asv.long, "ASV_per_sample-id_2000_gutsponge.csv")

#Prepare dataframe
asv.s.table.f<- asv.taxon.s.f %>%
  select(-c(Domain, Phylum, Class, Order, Family, Genus, Species, Consensus))

asv.s.table.f <- column_to_rownames(asv.s.table.f, var = colnames(asv.s.table.f)[1])

asv.s.table.df <- as.data.frame(asv.s.table.f)

t.asv.s.table<- t(asv.s.table.df)

##Calculate rarefaction curve
rare_curve <- rarecurve(t.asv.s.table, step = 100)

sample_names<- row.names(t.asv.s.table)

rare <- lapply(rare_curve, function(x){
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names

rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

#missing_samples <- setdiff(meta.s.f$sample_id, unique(rare$sample))
#print(missing_samples)

rare<- map_dfr(meta.s.f$sample_id, function(x){ #loop over samples
  z <- rare[rare$sample == x,] #subset rare according to sample 
  Type <- meta.s.f$type[meta.s.f$sample_id == x] #subset groupings according to sample, if more than one grouping repeat for all
  z <- data.frame(z, Type) #make a new data frame with the subsets
  return(z)
})

#Plotting rarefaction curves - merged

rare_curve<- ggplot(data = rare)+
  geom_line(aes(x = raw.read, y = ASV, group = sample, color = Type))+
  geom_text(data = rare %>% #here we need coordinates of the labels
              group_by(sample) %>% #first group by samples
              summarise(max_ASV = max(ASV), #find max OTU
                        max_raw = max(raw.read)), #find max raw read
            aes(x = max_raw, y = max_ASV, label = NA), check_overlap = T, hjust = 0)+
  scale_colour_manual(values = c("Risor skin" = "#CC79A7", "sponge" = "#D55E00"))+
  scale_x_continuous(labels =  scales::scientific_format())+
  theme_bw()

#ggsave("rare_curve_2000_skin_sponge_final.pdf", rare_curve, width = 8, height = 6, useDingbats = FALSE)
