#############################################
# calculer BC distance for tree individuals #
#############################################

library(tidyverse)
library(phyloseq)
library(vegan)

# rarefied ps_object
ps_ITS <- readRDS("~/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere/2023/out/ps_ITS.rds")
ps.rarefied.ITS = rarefy_even_depth(ps_ITS, rngseed=1, sample.size=0.9*min(sample_sums(ps_ITS)), replace=F)

df.meta <- data.frame(sample_data(ps.rarefied.ITS)) %>%
  rownames_to_column(var = "Sample") %>% 
  as_tibble() %>% 
  mutate(newID = str_extract(Sample, "2023-[1-3]-(...)-(..)-(..)")) %>% 
  select(-Sample, -Concentration, -sample, -seq, -replicate) %>% 
  unique

################################################  
df.otu <- as.data.frame(ps.rarefied.ITS@otu_table) %>%
  rownames_to_column(var = "Sample")

  
# stratégie, extraire les noms qui sont communs au réplicats pour les traiter comme 1 seul sample par la suite
df.otu$Sample_transform <- str_extract(df.otu$Sample, "2023-[1-3]-(...)-(..)-(..)") # extract the name from each ex: 
# for orchard_cultivar_type

# counts of réplicats in a table
sample_counts <- table(df.otu$Sample_transform)

# identify samples that appear more than once
non_single_replicats <- names(sample_counts[sample_counts > 1]) 
single_replicats <- names(sample_counts[sample_counts == 1]) # not interesting if =1 for BC, there will be no dissimilarity with himself

# subset the data to keep only the non-unique samples
df.otu$Sample_transform_not_unique <- df.otu$Sample_transform %in% non_single_replicats

# filter the dataframe to keep only the rows where the sample is not unique
# en effet les samples qui sont seuls auront une BC inexistante
df_non_unique <- df.otu[df.otu$Sample_transform_not_unique, ]

# split the data frame into smaller data frames - into a list of df
# each df is a sample with minimum 2 replicats
dfs_list <- split(df_non_unique, df_non_unique$Sample_transform)

# calculer BC distance en excluant les colonnes non-numériques
calculate_bray_curtis <- function(df) {
  numeric_df <- df %>% select(-Sample, -Sample_transform)  # exclude non-numeric columns for vegdist
  vegdist(numeric_df, method = "bray")
}

# apply the function to each data frame in the list
BC.rep <- lapply(dfs_list, calculate_bray_curtis)



# function for mean
calculate_mean <- function(df) {
  avg <- df %>% as.vector %>% mean   # calculate the mean
  return(avg)
}

# function for standard deviation
calculate_sd <- function(df) {  
  standard_deviation <- sd(df %>% as.vector())   # calculate the standard deviation of all column means
  return(standard_deviation)
}

mean.BC.sample <- map_dfr(BC.rep, function(dataframe) {calculate_mean(dataframe)}) %>% 
  pivot_longer(cols = everything(),
               names_to = 'newID',
               values_to = 'Mean.BC')

sd.BC.sample <- map_dfr(BC.rep, function(dataframe) calculate_sd(dataframe)) %>% 
  pivot_longer(cols = everything(),
               names_to = 'newID',
               values_to = 'sd.BC')

# big df
mean.sd <- left_join(mean.BC.sample, sd.BC.sample, by = "newID")
# NA = 2 sample


meta.mean.BC <- left_join(mean.sd, df.meta, by = "newID")

# plot
ggplot(meta.mean.BC) +
  geom_histogram(aes(x = Mean.BC), bins = 30)

ggplot(meta.mean.BC) +
  geom_boxplot(aes(x=time, y=Mean.BC))
