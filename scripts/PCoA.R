library(phyloseq)
library(vegan)
library(tidyverse)
library(patchwork)
library(ggpubr)
source("https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R")

ps_ITS <- readRDS("~/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere/2023/out/ps_ITS.rds")
ps.rarefied.ITS = rarefy_even_depth(ps_ITS, rngseed=1, sample.size=0.9*min(sample_sums(ps_ITS)), replace=F)

# condition list that you want for your PCoA
flower_practice <- 'Flower'

leafs_practice <- 'Leaf'
leafs_practice_May <- c('Leaf', 'May')
leafs_practice_July <- c('Leaf', 'July')
leafs_practice_August <- c('Leaf', 'August')

leafs_HC_SP_practice <- c('Spartan', 'Honeycrisp')
leafs_HC_SP_practice_May <- c('Spartan', 'Honeycrisp', 'Leaf', 'May')
leafs_HC_SP_practice_July <- c('Spartan', 'Honeycrisp', 'Leaf', 'July')
leafs_HC_SP_practice_August <- c('Spartan', 'Honeycrisp', 'Leaf', 'August')

conditions <- list(flower_practice,leafs_practice, leafs_practice_May,
                   leafs_practice_July, leafs_practice_August, leafs_HC_SP_practice,
                   leafs_HC_SP_practice_May, leafs_HC_SP_practice_July, leafs_HC_SP_practice_August)

# function to subset ps object in function of condition (hard coded)
subset_ps_by_condition <- function(ps, condition) {
  
  # transform in sample_data
  samdf <- as(phyloseq::sample_data(ps), "data.frame")

  # initialize filter using sample names
    # we will create a vector of the parameters we want
  keep_samples <- rownames(samdf)
  
  # apply cultivar filter (OR logic)
  if (any(c("Spartan", "Honeycrisp") %in% condition)) {
    cultivars <- condition[condition %in% c("Spartan", "Honeycrisp")]
    keep_samples <- keep_samples[samdf[keep_samples, "cultivar"] %in% cultivars]
  }
  
  # apply type filter
  if ("Leaf" %in% condition) {
    keep_samples <- keep_samples[samdf[keep_samples, "type"] == "Leaf"]
  }
  if ("Flower" %in% condition) {
    keep_samples <- keep_samples[samdf[keep_samples, "type"] == "Flower"]
  }
  
  # apply time filter
  if ("May" %in% condition) {
    keep_samples <- keep_samples[samdf[keep_samples, "time"] == "May"]
  }
  if ("July" %in% condition) {
    keep_samples <- keep_samples[samdf[keep_samples, "time"] == "July"]
  }
  if ("August" %in% condition) {
    keep_samples <- keep_samples[samdf[keep_samples, "time"] == "August"]
  }
  
  # safely subset and prune with new list
  ps_sub <- phyloseq::prune_samples(keep_samples, ps) %>% 
    phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0, .)
  
  return(ps_sub)
}


# apply with imap
  # results of subset
subset_ps_list <- purrr::imap(conditions, ~ subset_ps_by_condition(ps.rarefied.ITS, .x))

# name results based on conditions
names(subset_ps_list) <- sapply(conditions, paste, collapse = "_")

## function to prepare the data for plotting
PCoA.prep <- function(ps) {
  # extract counts table:
  counts_table <- ps %>% 
    vst_ps_to_mx() # variance-stabilizing transformation
  
  # compute Bray-Curtis dissimilarity
  PCOA <- capscale(counts_table ~ 1, distance = 'bray') # Calculate distances
  
  # first three Eigenvalues
  eig <- round(PCOA$CA$eig[1:3] / sum(PCOA$CA$eig), 2)
  
  # extract the firt two coordinates
  pcoa.df <- data.frame(sample_data(ps))
  pcoa.df$PCo1 <- scores(PCOA)$sites[, 1]
  pcoa.df$PCo2 <- scores(PCOA)$sites[, 2]
  
  # add eigenvalues as new columns (same for all rows)
  pcoa.df$Eig1 <- eig[1]
  pcoa.df$Eig2 <- eig[2]
  pcoa.df$Eig3 <- eig[3]
  
  return(pcoa.df)
}

# apply to list of phyloseq object
prep.data.PCoA <- imap(subset_ps_list, function(subset_ps, .y) {PCoA.prep(subset_ps) })



# function to plot PCoA in function of PCoA.df 
plotting.PCoA <- function(PCoA.df) {
  # define site_shapes
  site_shapes <- c(
    "ASB" = 21,
    "COM" = 22,
    "MIB" = 23,
    "MIC" = 24,
    "PMB" = 25, 
    "VBS" = 8
  )
  
  PCoA.df %>% 
    ggplot(aes(x = PCo1, y = PCo2, colour = practice)) +
    stat_ellipse(aes(fill = practice), 
                 level = 0.95, 
                 geom = 'polygon', 
                 alpha = 0.2) +
    geom_point(aes( shape = site, fill = practice), color = "black", size = 3) +
    scale_shape_manual(values = site_shapes) +
    labs(x = paste("PCo1:", 100*PCoA.df$Eig1, "%"),
         y = paste("PCo2:", 100*PCoA.df$Eig2, "%"),
         colour = "Practice",
         fill = "Practice",
         shape = "Site") +
    theme_light()
}

plot.PCoA <- imap(prep.data.PCoA, function(PCoA.df.subset, .y) {plotting.PCoA(PCoA.df.subset)})

# exports all plot individually 
for (name in names(plot.PCoA)) {
  ggsave(filename = paste0("./final/", name, ".png"), plot = plot.PCoA[[name]], units = "px", width = 1000, height = 750, dpi = 200)
}

# delete legend from Flower because it is different and it will cause problem in the patchwork
plot.PCoA$Flower <- plot.PCoA$Flower + 
  guides(color = "none", fill = "none", shape = "none")

legend <- get_legend(plot.PCoA$Leaf) %>% as_ggplot()


# combine plots (1 legend of Site)
combined_plot <- 
  (plot.PCoA$Flower | plot.PCoA$Leaf | plot.PCoA$Spartan_Honeycrisp) / 
  (plot.PCoA$Leaf_May | plot.PCoA$Leaf_July | plot.PCoA$Leaf_August) /
  (plot.PCoA$Spartan_Honeycrisp_Leaf_May | plot.PCoA$Spartan_Honeycrisp_Leaf_July | plot.PCoA$Spartan_Honeycrisp_Leaf_August) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = NULL,
    tag_suffix = NULL,
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 0)
    )
  ) & theme(legend.text = element_text(size = 13),
                   legend.title = element_text(size = 16))

# export all plots as combined in a patchwork
ggsave("./final/PCoA_combined.png", combined_plot, units = "px", width = 3000, height = 1600, dpi = 200)

