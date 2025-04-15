library(pacman)
p_load(phyloseq, vegan, tidyverse, patchwork, ggpubr, DESeq2)
source("https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/rarefy_even_depth2.R")
source("https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R")


ps_ITS <- readRDS("~/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere/2023/out/ps_ITS.rds")
ps.rarefied.ITS = rarefy_even_depth(ps_ITS, rngseed=1, sample.size=0.9*min(sample_sums(ps_ITS)), replace=F)
# raréfaction à 2257 séquences
# rngseed: le nombre donner à set_seed pour un processus aléatoire qui est pas un processus aléatoire. 

#################
### perMANOVA ###
#################
# Pre-calculate the sample size to avoid repeated calculations
target_size <- 0.9 * min(phyloseq::sample_sums(ps_ITS))

# Initialize list with exact length for better memory management
list_rarefaction_seed_1_100 <- vector("list", 100)

# Use lapply for more efficient iteration
list_rarefaction_seed_1_100 <- lapply(1:100, function(seed) {
  rarefy_even_depth2(
    physeq = ps_ITS,
    rngseed = seed,
    sample.size = target_size,
    replace = FALSE,
    ncores = 7,
    verbose = TRUE
  )
})

list_rarefaction_seed_1_101 <- append(list_rarefaction_seed_1_100, ps_ITS)

# deseq list
deseq_physeq_list <- lapply(list_rarefaction_seed_1_101, function(ps) {
  phyloseq_to_deseq2(ps, ~orchard) # parameter to change...
})

data_transformation_permanova <- function(deseq_object) {
  
  # create function for geometric mean (taken from online website of phyloseq)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
  # estimate factor size
  deseq = estimateSizeFactors(deseq_object, geoMeans = apply(counts(deseq_object), 1, gm_mean))
  
  vst <- DESeq2::varianceStabilizingTransformation(deseq, blind=T)
  vst.mat <- SummarizedExperiment::assay(vst) # extract transformed asv table
  t.vst.mat<-t(vst.mat)
  t.vst.mat[which(t.vst.mat<0)]<-0
  
  t.pcoa.vst <- vegdist(t.vst.mat, distance="bray")
  
  return(t.pcoa.vst)
}

# prepared data
data_permanova_list <- lapply(deseq_physeq_list, function(deseq_object) {
  data_transformation_permanova(deseq_object)
})

# sample data
sample.df <- as(sample_data(list_rarefaction_seed_1_100[[1]]), "data.frame") %>% 
  mutate(lost_percentage = (sample.df$seqDepth - 2257) / sample.df$seqDepth *100)

mean(sample.df$lost_percentage)
sd(sample.df$lost_percentage)

# number of permutations
perm <- how(nperm = 999)

setBlocks(perm) <- with(sample.df, TreeID)

permanova.list <- lapply(data_permanova_list, function(prepared_data) {
  vegan::adonis2(prepared_data ~ practice + time + cultivar + type, 
                 data = sample.df,  
                 permutations = perm,
                 parallel = 6,
                 by = "margin")
})

perm.out <- imap(permanova.list, function(permanova_results, .y) {
  # initial df
  df <- data.frame(
    Variable = rownames(permanova_results),
    R2 = permanova_results$R2, 
    Seed = .y) %>% 
    filter(!Variable %in% c("Total", "Residual")) %>% 
    mutate(Variable = as.factor(Variable))
  return(df)
}) %>% list_rbind()

# vérifier si la raréfaction a un effet sur la beta-diversity (perMANOVA) en comparant 100 objet raréfié avec un objet ps non-raréfié
perm.out %>% 
  ggplot(aes(x = Variable, y = R2)) +
  geom_jitter(aes(color = ifelse(Seed == 101, "red", "black")), width = 0.2, height = 0, size = 2) +
  guides(color="none")
# conclusion: time et practice, sont différents!
# R2: proportion de variance expliqués

# est-ce que ces variables (time et practice) qui sont différentes dans les permanova sont confondues avec la profondeur de séquençage
# la seule différence entre raréfction et non-raréfié = mettre la profondeur de séquençage au même niveau
# est-ce la profondeur de séquençage confound notre modèle
sample_data(ps_ITS) %>% ggplot(aes(x=time, y=seqDepth)) + geom_boxplot()
sample_data(ps_ITS) %>% ggplot(aes(x=practice, y=seqDepth)) + geom_boxplot()
sample_data(ps_ITS) %>% ggplot(aes(x=type, y=seqDepth)) + geom_boxplot()

# une seule perMANOVA avec le seqdepth
perm.non.rarefied <- adonis2(data_permanova_list[[101]] ~ practice + time + cultivar + type + seqDepth, 
           data = sample.df,  
           permutations = perm,
           parallel = 6,
           by = "margin")


################ 
##### PCoA #####
################

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
subset_ps_list <- purrr::imap(conditions, ~ subset_ps_by_condition(ps_ITS, .x))

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
    "VBS" = 4
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

# Highlight some sites!
plotting.PCoA.highlight <- function(PCoA.df, site_1, site_2) {
  site_shapes <- c(
    "ASB" = 21,
    "COM" = 22,
    "MIB" = 23,
    "MIC" = 24,
    "PMB" = 25, 
    "VBS" = 4
  )
  
  # Create a new column to identify the two highlighted sites
  PCoA.df <- PCoA.df %>%
    mutate(highlight = ifelse(site %in% c(site_1, site_2), site, "Other"))
  
  PCoA.df %>% 
    ggplot(aes(x = PCo1, y = PCo2)) +
    # Colored ellipses ONLY for highlighted sites
    stat_ellipse(
      data = subset(PCoA.df, highlight != "Other"),
      aes(fill = highlight, color = highlight), 
      level = 0.95, 
      geom = 'polygon', 
      alpha = 0.2
    ) +
    # All points (grey for non-highlighted, colored for highlighted)
    geom_point(
      aes(shape = site, fill = highlight),
      color = "black",  # outline color
      size = 3
    ) +
    scale_shape_manual(values = site_shapes) +
    scale_fill_manual(
      values = c(
        setNames("slateblue1", site_1),  # orange for site_1
        setNames("#0072B2", site_2),  # blue for site_2
        "Other" = "grey70"            # grey for others
      )
    ) +
    scale_color_manual(
      values = c(
        setNames("slateblue1", site_1),  # orange for site_1
        setNames("#0072B2", site_2)   # blue for site_2
      )
    ) +
    labs(
      x = paste("PCo1:", 100*PCoA.df$Eig1[1], "%"),
      y = paste("PCo2:", 100*PCoA.df$Eig2[1], "%"),
      shape = "Site",
      color = "Highlight"
    ) +
    theme_light() +
    guides(
      fill = "none"   # Hide color legend
    )
}

# plot list: conditions <- flower_practice,leafs_practice, leafs_practice_May,
# leafs_practice_July, leafs_practice_August, leafs_HC_SP_practice,
# leafs_HC_SP_practice_May, leafs_HC_SP_practice_July, leafs_HC_SP_practice_August)

plot.PCoA <- imap(prep.data.PCoA, function(PCoA.df.subset, .y) {plotting.PCoA(PCoA.df.subset)}); plot.PCoA

# subset of flowers with COM and PMB highlighted
plot.PCoA.GREY.COM.PMB <- imap(prep.data.PCoA["Flower"], function(PCoA.df.subset, .y) {plotting.PCoA.highlight(PCoA.df.subset, "COM", "PMB")}); plot.PCoA.GREY.COM.PMB

# subset of leafs at different time points with MIC and MIB highlighted
plot.PCoA.GREY.MIB.MIC <- imap(prep.data.PCoA[c("Leaf_May", "Leaf_July", "Leaf_August")], function(PCoA.df.subset, .y) {plotting.PCoA.highlight(PCoA.df.subset, "MIB", "MIC")}); plot.PCoA.GREY.MIB.MIC

# exports all plot individually 
#for (name in names(plot.PCoA)) {
  #ggsave(filename = paste0("./final/", name, ".png"), plot = plot.PCoA[[name]], units = "px", width = 1000, height = 750, dpi = 200)
#}

# delete legend from Flower because it is different and it will cause problem in the patchwork
#plot.PCoA$Flower <- plot.PCoA$Flower + 
  #guides(color = "none", fill = "none", shape = "none")

###################
#### Figure 1 #####
###################

figure1 <- 
  (plot.PCoA$Flower | plot.PCoA.GREY.COM.PMB$Flower) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = NULL,
    tag_suffix = NULL
  ) & theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 0),
             legend.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 13),
            legend.title = element_text(size = 16))

ggsave("./final/PCoA/FIG1.png", units = "px", width = 1650, height = 750, dpi = 225) 

###################
#### Figure 2 #####
###################
# leafs plot across different time points AND leafs plot across different time points - MIB and MIC highlighted

figure2 <- 
  (plot.PCoA$Leaf_May | plot.PCoA$Leaf_July | plot.PCoA$Leaf_August) /
  (plot.PCoA.GREY.MIB.MIC$Leaf_May | plot.PCoA.GREY.MIB.MIC$Leaf_July | plot.PCoA.GREY.MIB.MIC$Leaf_August) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = NULL,
    tag_suffix = NULL
  ) & theme(plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 0),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 13),
                   legend.title = element_text(size = 16))

ggsave("./final/PCoA/FIG2.png", units = "px", width = 1650, height = 1000, dpi = 200) 

#####################
## SUPPLEMENTARIES ##
#####################
## patchwork all leafs only Honeycrisp and Spartan and 3 different times points only Honeycrisp and Spartan

supplementary_combined_plot <- 
  (plot.PCoA$Spartan_Honeycrisp | plot.PCoA$Spartan_Honeycrisp_Leaf_May)  / 
  (plot.PCoA$Spartan_Honeycrisp_Leaf_July | plot.PCoA$Spartan_Honeycrisp_Leaf_August) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    tag_levels = 'A',
    tag_prefix = NULL,
    tag_suffix = NULL
  ) & theme(plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 0),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 13),
            legend.title = element_text(size = 16))

ggsave("./final/PCoA/supplementary_HC_SP.png", units = "px", width = 1250, height = 1000, dpi = 200) 

# export all plots as combined in a patchwork
#ggsave("./final/PCoA_combined.png", combined_plot, units = "px", width = 3000, height = 1600, dpi = 200)

