#######################
###### PHYLOSEQ #######
#######################

library(tidyverse)
library(phyloseq)
library(permute)
library(lattice)
library(vegan) # sanity check | diversity measure
library(cowplot)
library(ggpubr) # plot + stats
library(rstatix) # statistics
library(car) # levene test

setwd("/Users/alexisroy/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere")

ps_16S <- readRDS("~/2023/out/ps_16S.rds")

##############################
##### DIVERSITY METRICS #####
#############################

## FUNCTIONS ##

# transform the phyloseq object in data.frame
sam_data_as_tibble <- function(ps) {
  require('dplyr', 'phyloseq')
  sample_data(ps) %>%
    rownames_to_column('Sample') %>% 
    as.tibble
}


# preparing the rarefaction curve
otu_matrix_t <- t(as(otu_table(ps_16S), "matrix"))


# rarefaction curve
rarecurve(otu_matrix_t, step = 50, label = FALSE, cex = 0.5, xlim = c(0, 10000))

# raréfaction 
ps.rarefied = rarefy_even_depth(ps_16S, rngseed=1, sample.size=0.9*min(sample_sums(ps_16S)), replace=F)
# 0,9: pour réchantilloner le petit sampling. 




# faire fondre l'objet phyloseq
phylo.df <- psmelt(ps_16S) # pour comparaison
phylo.rarefied.df <- psmelt(ps.rarefied)

dim(phylo.df) # voir différences entre raréfié et non raréfié
dim(phylo.rarefied.df)


# calculate diversity metric in a table
alpha.diversity.metrics.df <- phylo.rarefied.df %>%
  filter(Abundance>0) %>% # filtrer les lignes (ASV) de 0
  group_by(Sample) %>% # grouper par sample
  mutate(Richness = n(), # unique ASV par sample
         RelativeAbundance = Abundance / sum(Abundance), 
         ShannonComponent = ifelse(RelativeAbundance > 0, 
                                   -RelativeAbundance * log(RelativeAbundance), # NOTE: fonctionnement du ifelse; ifelse(test, yes, no)
                                   0), # calculer pour chaque sample le (abondance relative * log(abondance relative)) terme
         SimpsonIndex = (sum(RelativeAbundance^2)), 
         ) %>%
  summarise(
    time = dplyr::first(time),
    orchard = dplyr::first(orchard),
    type = dplyr::first(type),
    cultivar = dplyr::first(cultivar),
    Richness = dplyr::first(Richness), # Preserver mes variable avec dplyr::first(column_name)
    ShannonIndex = sum(ShannonComponent),  
    expShannonIndex = exp(sum(ShannonIndex)), # formule pour expShannonIndex
    SimpsonIndex = dplyr::first(SimpsonIndex),
    InverseSimpsonIndex = 1 / sum(RelativeAbundance^2), # formule pour InverseSimpsonIndex
  ) %>% 
  ungroup() 

View(alpha.diversity.metrics.df) # ce tableau ne contient pas les contrôles

# sanity check
# extraire la table d'OTU de l'objet phyloseq raréfié
otu_table_matrix <- as(otu_table(ps.rarefied), "matrix")

# Calculate Shannon diversity
shannon_diversity <- diversity(otu_table_matrix, index = "shannon")

# Calculate Simpson diversity
simpson_diversity <- diversity(otu_table_matrix, index = "simpson") 

# Calculate Inverse Simpson diversity
invsimpson_diversity <- diversity(otu_table_matrix, index = "invsimpson")




##############################
########## PLOTS ############
#############################

# mettre en ordre inverse alphabétique la variable "time"
alpha.diversity.metrics.df <- alpha.diversity.metrics.df %>%
  mutate(time = factor(time, levels = rev(sort(unique(time)))))


# modifer le plot - BRILLANT
R.S.S.index.long <- alpha.diversity.metrics.df %>% 
  pivot_longer(cols = c(Richness, ShannonIndex, expShannonIndex, SimpsonIndex, InverseSimpsonIndex),
               names_to = "DiversityMeasure",
               values_to = "DiversityIndex") 
  

# filter for the normal diversity values
R.S.S.index.normal <- R.S.S.index.long %>%
  filter(DiversityMeasure %in% c("Richness", "ShannonIndex", "SimpsonIndex")) 

# filter for the Hill diversity values
R.S.S.index.hill <- R.S.S.index.long %>%
  filter(DiversityMeasure %in% c("Richness", "expShannonIndex", "InverseSimpsonIndex")) 


# mettre en ordre les mesures de diversité
R.S.S.index.normal$DiversityMeasure <- factor(R.S.S.index.normal$DiversityMeasure,
                                              levels = c("Richness", "ShannonIndex", "SimpsonIndex"))

R.S.S.index.hill$DiversityMeasure <- factor(R.S.S.index.hill$DiversityMeasure,
                                              levels = c("Richness", "expShannonIndex", "InverseSimpsonIndex"))



#####################
  # STATISTIQUES #
#####################

t_test <- t_test(alpha.diversity.metrics.df, expShannonIndex  ~ time)

# Vérification des conditions du test de t
# normalité des résidus; homoscédasticité des variances; indépendance des observation

  # normalité
  model <- lm(expShannonIndex ~ time, data = alpha.diversity.metrics.df) # ok?
  residuals <- residuals(model)
  hist(residuals) # visualisation
  
  # variances
  leveneTest(expShannonIndex ~ time, data = alpha.diversity.metrics.df)
  # if p-value < 0.05 -> variances are not equal
  
  # indépendance
  plot(model$fitted.values, residuals)
  abline(h = 0, col = "red")


# GENERAL
# Richness
kruskal.test(Richness ~ time, data = alpha.diversity.metrics.df) # is there a difference ? yes: p-value < 0.05. 
pairwise.wilcox.test(alpha.diversity.metrics.df$Richness, alpha.diversity.metrics.df$time,
                     p.adjust.method = "BH") # finding the difference between the groups

# Shannon index
kruskal.test(ShannonIndex ~ time, data = alpha.diversity.metrics.df) # is there a difference ? yes: p-value < 0.05. 
pairwise.wilcox.test(alpha.diversity.metrics.df$ShannonIndex, alpha.diversity.metrics.df$time,
                     p.adjust.method = "BH") # finding the difference between the groups. post-hoc test.

# Exponential Shannon index
kruskal.test(expShannonIndex ~ time, data = alpha.diversity.metrics.df) # is there a difference ? yes: p-value < 0.05. 
pairwise.wilcox.test(alpha.diversity.metrics.df$expShannonIndex, alpha.diversity.metrics.df$time,
                     p.adjust.method = "BH") # finding the difference between the groups. post-hoc test.

# Simpson index
kruskal.test(SimpsonIndex ~ time, data = alpha.diversity.metrics.df) # is there a difference ? yes: p-value < 0.05. 
pairwise.wilcox.test(alpha.diversity.metrics.df$SimpsonIndex, alpha.diversity.metrics.df$time,
                     p.adjust.method = "BH") # finding the difference between the groups. post-hoc test.

# Inverse Simpson Index
kruskal.test(InverseSimpsonIndex ~ time, data = alpha.diversity.metrics.df) # is there a difference ? yes: p-value < 0.05. 
pairwise.wilcox.test(alpha.diversity.metrics.df$InverseSimpsonIndex, alpha.diversity.metrics.df$time,
                     p.adjust.method = "BH") # finding the difference between the groups. post-hoc test.

# OVERALL
compare_means(c(Richness, ShannonIndex, SimpsonIndex ) ~ time,  data = alpha.diversity.metrics.df)
compare_means(c(Richness, expShannonIndex, InverseSimpsonIndex ) ~ time,  data = alpha.diversity.metrics.df)



# TIME TO PLOT
  # normal values
# Create the plot

my_comparisons <- list( c("May", "July"), c("July", "August"),c("May", "August") ) # for stats

plot.1.general <- ggplot(R.S.S.index.normal, aes(x = time, y = DiversityIndex)) +
  geom_boxplot(aes(fill = time), alpha = 0.25, outliers = FALSE) +
  geom_jitter(aes(shape = orchard, color = cultivar), width = 0.15, alpha = 0.5, size = 0.75) +
  labs(
    x = "Time", 
    y = "Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
    subtitle = "Data from 2023", 
    shape = "Orchard", 
    color = "Cultivar", 
    fill = "Time"
  ) +
  facet_wrap(~DiversityMeasure, scales = "free_y", ncol = 3,
             labeller = as_labeller(c("Richness" = "Richness", 
                                      "ShannonIndex" = "Shannon Index", 
                                      "SimpsonIndex" = "Simpson Index"))) +
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif") + 
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom")  
  
plot.1.general

# time to plot 
  # hill values
plot.2.general <- ggplot(R.S.S.index.hill, aes(x = time, y = DiversityIndex)) +
  geom_boxplot(aes(fill = time), alpha = 0.25, outliers = FALSE) +
  geom_jitter(aes(shape = orchard, color = cultivar), width = 0.15, alpha = 0.5, size = 0.75) +
  labs(
    x = "Time", 
    y = "Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
    subtitle = "Data from 2023", 
    shape = "Orchard", 
    color = "Cultivar", 
    fill = "Time"
  ) +
  facet_wrap(~ DiversityMeasure,
             labeller = as_labeller(c("Richness" = "Richness", 
                                      "expShannonIndex" = "Exponential Shannon Index", 
                                      "InverseSimpsonIndex" = "Inverse Simpson Index"))) +
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif") +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom")  

plot.2.general

##########################
## Orchard Specifically ##
##########################
# INUTILE lol
## NORMAL Values
plot.1.orchard <- ggplot(R.S.S.index.normal, aes(x = time, y = DiversityIndex)) +
  geom_boxplot(aes(fill = time), alpha = 0.25, outliers = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 0.75) +
  labs(
    x = "Time", 
    y = "Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
    subtitle = "Data from 2023",
    fill = "Time"
  ) +
  facet_wrap(~ orchard + DiversityMeasure, scales = "free_y",  ncol = 3, 
             labeller = labeller(
               DiversityMeasure = c(
                 "Richness" = "Richness", 
                 "ShannonIndex" = "Shannon Index",
                 "SimpsonIndex" = "Simpson Index"
               )
             )) +
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif") +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 10)  # Adjust the size if needed
  )  

plot.1.orchard


## HILL NUMBERS
plot.2.orchard <- ggplot(R.S.S.index.hill, aes(x = time, y = DiversityIndex)) +
  geom_boxplot(aes(fill = time), alpha = 0.25, outliers = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 0.75) +
  labs(
    x = "Time", 
    y = "Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
    subtitle = "Data from 2023",
    fill = "Time"
  ) +
  facet_wrap(~ orchard + DiversityMeasure, ncol = 3, 
             labeller = labeller(
               DiversityMeasure = c(
                 "Richness" = "Richness", 
                 "expShannonIndex" = "Exponential Shannon Index",
                 "InverseSimpsonIndex" = "Inverse Simpson Index"
               )
             )) +
  stat_compare_means(comparisons = my_comparisons, label =  "p.signif") +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 10)  # Adjust the size if needed
  ) 

plot.2.orchard


