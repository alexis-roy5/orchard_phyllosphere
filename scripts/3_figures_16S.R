#######################
## FIGURES & METRICS ##
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
library(RColorBrewer)
library(ggdist)


# setwd("/Users/alexisroy/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere")

ps_16S <- readRDS("./2023/out/ps_16S.rds")

##############################
##### DIVERSITY METRICS #####
#############################

## FUNCTIONS ##

# transform the phyloseq object in data.frame
sam_data_as_tibble <- function(ps) {
  require('dplyr', 'phyloseq')
  data.frame(sample_data(ps)) %>%
    rownames_to_column('Sample') %>% 
    tibble
}


# preparing the rarefaction curve
otu_matrix_t <- t(as(otu_table(ps_16S), "matrix"))


# rarefaction curve
rarecurve(otu_matrix_t, step = 50, label = FALSE, cex = 0.5, xlim = c(0, 10000))

# raréfaction 
ps.rarefied = rarefy_even_depth(ps_16S, rngseed=1, sample.size=0.9*min(sample_sums(ps_16S)), replace=F)
# 0,9: pour réchantilloner le petit sampling. 




# faire fondre l'objet phyloseq
phylo.df <- psflashmelt(ps_16S) # pour comparaison
phylo.rarefied.df <- psflashmelt(ps.rarefied)

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
    practice = dplyr::first(practice),
    type = dplyr::first(type),
    cultivar = dplyr::first(cultivar),
    RelativeAbundance = dplyr::first(RelativeAbundance),
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
  shapiro.test(residuals(model)) # shapiro.wilk
  
  
  # variances
  leveneTest(expShannonIndex ~ time, data = alpha.diversity.metrics.df)
  # if p-value < 0.05 -> variances are not equal
  
  # indépendance
  plot(model$fitted.values, residuals)
  abline(h = 0, col = "red")

# si les conditions du test ne sont pas remplies, on devrait faire un Wilcoxon-Mann-Whitney test -> dans le plot ...



# TIME TO PLOT
  # normal values
# Create the plot
  
nb.sample.months <- alpha.diversity.metrics.df %>%
  count(time) 
  
my_comparisons <- list( c("May", "July"), c("July", "August"),c("May", "August") ) # for stats

plot.1.general <- ggplot(R.S.S.index.normal, aes(x = time, y = DiversityIndex)) +
  geom_boxplot(alpha = 0.25, outlier.shape = NA, size = 0.5) +  
  geom_jitter(aes(color = cultivar), width = 0.15, size = 1.5) + 
  labs(
    x = NULL,
    y = "Alpha Diversity Index",
    subtitle = "Data from 2023",
    color = "Cultivar",
    fill = "Time"
  ) +
  facet_wrap(~DiversityMeasure, scales = "free_y", ncol = 3,
             labeller = as_labeller(c("Richness" = "Richness", 
                                      "ShannonIndex" = "Shannon Index", 
                                      "SimpsonIndex" = "Simpson Index"))) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", 
                     method = 'wilcox.test', p.adjust.method = "bonferroni") + 
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),  # Supprime le quadrillage principal
    panel.grid.minor = element_blank()  # Supprime le quadrillage mineur
  ) +
  scale_color_brewer(palette = "Dark2")

plot.1.general

ggsave("./2023/out/alpha_orchard.png", plot = plot.1.general, units = "px", width = 3000, height = 2000, dpi = 260)


### HILL CURVE ###
# Define the HILL function
HILL <- function(dataframe.abundance, q) {
  p <- dataframe.abundance$Richness / sum(dataframe.abundance$Richness)  # Proportional abundances
  if (q == 1) {
    # Shannon diversity (q = 1)
    return(exp(-sum(p * log(p))))
  } else {
    # General case for q != 1
    return((sum(p^q))^(1 / (1 - q)))
  }
}

# Example species abundance data
df.abundance <- phylo.rarefied.df %>% 
  filter(Abundance > 0) %>% 
  reframe(Sample = Sample,
          Richness = Abundance,
          ASV = OTU,
          Time = time)

# Create a sequence of q values from 0 to 2
q_values <- seq(0, 2, by = 0.01)


# function to calculate Hill numbers for a SINGLE sample
calculate_hill_for_sample <- function(sample_data, q_values) {
  results <- sapply(q_values, function(q) {
    HILL(sample_data, q)
  })
  return(data.frame(q = q_values, Hill = results))
}

# function for each sample using group_modify
hill_results <- df.abundance %>%
  group_by(Sample) %>% 
  group_modify(~ calculate_hill_for_sample(.x, q_values))

# adding time
hill_results_time <- df.abundance %>%
  select(Time, Sample) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  left_join(hill_results, ., by = "Sample")

plot.hill.curve <- ggplot(hill_results_time %>% 
                            group_by(Sample) %>% 
                            mutate(Time = factor(Time, levels = c("May", "July", "August"))), 
                          aes(x = q, y = Hill)) +
  geom_point(size = 0.1) +
  geom_smooth()+
  facet_wrap(~Time)

table(is.na(hill_results_time$Time))


## HILL NUMBERS - boxplot
plot.2.general.box <- ggplot(R.S.S.index.hill, aes(x = DiversityMeasure, y = DiversityIndex)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(aes(color = cultivar), width = 0.15, alpha = 0.5, size = 0.75) +
  labs(
    x = NULL, 
    y = "Alpha Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity (Hill Numbers) in Three Different Months (n = 349)", 
    subtitle = "Data from 2023", 
    color = "Cultivar"
  ) +
  facet_wrap(~ time) +
  scale_x_discrete(labels = c(
    "Richness" = "Richness",
    "expShannonIndex" = "Exponential Shannon Index",
    "InverseSimpsonIndex" = "Inverse Simpson Index"
  )) + 
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_brewer(palette = "Dark2")

plot.2.general.box

ggsave("./2023/out/alpha_orchard_hill.png", plot = plot.2.general.box, units = "px", width = 3000, height = 2000, dpi = 260)

### BIO/CONV ###

plot.3.bio <- ggplot(R.S.S.index.normal, aes(x = time, y = DiversityIndex, color = practice)) +
  geom_boxplot(outlier.shape = NA, size = 0.5) +  
  geom_jitter(width = 0.1, height = 0, size = 1.5) + 
  labs(
    x = NULL,
    y = "Alpha Diversity Index",
    subtitle = "Data from 2023",
    color = "Practice"
  ) +
  facet_wrap(~DiversityMeasure, scales = "free_y", ncol = 3,
             labeller = as_labeller(c("Richness" = "Richness", 
                                      "ShannonIndex" = "Shannon Index", 
                                      "SimpsonIndex" = "Simpson Index"))) +
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2")


plot.3.bio

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
    y = "Alpha Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
    subtitle = "Data from 2023",
    fill = "Time"
  ) +
  facet_grid(DiversityMeasure ~ orchard , scales = "free_y",
             labeller = labeller(
               DiversityMeasure = c(
                 "Richness" = "Richness", 
                 "ShannonIndex" = "Shannon Index",
                 "SimpsonIndex" = "Simpson Index"
               )
             )) +
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
    y = "Alpha Diversity Index",
    title = "Orchard Phyllosphere Alpha Diversity in Three Different Months", 
    subtitle = "Data from 2023",
    fill = "Time"
  ) +
  facet_grid(DiversityMeasure ~ orchard,  
             labeller = labeller(
               DiversityMeasure = c(
                 "Richness" = "Richness", 
                 "expShannonIndex" = "Exponential Shannon Index",
                 "InverseSimpsonIndex" = "Inverse Simpson Index"
               )
             )) +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    strip.text = element_text(size = 10)  # Adjust the size if needed
  ) 

plot.2.orchard


# Augmente la taille de base des textes
theme_para

theme_para <-  theme(
  plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),  # Titre plus grand
  axis.text = element_text(size = 14),  # Texte des axes
  axis.title = element_text(size = 16),  # Titres des axes
  legend.text = element_text(size = 14),  # Texte de la légende
  legend.title = element_text(size = 14, face = "bold"),  # Titre de la légende plus visible
  plot.caption = element_text(size = 14, hjust = 0.5, vjust = 3, face = "italic"),  # Caption plus lisible et centrée
  strip.text = element_text(size = 14, face = "bold"),  # Titres des facettes
  legend.position = "bottom",
  panel.grid.major = element_line(size = 0.8),  # Grilles plus épaisses
  panel.grid.minor = element_line(size = 0.5),  # Grilles secondaires
  axis.line = element_line(size = 1.2),  # Axes plus épais
  panel.border = element_rect(size = 1.5))  