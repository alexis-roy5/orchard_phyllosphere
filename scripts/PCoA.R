library(vegan)
library(tidyverse)

ps_ITS <- readRDS("~/Documents/1_Université/Stages/Labo_ILL/orchard_phyllosphere/2023/out/ps_ITS.rds")
ps.rarefied.ITS = rarefy_even_depth(ps_ITS, rngseed=1, sample.size=0.9*min(sample_sums(ps_ITS)), replace=F)



# Calculer les distances (ou dissimilarités dans ce cas-ci)
d <- vegdist(ps.rarefied.ITS@otu_table, method = "bray")
View(d) # une matrice de distances entre paires d'échantillons

pcoa_res <- cmdscale(d, k = 3, eig = TRUE)
View(pcoa_res$point) #coordonnées des trois premiers axes (t'aurais aussi pu mettre k = length(d)-1 pour avoir tous les axes mais on s'en fout un peu, toutes les eigenvalues vont être disponibles pareil)

pcoa_coords <- pcoa_res$points %>%
  as.data.frame() %>%
  rename_with(~ paste0("PCo", seq_along(.))) %>% # nommer les colonnes PCo1, PCo2, ...
  select(PCo1, PCo2, PCo3) %>% 
  bind_cols(ps.rarefied.ITS@sam_data %>% as.matrix())

# Plot PCoA !  tu peux changer la variable type par ce que tu veux visualiser.
# PCoAs. Flowers vs leaves
all.may <- ggplot(pcoa_coords %>% filter(time == "May"), aes(x = PCo1, y = PCo2, colour = type)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2) +
  theme_light() ; all.may

HC.vs.SP.may <- ggplot(pcoa_coords %>% filter(time == "May" & cultivar %in% c("Spartan", "Honeycrisp")), aes(x = PCo1, y = PCo2, colour = cultivar)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = cultivar)) +
  theme_light() ; HC.vs.SP.may 

# PCoAs. Only leaves. 
leafs.time <- ggplot(pcoa_coords %>% filter(type == "Leaf"), aes(x = PCo1, y = PCo2, colour = time)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = time)) +
  theme_light() ; leafs.time

leafs.site <- ggplot(pcoa_coords %>% filter(type == "Leaf"), aes(x = PCo1, y = PCo2, colour = site)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = site)) +
  theme_light() ; leafs.site

# PCoAs. HC vs SP. 
HC.vs.SP.time <- ggplot(pcoa_coords %>% filter(cultivar == "Honeycrisp" | cultivar == "Spartan"), aes(x = PCo1, y = PCo2, colour = time)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = time)) +
  theme_light() ; HC.vs.SP.time

HC.vs.SP.site <- ggplot(pcoa_coords %>% filter(cultivar == "Honeycrisp" | cultivar == "Spartan"), aes(x = PCo1, y = PCo2, colour = site)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = site)) +
  theme_light() ; HC.vs.SP.site


# PCoAs. Conventional vs Organic
all.org.vs.conv <- ggplot(pcoa_coords, aes(x = PCo1, y = PCo2, colour = practice, shape = time)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = practice)) +
  theme_light() ; all.org.vs.conv 

HC.vs.SP.org.vs.conv <- ggplot(pcoa_coords %>% filter(cultivar == "Honeycrisp" | cultivar == "Spartan"), aes(x = PCo1, y = PCo2, colour = practice)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = practice)) +
  theme_light() ; HC.vs.SP.org.vs.conv

B1.vs.B2.org.vs.conv <- ggplot(pcoa_coords %>% filter(code == "B1" | code == "B2"), aes(x = PCo1, y = PCo2, colour = practice)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = practice)) +
  theme_light() ; B1.vs.B2.org.vs.conv

D1.vs.D2.org.vs.conv <- ggplot(pcoa_coords %>% filter(code == "D1" | code == "D2"), aes(x = PCo1, y = PCo2, colour = practice)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = practice)) +
  theme_light() ; D1.vs.D2.org.vs.conv

B1.B2.D1.D2.org.vs.conv <- ggplot(pcoa_coords %>% filter(code == "B1" | code == "B2" | code == "D1" | code == "D2"), 
                                  aes(x = PCo1, y = PCo2, colour = practice)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = practice)) +
  theme_light() ; B1.B2.D1.D2.org.vs.conv

# PCoAs. Susceptibility. 
A.B1.C.susceptibility <- ggplot(pcoa_coords %>% filter(code == "A" | code == "B1" | code == "C"), 
                                  aes(x = PCo1, y = PCo2, colour = susceptibility)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = susceptibility)) +
  theme_light() ; A.B1.C.susceptibility

flowers.susceptibility <- ggplot(pcoa_coords %>% filter(type == "Flower"), 
                                aes(x = PCo1, y = PCo2, colour = susceptibility)) +
  geom_point() +
  stat_ellipse(level = 0.95, geom = 'polygon', alpha = 0.2, aes(fill = susceptibility)) +
  theme_light() ; flowers.susceptibility

