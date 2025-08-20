library(pacman)
p_load(phyloseq, tidyverse, ANCOMBC, RColorBrewer, patchwork, magrittr, update = FALSE)
ps <- read_rds('2023/out/ps_ITS.rds')

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')

taxLvl <- 'Genus'
ps_gen <- tax_glom2(ps, taxrank = taxLvl)

# # ANCOM model
# Execute this once then save it below
ancom_out <- ancombc2(
  data = ps_gen,
  #tax_level= "taxLvl", # we do it ourselves above to retain the taxLvl level ids
  prv_cut = 0.10,
  fix_formula="practice + cultivar + type",
  rand_formula = "(1 | TreeID)",
  group = "practice", # specify group if >=3 groups exist, allows structural zero detection
  struc_zero = TRUE,
  alpha = 0.01,
  pairwise = TRUE,
  verbose = TRUE,
  n_cl = 8 # cores for parallel computing
)

# write_rds(ancom_out, '2023/data/ITS/DAA/ancom_Genus_withType.rds')
ancom_out <- read_rds('2023/data/ITS/DAA/ancom_Genus.rds')

# Parse output 
signif_threshold <- 0.01

# Filter table
ancom_filtered <- tibble(ancom_out$res) %>%
  dplyr::select(-starts_with('W_'), -starts_with('p_'), -starts_with('diff_')) %>% 
  # Remove taxon for which no comparison passed the ss 
  filter(!if_all(.cols = contains('passed_ss_'), 
                 .fns = ~ .x == FALSE)) %>%
  # Then only keep taxa where at least one has q < 0.01
  filter(!if_all(.cols = contains('q_practice'), 
                 .fns = ~ .x > signif_threshold)) %>% 
  rowwise() %>% 
  # Remove taxon where no comparison is both <= signif_threshold & passed_ss TRUE
  filter(any(
    c_across(contains('q_')) <= signif_threshold
    & c_across(contains('passed_ss_')) == TRUE
  )) %>% ungroup()


# Format for plotting
plot_ancom.df <- ancom_filtered %>% 
  select(taxon, all_of(contains('practice'))) %>% 
  # # long format 
  pivot_longer(cols = -taxon, 
               names_to = c(".value", "Group"), 
               names_pattern = "(lfc|se|q|passed_ss)_(.+)", 
               values_drop_na = TRUE) %>% 
  # lfc become 0 when q > threshold for plotting purposes
  mutate(
    across(c(lfc,se), ~ case_when(q > signif_threshold ~ 0, TRUE ~ .x)),
    Group = factor(
      case_when(lfc > 0 ~ 'Organic',
                TRUE ~ 'Conventional')
    ), q=q, .keep = 'unused'
  ) %>% 
  left_join(ps_gen %>% # identifier \ species association table
              tax_table() %>% data.frame() %>% 
              select(taxLvl, Phylum, Class, Order) %>% tibble(),
            join_by(taxon == !!sym(taxLvl))) %>% 
  filter(lfc!=0) 

# Arrange taxa by descendnig lfc for waterfall
taxon_order <- plot_ancom.df %>% 
  arrange(desc(lfc)) %>% 
  pull(taxon)

plot_ancom.df %<>%
  mutate(taxon = factor(taxon, levels = taxon_order))

# alternating Grey-white background for plot
bg_waterfall_data <- plot_ancom.df %>%
  distinct(taxon) %>% # Get unique taxa
  mutate(
    taxon_index = as.numeric(as.factor(taxon)),
    bg_color_id = (taxon_index %% 2) # 0 for pale grey, 1 for white
  )
bg_colors <- c("0" = "grey90", "1" = "white")

# Coloured tiles by higher taxonomic rank
tiles <- plot_ancom.df %>% 
  ggplot(aes(x = '1', y = taxon, fill = Phylum)) +
  geom_tile(width = 2) + theme_minimal()+
  theme(axis.text.x = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_brewer(palette = 'Set2')

# Main plot
main_plot <- plot_ancom.df %>% 
  ggplot(aes(x = lfc, y = taxon, fill = Group)) +
  #grey tile behind waterfall
  geom_rect(data = bg_waterfall_data, 
            aes(xmin = -Inf, xmax = Inf, 
                ymin = as.numeric(taxon) - 0.5, 
                ymax = as.numeric(taxon) + 0.5), 
            fill = bg_colors[factor(bg_waterfall_data$bg_color_id)], # DIRECTLY SET FILL
            alpha = 0.5, 
            inherit.aes = FALSE) +
  geom_col(width = 0.5) +
  geom_errorbar(aes(x = lfc,
                    xmin = lfc - se, 
                    xmax = lfc + se), 
                width = 0.2,
                linewidth = 0.3,
                position = position_dodge(width = 0.9), 
                color = "grey20") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.2) +
  
  scale_fill_manual(values = c('#F28E2B', '#499894')) +
  theme_minimal()  +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = 'Log2 fold changes in absolute abundances relative to organic',
       fill = 'Enriched in')

# Patchwork
tiles + main_plot +
  plot_layout(design = 'ABBBBBBBBBBBBBB', guides = 'collect') 

ggsave('2023/out/DAA_genus_phylum.pdf', 
       bg = 'white', width = 1600, height = 2000, 
       units = 'px', dpi = 220)

