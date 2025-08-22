library(pacman)
p_load(phyloseq, tidyverse, ANCOMBC, RColorBrewer, patchwork, magrittr, 
       update = FALSE)

source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/tax_glom2.R')
source('https://raw.githubusercontent.com/jorondo1/misc_scripts/refs/heads/main/community_functions.R')


# DAA TEST  ---------------------------------------------------------------

ps <- read_rds('2023/out/ps_ITS.rds')
taxLvl <- 'Genus'

ps.flower <- subset_samples(
  ps, type == 'Flower'
)  %>% tax_glom2(taxrank = taxLvl) # tax_glom2 will remove taxa absent from all samples

ps.leaf <- subset_samples(
  ps, type == 'Leaf'
) %>% tax_glom2(taxrank = taxLvl)

# # ANCOM model
# Execute this once then save it below
# ancom_flower.out <- ancombc2(
#   data = ps.flower,
#   prv_cut = 0.10,
#   fix_formula="practice + cultivar",
#   #rand_formula = "(1 | TreeID)", # Not possible with this subset
#   group = "practice", # specify group if >=3 groups exist, allows structural zero detection
#   struc_zero = TRUE,
#   alpha = 0.01,
#   pairwise = TRUE,
#   verbose = TRUE,
#   n_cl = 8 # cores for parallel computing
# )
# 
# ancom_leaf.out <- ancombc2(
#   data = ps.leaf,
#   prv_cut = 0.10,
#   fix_formula="practice + cultivar",
#   #  rand_formula = "(1 | TreeID)",
#   group = "practice", # specify group if >=3 groups exist, allows structural zero detection
#   struc_zero = TRUE,
#   alpha = 0.01,
#   pairwise = TRUE,
#   verbose = TRUE,
#   n_cl = 8 # cores for parallel computing
# )

# write_rds(ancom_leaf.out, '2023/data/ITS/DAA/ancom_Genus_leaf.rds')
# write_rds(ancom_flower.out, '2023/data/ITS/DAA/ancom_Genus_flower.rds')


# PROCESSING FUNCTIONS ----------------------------------------------------

# Filter table
filter_ancom <- function(ancom_out,
                         signif_threshold, 
                         num_tests) { 
  # num_tests is a manual Holm correction, because ancom will only have corrected 
  # within each test, not accounting for the fact we did two separate tests.
  
  tibble(ancom_out$res) %>%
    dplyr::select(-starts_with('W_'), -starts_with('diff_'),
                  -starts_with('q_') # we ignore the built-in p-val corrections
                  ) %>% 
    mutate(q_practiceOrganic = p.adjust(p_practiceOrganic, 'holm', num_tests)) %>% 
    # Remove taxon for which no comparison passed the ss 
    filter(!passed_ss_practiceOrganic == FALSE) %>%
    # Then only keep taxa where at least one has q < 0.01
    filter(!q_practiceOrganic > signif_threshold) %>% 
    rowwise() %>% 
    # Remove taxon where no comparison is both <= signif_threshold & passed_ss TRUE
    filter(any(
      q_practiceOrganic <= signif_threshold
      & passed_ss_practiceOrganic == TRUE
    )) %>% ungroup()
}

# Format for plotting
ancom_plot_data <- function(ancom_filtered, ps) { 
  
  plot_ancom.df <- ancom_filtered %>% 
    dplyr::select(taxon, all_of(contains('practice'))) %>% 
    # # long format 
    pivot_longer(cols = -taxon, 
                 names_to = c(".value", "Group"), 
                 names_pattern = "(lfc|se|q|passed_ss)_(.+)", 
                 values_drop_na = TRUE) %>% 
    # lfc become 0 when q > threshold for plotting purposes
    mutate(
      Group = factor(
        case_when(lfc > 0 ~ 'Organic',
                  TRUE ~ 'Conventional')
      )
    ) %>% 
    left_join(ps %>% # identifier \ species association table
                tax_table() %>% data.frame() %>% 
                select(all_of(taxLvl), Phylum, Class, Order) %>% tibble(),
              join_by(taxon == !!sym(taxLvl))) %>% 
    filter(lfc!=0) %>% 
    mutate(taxon = str_replace_all(taxon, "_", " ")) # for *_gen_Incertae_sedis
  
  # Arrange taxa by descendnig lfc for waterfall
  taxon_order <- plot_ancom.df %>% 
    arrange(desc(lfc)) %>% 
    pull(taxon)
  
  plot_ancom.df %<>%
    mutate(taxon = factor(taxon, levels = taxon_order))
}

# Main plotting function
plot_ancom <- function(plot_data, taxLvl_palette, tile_rank) {
  
  plot_data <- droplevels(plot_data)
  
  # create alternating Grey-white background for plot
  bg_waterfall_data <- plot_data %>%
    distinct(taxon) %>% # Get unique taxa
    mutate(
      taxon_index = as.numeric(as.factor(taxon)),
      bg_color_id = (taxon_index %% 2) # 0 for pale grey, 1 for white
    )
  
  bg_colors <- c("0" = "grey90", "1" = "white")
  
  # If too many tested taxon levels for the palette, expand it
  tile_rank_N <- length(unique(plot_data[[tile_rank]]))
  max_colours <- RColorBrewer::brewer.pal.info[taxLvl_palette,'maxcolors']
  taxLvl_palette <- colorRampPalette(brewer.pal(max_colours,taxLvl_palette))(tile_rank_N)

  # Coloured tiles by higher taxonomic rank
  tiles <- plot_data %>% 
    ggplot(aes(x = '1', y = taxon, fill = !!sym(tile_rank))) +
    geom_tile(width = 2) + theme_minimal()+
    theme(axis.text.x = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_text(hjust = 1)) +
    scale_fill_manual(values = taxLvl_palette)
  
  # Main plot
  main_plot <- plot_data %>% 
    ggplot(aes(x = lfc, y = taxon, fill = Group)) +
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
    # Enrichment colour code :
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
}

# EXECUTION ---------------------------------------------------------------

# Filter ancom out
ancom_leaf.out <- read_rds('2023/data/ITS/DAA/ancom_Genus_leaf.rds')
ancom_flower.out <- read_rds('2023/data/ITS/DAA/ancom_Genus_flower.rds')

# Account for both series of tests in pval correction: 
num_tests <- nrow(ancom_flower.out$res)+nrow(ancom_leaf.out$res)

leaf_filtered <- filter_ancom(
  ancom_out = ancom_leaf.out, 
  signif_threshold = 0.01, num_tests = num_tests)

flower_filtered <- filter_ancom(
  ancom_out = ancom_flower.out, 
  signif_threshold = 0.01, num_tests = num_tests)

# Prep plot data
leaf_plot_data <- ancom_plot_data(leaf_filtered, ps.leaf)
flower_plot_data <- ancom_plot_data(flower_filtered, ps.flower)

# PLOTS : 

leaf_plot_data %>% 
  filter(lfc>=1 | lfc <= -1) %>% 
  plot_ancom(taxLvl_palette = 'Set3',
           tile_rank = 'Phylum')
ggsave('2023/out/DAA_leaf_genus_phylum.pdf', 
       bg = 'white', width = 1600, height = 2000, 
       units = 'px', dpi = 220)
theme(panel.grid = )
plot_ancom(flower_plot_data, 
           taxLvl_palette = 'Set3',
           tile_rank = 'Phylum')
ggsave('2023/out/DAA_flower_genus_phylum.pdf', 
       bg = 'white', width = 1600, height = 2000, 
       units = 'px', dpi = 220)

leaf_plot_data %>% 
  filter(lfc>=1 | lfc <= -1) %>% 
  plot_ancom(taxLvl_palette = 'Set1',
           tile_rank = 'Class')
ggsave('2023/out/DAA_leaf_genus_class.pdf', 
       bg = 'white', width = 1600, height = 2000, 
       units = 'px', dpi = 220)

plot_ancom(flower_plot_data, 
           taxLvl_palette = 'Set1',
           tile_rank = 'Class')
ggsave('2023/out/DAA_flower_genus_class.pdf', 
       bg = 'white', width = 1600, height = 2000, 
       units = 'px', dpi = 220)

# Raw data table as html
flower_plot_data %>% 
  select(taxon, lfc, se, q) %>% 
  mutate(across(where(is.numeric), ~round(.x, 4)),
         q = case_when(q < 0.0001 ~ '< 0.0001',
                       TRUE ~ as.character(q))) %>% 
  kable("html",
        caption = '') %>%
  kable_styling(full_width = FALSE) %>% 
  save_kable(file = '2023/out/DAA_table_flower.html')
  
leaf_plot_data %>% 
  select(taxon, lfc, se, q) %>% 
  mutate(across(where(is.numeric), ~round(.x, 4)),
         q = case_when(q < 0.0001 ~ '< 0.0001',
                     TRUE ~ as.character(q))) %>% 
  kable("html",
        caption = '') %>%
  kable_styling(full_width = FALSE) %>% 
  save_kable(file = '2023/out/DAA_table_leaf.html')

