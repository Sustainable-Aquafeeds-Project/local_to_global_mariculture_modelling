library(tidyverse)
library(cowplot)

# Data structure -----------------------------------------------------------------------------------------------------------------
# > str(CN_vals_forplot)
# tibble [99,000 × 17] (S3: tbl_df/tbl/data.frame)
#  $ farm_ID      : int [1:99000] 229 229 229 229 1854 1854 1854 1854 1855 1855 ...
#  $ feed_group   : Factor w/ 3 levels "MD","PD","AD": 1 1 1 1 1 1 1 1 1 1 ...
#  $ feed_type    : Factor w/ 3 levels "Minimum","Mean",..: 2 2 2 2 2 2 2 2 2 2 ...
#  $ measure      : Factor w/ 2 levels "excr","total_lost": 1 1 2 2 1 1 2 2 1 1 ...
#  $ mean         : num [1:99000] 9.88e+12 9.19e+13 1.27e+13 1.19e+14 1.89e+09 ...
#  $ mean_per_biom: num [1:99000] 9.88e+12 6.03e-02 8.35e-03 7.79e-02 1.89e+09 ...
#  $ nutrient     : Factor w/ 2 levels "C","N": 2 1 2 1 2 1 2 1 2 1 ...
#  $ iso3c        : chr [1:99000] "CHL" "CHL" "CHL" "CHL" ...
#  $ country      : Factor w/ 10 levels "Australia","Canada",..: 3 3 3 3 6 6 6 6 6 6 ...
#  $ species_group: chr [1:99000] "Salmonidae fish" "Salmonidae fish" "Salmonidae fish" "Salmonidae fish" ...
#  $ fao_code     : chr [1:99000] "87" "87" "87" "87" ...
#  $ hemisphere   : chr [1:99000] "S" "S" "S" "S" ...
#  $ model_name   : chr [1:99000] "atlantic_salmon" "atlantic_salmon" "atlantic_salmon" "atlantic_salmon" ...
#  $ geometry     :sfc_POINT of length 99000; first list element:  'XY' num [1:2] -73.7 -42.9
#  $ longitude    : num [1:99000] -73.7 -73.7 -73.7 -73.7 -24 ...
#  $ latitude     : num [1:99000] -42.9 -42.9 -42.9 -42.9 65.6 ...
#  $ region       : Factor w/ 7 levels "Australia","Chile",..: 2 2 2 2 3 3 3 3 3 3 ...


# Total lost only ----------------------------------------------------------------------------------------------------------------
CN_vals_total_lost <- CN_vals_forplot %>%
  filter(measure == "total_lost") %>% 
  select(feed_group, feed_type, mean_per_biom, nutrient, region)

nutrient_labels <- CN_vals_total_lost %>%
  distinct(region, nutrient) %>% 
  mutate(label = c(letters[1:14]))

p_CN_total_lost_1 <- CN_vals_total_lost %>% 
  mutate(mean_per_biom = set_units(mean_per_biom, "kg t_fish-1") %>% drop_units()) %>% 
  ggplot(aes(x = feed_type, y = mean_per_biom, fill = feed_group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = feed_pal) +
  facet_grid(
    rows = vars(nutrient), 
    cols = vars(region),
    scales = "free"
  ) +
  geom_text(
    data = nutrient_labels, 
    aes(x = -Inf, y = Inf, label = label),
    hjust = -1, vjust = 1.5,
    size = 6, fontface = "bold",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Feed digestibility",
    y = expression("Total loading (kg"~t[fish]^{-1}*")")
  ) +
  prettyplot_600() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank()
    )

# Different wrapping
regions <- levels(CN_vals_total_lost$region)
ls <- regions %>% 
  map(function(reg) {
    reg_df <- CN_vals_total_lost %>% 
      filter(region == reg) %>% mutate(regions = droplevels(region))
    p_CN_total_lost_1 %+% reg_df

  reg_df %>% 
    mutate(mean_per_biom = set_units(mean_per_biom, "kg t_fish-1") %>% drop_units()) %>% 
    ggplot(aes(x = feed_type, y = mean_per_biom, fill = feed_group)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = feed_pal) +
    facet_grid(
      rows = vars(nutrient), 
      scales = "free"
    ) +
    labs(
      x = "Feed digestibility",
      y = expression("Total loading (kg"~t[fish]^{-1}*")")
    ) +
    prettyplot_600() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank()
      )
  })
ls[[1]]


