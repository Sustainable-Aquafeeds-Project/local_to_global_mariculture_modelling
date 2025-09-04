library(magrittr)
library(stringr)
library(tidyr)
library(dplyr)
library(arrow)
library(ggplot2)
library(units)
library(conflicted)
conflicts_prefer(dplyr::select(), dplyr::filter())

fixnum <- function(n, digits = 4) {str_flatten(c(rep("0", digits-nchar(as.character(n))), as.character(n)))}
fixup <- function(df) {
  df %>% 
    mutate(sd_rel = sd/mean) %>% 
    group_by(farm_ID, feed, days) %>% 
    reframe(mean = sum(mean),
            sd_rel = sum(sd_rel)) %>% 
    mutate(sd = sd_rel * mean,
           feed = factor(feed, levels = c("past", "reference", "future")))
}

this_species <- "atlantic_salmon"
prod_path <- file.path("data", this_species, "data_products")
fig_path <- file.path("data", this_species, "figures")

farm_IDs <- targets::tar_read(farm_IDs, store = "05_targets_farm")
farm_static_data <- targets::tar_read(farm_static_data, store = "05_targets_farm")
farm_outs <- list.files(file.path(prod_path, "model_outputs_farm"), full.names = T)

# Farm biomass plots
for (fid in farm_IDs) {
  farm_biom <- farm_outs %>% 
    str_subset("biomass") %>% 
    str_subset(fixnum(fid)) %>% 
    read_parquet() %>% 
    fixup() %>% 
    mutate(
      mean = set_units(set_units(mean, "g"), "t"),
      sd = set_units(set_units(sd, "g"), "t")
    )
  farm_t <- farm_static_data %>% 
    filter(farm_id == fid) %>% 
    pull(tonnes_per_farm) %>% 
    set_units("t")
  
  p <- ggplot(farm_biom, aes(x = days, y = mean, ymin = mean-sd, ymax = mean+sd, colour = feed, fill = feed)) +
    geom_ribbon(alpha = 0, linetype = "dotted", linewidth = 0.6) +
    geom_line(linewidth = 0.75) +
    geom_hline(aes(yintercept = farm_t), linetype = "dashed") +
    theme_classic() +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0,550,150)) +
    theme(legend.position = "none",
          aspect.ratio = 0.75) +
    labs(x = "Day of production", y = "Farm biomass")
  
  ggsave(
    plot = p, 
    filename = file.path(fig_path, "farm_biom", paste0("farm_biom_", fixnum(fid), ".png")),
    height = 150, width = 200, units = "mm"
  )
}

# Farm input plots - carbs
for (fid in farm_IDs) {
  farm_biom <- farm_outs %>% str_subset("biomass") %>% str_subset(fixnum(fid)) %>% read_parquet() %>% 
    fixup() %>% 
    mutate(mean = set_units(set_units(mean, "g"), "t"),
           sd = set_units(set_units(sd, "g"), "t"))
  excre <- farm_outs %>% str_subset("C_excr") %>% str_subset(fixnum(fid)) %>% read_parquet() 
  uneat <- farm_outs %>% str_subset("C_uneat") %>% str_subset(fixnum(fid)) %>% read_parquet()
  input <- merge(excre, uneat, by = c("farm_ID", "feed", "days")) %>% 
    mutate(mean = mean.x + mean.y,
           sd = mean * ((sd.x/mean.x) + (sd.y/mean.y))) %>% 
    select(-c(mean.x, mean.y, sd.x, sd.y)) %>% 
    fixup() %>% 
    mutate(mean = set_units(set_units(mean, "g"), "kg"),
           sd = set_units(set_units(sd, "g"), "kg")) %>% 
    merge(farm_biom, by = c("farm_ID", "feed", "days")) %>% 
    mutate(mean = mean.x/mean.y,
           sd = (sd_rel.x/sd_rel.y) * mean,
           mean = set_units(mean, "kg/t"),
           sd = set_units(sd, "kg/t")) %>% 
    select(-c(mean.x, mean.y, sd.x, sd.y, sd_rel.x, sd_rel.y))
  
  p <- ggplot(input, aes(x = days, y = mean, ymin = mean, ymax = mean+sd, colour = feed, fill = feed)) +
    geom_ribbon(alpha = 0, linetype = "dotted", linewidth = 0.6) +
    geom_line(linewidth = 0.75) +
    theme_classic() +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0,550,150)) +
    theme(legend.position = "none",
          aspect.ratio = 0.75) +
    labs(x = "Day of production", y = "Total excreted carbohydrates")
  
  ggsave(
    plot = p, 
    filename = file.path(fig_path, "input_C", paste0("input_C_", fixnum(fid), ".png")),
    height = 150, width = 200, units = "mm"
  )
}

# Farm input plots - lipids
for (fid in farm_IDs) {
  farm_biom <- farm_outs %>% str_subset("biomass") %>% str_subset(fixnum(fid)) %>% read_parquet() %>% 
    fixup() %>% 
    mutate(mean = set_units(set_units(mean, "g"), "t"),
           sd = set_units(set_units(sd, "g"), "t"))
  excre <- farm_outs %>% str_subset("L_excr") %>% str_subset(fixnum(fid)) %>% read_parquet() 
  uneat <- farm_outs %>% str_subset("L_uneat") %>% str_subset(fixnum(fid)) %>% read_parquet()
  input <- merge(excre, uneat, by = c("farm_ID", "feed", "days")) %>% 
    mutate(mean = mean.x + mean.y,
           sd = mean * ((sd.x/mean.x) + (sd.y/mean.y))) %>% 
    select(-c(mean.x, mean.y, sd.x, sd.y)) %>% 
    fixup() %>% 
    mutate(mean = set_units(set_units(mean, "g"), "kg"),
           sd = set_units(set_units(sd, "g"), "kg")) %>% 
    merge(farm_biom, by = c("farm_ID", "feed", "days")) %>% 
    mutate(mean = mean.x/mean.y,
           sd = (sd_rel.x/sd_rel.y) * mean,
           mean = set_units(mean, "kg/t"),
           sd = set_units(sd, "kg/t")) %>% 
    select(-c(mean.x, mean.y, sd.x, sd.y, sd_rel.x, sd_rel.y))
  
  p <- ggplot(input, aes(x = days, y = mean, ymin = mean, ymax = mean+sd, colour = feed, fill = feed)) +
    geom_ribbon(alpha = 0, linetype = "dotted", linewidth = 0.6) +
    geom_line(linewidth = 0.75) +
    theme_classic() +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0,550,150)) +
    theme(legend.position = "none",
          aspect.ratio = 0.75) +
    labs(x = "Day of production", y = "Total excreted lipids")
  
  ggsave(
    plot = p, 
    filename = file.path(fig_path, "input_L", paste0("input_L_", fixnum(fid), ".png")),
    height = 150, width = 200, units = "mm"
  )
}

# Farm input plots - proteins
for (fid in farm_IDs) {
  farm_biom <- farm_outs %>% str_subset("biomass") %>% str_subset(fixnum(fid)) %>% read_parquet() %>% 
    fixup() %>% 
    mutate(mean = set_units(set_units(mean, "g"), "t"),
           sd = set_units(set_units(sd, "g"), "t"))
  excre <- farm_outs %>% str_subset("P_excr") %>% str_subset(fixnum(fid)) %>% read_parquet() 
  uneat <- farm_outs %>% str_subset("P_uneat") %>% str_subset(fixnum(fid)) %>% read_parquet()
  input <- merge(excre, uneat, by = c("farm_ID", "feed", "days")) %>% 
    mutate(mean = mean.x + mean.y,
           sd = mean * ((sd.x/mean.x) + (sd.y/mean.y))) %>% 
    select(-c(mean.x, mean.y, sd.x, sd.y)) %>% 
    fixup() %>% 
    mutate(mean = set_units(set_units(mean, "g"), "kg"),
           sd = set_units(set_units(sd, "g"), "kg")) %>% 
    merge(farm_biom, by = c("farm_ID", "feed", "days")) %>% 
    mutate(mean = mean.x/mean.y,
           sd = (sd_rel.x/sd_rel.y) * mean,
           mean = set_units(mean, "kg/t"),
           sd = set_units(sd, "kg/t")) %>% 
    select(-c(mean.x, mean.y, sd.x, sd.y, sd_rel.x, sd_rel.y))
  
  p <- ggplot(input, aes(x = days, y = mean, ymin = mean, ymax = mean+sd, colour = feed, fill = feed)) +
    geom_ribbon(alpha = 0, linetype = "dotted", linewidth = 0.6) +
    geom_line(linewidth = 0.75) +
    theme_classic() +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(0,550,150)) +
    theme(#legend.position = "none",
          aspect.ratio = 0.75) +
    labs(x = "Day of production", y = "Total excreted protein")
  
  ggsave(
    plot = p, 
    filename = file.path(fig_path, "input_P", paste0("input_P_", fixnum(fid), ".png")),
    height = 150, width = 200, units = "mm"
  )
}


