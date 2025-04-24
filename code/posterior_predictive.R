library(tidyverse)
library(janitor)
library(brms)
library(ggthemes)
library(viridis)

theme_set(theme_default())

brm_small2 = readRDS(file = "models/brm_small2.rds")

post_check = pp_check(brm_small2, type = "dens_overlay_grouped", group = "gene")

post_check + 
  facet_wrap(~group, scales = "free_y")
