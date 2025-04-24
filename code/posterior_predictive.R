library(tidyverse)
library(janitor)
library(brms)
library(ggthemes)
library(viridis)

theme_set(theme_default())

brm_small2 = readRDS(file = "models/brm_small2.rds")

post_check = pp_check(brm_small2, type = "hist")

ggsave(post_check + theme(axis.text.x = element_text(size = 8)), 
       file = "plots/post_check.jpg", width = 6.5, height = 5)
