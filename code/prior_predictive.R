library(tidyverse)
library(janitor)
library(readxl)
library(brms)
library(ggthemes)
library(viridis)

theme_set(theme_default())

brm_small2 = readRDS(file = "models/brm_small2.rds")
brm_small2_prior = readRDS(file = "models/brm_small2_prior.rds")

# plot model --------------------------------------------------------------
library(tidybayes)
library(ggridges)

posts = brm_small2$data %>% 
  distinct(carers, cat_gene, gene) %>% 
  mutate(h_prop_meth_s = 0) %>% 
  add_epred_draws(brm_small2, re_formula = ~ (1 + carers|cat_gene/gene)) %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  # pivot_wider(names_from = carers, values_from = .epred) %>%
  # mutate(diff = `3+` - `2`) %>%
  group_by(gene) %>%
  mutate(median = median(.epred),
         prior_post = "Posterior") 

priors = brm_small2_prior$data %>% 
  distinct(carers, cat_gene, gene) %>% 
  mutate(h_prop_meth_s = 0) %>% 
  add_epred_draws(brm_small2_prior, re_formula = ~ (1 + carers|cat_gene/gene)) %>% 
  ungroup() %>% 
  select(-.row, -.chain, -.iteration) %>% 
  # pivot_wider(names_from = carers, values_from = .epred) %>%
  # mutate(diff = `3+` - `2`) %>%
  group_by(gene) %>%
  mutate(median = median(.epred),
         prior_post = "Prior") 

posts_control = posts %>% filter(cat_gene == "control") %>% 
  select(-cat_gene) %>% 
  expand_grid(cat_gene = c("stress", "social")) %>% 
  mutate(gene = paste0(gene, " (control)")) %>% 
  mutate(median = 0, # to ensure that it will plot at the bottom
         prior_post = "Posterior")

priors_control = priors %>% filter(cat_gene == "control") %>% 
  select(-cat_gene) %>% 
  expand_grid(cat_gene = c("stress", "social")) %>% 
  mutate(gene = paste0(gene, " (control)")) %>% 
  mutate(median = 0, # to ensure that it will plot at the bottom
         prior_post = "Prior")

posts_all = posts %>% filter(cat_gene != "control") %>% 
  bind_rows(posts_control)
priors_all = priors %>% filter(cat_gene != "control") %>% 
  bind_rows(priors_control)

posts_priors_all = bind_rows(posts_all, priors_all) %>% 
  mutate(cat_gene = str_to_sentence(cat_gene),
         carers = str_to_sentence(carers))

meth_plot_prior = posts_priors_all %>%
  # filter(prior_post == "Prior") %>% 
  ggplot(aes(x = .epred, y = reorder(gene, median), fill = carers, 
             alpha = prior_post)) +
  stat_slab(normalize = "groups", aes(group = interaction(carers, gene, prior_post))) +
  facet_wrap(~cat_gene, scales = "free_y") +
  # viridis::scale_fill_viridis() +
  ggthemes::scale_fill_colorblind() +
  scale_alpha_manual(values = c(0.8, 0.2)) +
  labs(x = "Proportion Methylated",
       y = "Gene") +
  guides(alpha = "none")

ggsave(meth_plot_prior, file = "plots/meth_plot_prior.jpg",
       width = 7, height = 5)

