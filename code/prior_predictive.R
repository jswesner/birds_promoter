library(tidyverse)
library(janitor)
library(readxl)
library(brms)
library(ggthemes)
library(viridis)

theme_set(theme_default())

brm_small2 = readRDS(file = "models/brm_small2.rds")

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

posts_control = posts %>% filter(cat_gene == "control") %>% 
  select(-cat_gene) %>% 
  expand_grid(cat_gene = c("stress", "social")) %>% 
  mutate(gene = paste0(gene, " (control)")) %>% 
  mutate(median = 0, # to ensure that it will plot at the bottom
         prior_post = "Posterior")

# simulate priors (checked this with brms sample_prior = "only" and it is the same. just faster by hand this way)
n = 1000
rexp_phi = 4
int_sd = 2
b_sd = 0.5

priors = brm_small2$data %>% 
  distinct(carers, gene, cat_gene) %>% 
  expand_grid(.draw = 1:n) %>% 
  mutate(intercept = rnorm(nrow(.), 0, int_sd),
         b = rnorm(nrow(.), 0, b_sd)) %>% 
  # expand_grid(carers = c(0, 1)) %>% 
  mutate(carers_n = case_when(carers == "2" ~ 0,
                              TRUE ~ 1)) %>% 
  mutate(.linpred = intercept + b*carers_n,
         .epred = inv_logit_scaled(.linpred)) %>% 
  rowwise() %>% 
  mutate(gene_int_sd = rexp(1, rexp_phi),
         gene_int_offset = rnorm(1, 0, gene_int_sd),
         gene_b_sd = rexp(1, rexp_phi),
         gene_b_offset = rnorm(1, 0, gene_b_sd)) %>% 
  mutate(cat_gene_int_sd = rexp(1, rexp_phi),
         cat_gene_int_offset = rnorm(1, 0, cat_gene_int_sd),
         cat_gene_b_sd = rexp(1, rexp_phi),
         cat_gene_b_offset = rnorm(1, 0, cat_gene_b_sd)) %>% 
  mutate(ind_int_sd = rexp(1, rexp_phi),
         ind_int_offset = rnorm(1, 0, ind_int_sd)) %>% 
  mutate(re_null_linpred = intercept + ind_int_offset + cat_gene_int_offset +
           gene_int_offset + (b + gene_b_offset + cat_gene_b_offset)*carers_n,
         re_null_epred = inv_logit_scaled(re_null_linpred)) %>% 
  select(carers, gene, cat_gene, .draw, re_null_epred) %>%
  rename(.epred = re_null_epred) %>% 
  group_by(gene) %>%
  mutate(median = median(.epred),
         prior_post = "Prior",
         h_prop_meth_s = 0)

priors_control = priors %>% filter(cat_gene == "control") %>% 
  select(-cat_gene) %>% 
  expand_grid(cat_gene = c("stress", "social")) %>% 
  mutate(gene = paste0(gene, " (control)")) %>% 
  mutate(median = 0, # to ensure that it will plot at the bottom
         prior_post = "Prior")

priors_all = priors %>% filter(cat_gene != "control") %>% 
  bind_rows(priors_control)


# combine and plot
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



