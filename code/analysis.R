library(tidyverse)
library(janitor)
library(readxl)
library(brms)
library(ggthemes)
library(viridis)

theme_set(theme_default())

# Gene: this is a gene code for each gene we’re looking in the promoter of
# Cat_gene: this is just social/stress (and one control gene- I wouldn’t be opposed to getting rid of some of the genes if there doesn’t seem to be anything (others I feel like have to be included regardless) but we’re hoping for nothing in the control gene, so that will stay no matter what)
# Site: this is just a reference for me to tell me where on the genome it is
# Distance: location of each CpG site relative to the start of the gene. This is probably important because what we’re looking for is methylation at each of these sites (the prop methylation is at this site specifically)
# % meth F/H: you can ignore these (and even delete them really)
# Individual: 😊
# f_prop_meth: this is the proportion of methylation at that site (distance) at fledging
# h_prop_meth: same, but at hatching
# group: this is the group that chick was raised by- should be a random factor since there is some replication here
# carers: number of carers in each group taking care of that chick- binary as 2 or more than 2

promoter_data <- read_excel("data/promoter_data.xlsx") %>% 
  clean_names() %>% 
  mutate(f_prop_meth = parse_number(f_prop_meth),
         h_prop_meth = parse_number(h_prop_meth)) %>% 
  ungroup %>% 
  mutate(distance_s = scale(distance),
         h_prop_meth_s = scale(h_prop_meth),
         f_prop_meth_s = case_when( f_prop_meth == 0 ~ 0.001,
                                    f_prop_meth == 1 ~ 0.999,
                            TRUE ~  f_prop_meth)) %>% 
  filter(!is.na(carers)) 

# brm_small2 = brm(f_prop_meth_s ~ carers + h_prop_meth_s + (1 + carers | cat_gene/gene) + (1 | individual),
#                  data = promoter_data,
#                  family = Beta(link = "logit"),
#                  prior = c(prior(normal(0, 2), class = Intercept),
#                            prior(normal(0, 0.5), class = b),
#                            prior(exponential(4), class = sd)),iter = 2000, chains = 4,
#                  cores = 4)

# saveRDS(brm_small2, file = "models/brm_small2.rds")

brm_small2 = readRDS(file = "models/brm_small2.rds")

brm_small2_prior = update(brm_small2, sample_prior = "only", newdata = brm_small2$data %>% 
                            group_by(gene, cat_gene, carers) %>% sample_n(1), iter = 500, chains = 1)

saveRDS(brm_small2_prior, file = "models/brm_small2_prior.rds")

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
  mutate(median = median(.epred)) 

posts_control = posts %>% filter(cat_gene == "control") %>% 
  select(-cat_gene) %>% 
  expand_grid(cat_gene = c("stress", "social")) %>% 
  mutate(gene = paste0(gene, " (control)")) %>% 
  mutate(median = 0) # to ensure that it will plot at the bottom

posts_all = posts %>% filter(cat_gene != "control") %>% 
  bind_rows(posts_control) %>% 
  mutate(cat_gene = str_to_sentence(cat_gene),
         carers = str_to_sentence(carers))

meth_plot = posts_all %>%
  ggplot(aes(x = .epred, y = reorder(gene, median), fill = carers)) +
  geom_density_ridges() +
  facet_wrap(~cat_gene, scales = "free_y") +
  # viridis::scale_fill_viridis() +
  ggthemes::scale_fill_colorblind() +
  labs(x = "Proportion Methylated",
       y = "Gene")

ggsave(meth_plot, file = "plots/meth_plot.jpg",
       width = 7, height = 5)


meth_plot_diff = posts_all %>% 
  pivot_wider(names_from = carers, values_from = .epred) %>%
  mutate(diff = `2` - `3+`) %>%
  group_by(gene) %>%
  mutate(median = median(diff)) %>%
  ggplot(aes(x = diff, y = reorder(gene, median), fill = median)) +
  geom_density_ridges() +
  facet_wrap(~cat_gene, scales = "free_y") +
  viridis::scale_fill_viridis(direction = 1) +
  # ggthemes::scale_fill_colorblind() +
  labs(x = "Difference in Methylation (`2` - `3+`)",
       y = "Gene") +
  geom_vline(xintercept = 0) +
  guides(fill = "none")

ggsave(meth_plot_diff, file = "plots/meth_plot_diff.jpg",
       width = 7, height = 5)


# summarize ---------------------------------------------------------------

median_table = posts_all %>% 
  group_by(carers, cat_gene, gene) %>% 
  median_qi(.epred) %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  rename(prop_methylated = .epred) %>% 
  select(-.width, -.point, -.interval)

diff_table = posts_all  %>% 
  pivot_wider(names_from = carers, values_from = .epred) %>%
  mutate(diff = `2` - `3+`) %>%  
  group_by(cat_gene, gene) %>% 
  median_qi(diff) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  rename(difference_2_minus_3plus = diff) %>% 
  select(-.width, -.point, -.interval)

prob_table = posts_all  %>% 
  pivot_wider(names_from = carers, values_from = .epred) %>%
  mutate(diff = `2` - `3+`) %>%  
  group_by(cat_gene, gene) %>%
  reframe(prob_diff = sum(diff>0)/max(.draw)) %>% 
  mutate(definition = "Probability that methylation is higher in 2 versus 3+ carers")

write_csv(median_table, file = "tables/median_table.csv")
write_csv(diff_table, file = "tables/diff_table.csv")
write_csv(prob_table, file = "tables/prob_table.csv")
