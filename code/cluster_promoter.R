library(brms)

brm_small2 = readRDS(file = "models/brm_small2.rds")
brm_small2 = update(brm_small2, iter = 2000, chains = 4)
saveRDS(brm_small2, file = "models/brm_small2.rds")