################################################################################
### Growth rate data and plasticity
################################################################################

### load packages

library(ggplot2); library(cowplot); library(cmdstanr); library(dplyr); library(tidyr)

### load growth rate data

growth_data <- read.csv("StartPop_TPC.csv")

### drop blanks

growth_data <- growth_data %>% filter(Genotype != 'blank')

### drop NA's and -Inf

growth_data <- growth_data %>% filter(!is.na(Growth.Rate.Hours) & is.finite(Growth.Rate.Hours))


### make some plots to look at the data

# growth rate in days

growth_data <- growth_data %>% mutate(growth_days = Growth.Rate.Hours*24)

ggplot(data = growth_data, aes(x = Temperature, y = growth_days, color = Genotype)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k = 6), se = FALSE)

### sort growth data by genotype so that things are consistent

growth_data <- growth_data %>% arrange(Genotype)

### need to make genotype a factor for things to work out alright

growth_data$Genotype <- as.factor(growth_data$Genotype)

### need to set up data for the Stan model

stan_data <- list(N = nrow(growth_data),
                  n_genotypes = 20,
                  growth = growth_data$growth_days,
                  temp = growth_data$Temperature,
                  genotype = as.numeric(growth_data$Genotype))

stanmodel <- 'Lactin2Hierarchical.stan'

mod <- cmdstan_model(stan_file = stanmodel)

fit <- mod$sample(
  data = stan_data, 
  chains = 1, 
  parallel_chains = 1,
  refresh = 100,
  iter_warmup = 10000,
  iter_sampling = 1000,
  init = list(list(mu_rho = 0.5, mu_Tmax = 40, mu_delta = -1, mu_lambda = -1))
)

print(n = 100, fit$summary())


draws_df <- fit$draws(format = "df")

exp(0.267*10) - exp(0.267*10 - ((1.1 - 10)/1)) + -1.57

hist(growth_data$growth_days)
is.finite(growth_data$Temperature)
