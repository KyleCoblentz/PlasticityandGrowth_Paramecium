################################################################################
### Growth rate data and plasticity
################################################################################

### load packages

library(ggplot2); library(cowplot); library(cmdstanr); library(dplyr); library(tidyr); library(rTPC); library(nls.multstart)

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

### fit lactin-2 model using rTPC to get reasonable initial values

start_vals <- get_start_vals(x = growth_data$Temperature, y = growth_data$growth_days, model_name = 'lactin2_1995')

lower_vals <- get_lower_lims(x = growth_data$Temperature, y = growth_data$growth_days, model_name = 'lactin2_1995')

upper_vals <- get_upper_lims(x = growth_data$Temperature, y = growth_data$growth_days, model_name = 'lactin2_1995')

# fit lactin-2 model to the entire dataset

nls_multstart(growth_days ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
              data = growth_data,
              iter = c(3,3,3,3),
              start_lower = start_vals - 10,
              start_upper = start_vals + 10,
              lower = lower_vals,
              upper = upper_vals,
              supp_errors = 'Y',
              convergence_count = FALSE)

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

### make function for initial values

init_fun <- function() {
  list(mu_rho = runif(1, 0.01, 0.9),
       rho_i = runif(20, 0.01, 0.9),
       mu_Tmax = runif(1, 20, 50),
       Tmax_i = runif(20, 20, 50),
       mu_delta = runif(1, 0.05, 10),
       delta_i = runif(20, 0.05, 10),
       mu_lambda = runif(1, -10, 1),
       lambda_i = runif(20, -10, 1))
}


fit <- mod$sample(
  data = stan_data, 
  chains = 1, 
  parallel_chains = 1,
  refresh = 100,
  iter_warmup = 10000,
  iter_sampling = 1000,
  init = init_fun)

print(n = 110, fit$summary())

draws_df <- fit$draws(format = "df")

plot(x = 1:1000, y = draws_df$`rho_i[7]`)

hist(growth_data$growth_days)
is.finite(growth_data$Temperature)
