################################################################################
### try to fit inidivudual genotypes to the lactin 2 model
################################################################################

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

ggplot(data = filter(growth_data, Genotype == '89'), aes(x = Temperature, y = growth_days, color = Genotype)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k = 6), se = FALSE)

### sort growth data by genotype so that things are consistent

growth_data <- growth_data %>% arrange(Genotype)

### fit lactin-2 model using rTPC to get reasonable initial values

growth_data_g89 <- filter(growth_data, Genotype == '89')

start_vals <- get_start_vals(x = growth_data_g89$Temperature, y = growth_data_g89$growth_days, model_name = 'lactin2_1995')

lower_vals <- get_lower_lims(x = growth_data_g89$Temperature, y = growth_data_g89$growth_days, model_name = 'lactin2_1995')

upper_vals <- get_upper_lims(x = growth_data_g89$Temperature, y = growth_data_g89$growth_days, model_name = 'lactin2_1995')

# fit lactin-2 model to G102

nls_multstart(growth_days ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
              data = growth_data_g89,
              iter = c(3,3,3,3),
              start_lower = start_vals - 10,
              start_upper = start_vals + 10,
              lower = lower_vals,
              upper = upper_vals,
              supp_errors = 'Y',
              convergence_count = FALSE)
### need to set up data for the Stan model

stan_data <- list(N = nrow(growth_data_g120),
                  growth = growth_data_g120$growth_days,
                  temperature = growth_data_g120$Temperature)

stanmodel <- 'Lactin2.stan'

mod <- cmdstan_model(stan_file = stanmodel)

### make function for initial values

init_fun <- function() {
  list(rho = runif(1, 0.01, 0.1),
       Tmax = runif(1, 20, 50),
       delta = runif(1, 0.05, 5),
       lambda = runif(1, -10, 1),
       sigma = runif(1, 0.0001, 1))
}


fit <- mod$sample(
  data = stan_data, 
  chains = 3, 
  parallel_chains = 3,
  refresh = 100,
  iter_warmup = 10000,
  iter_sampling = 1000,
  init = init_fun,
  adapt_delta = 0.99)

fit$summary()

draws_df <- fit$draws(format = "df")

plot(x = 1:1000, y = draws_df$delta_t[1:1000], type = 'l')
lines(x = 1:1000, y = draws_df$rho[1001:2000], col = 'red')
lines(x = 1:1000, y = draws_df$rho[2001:3000], col = 'green')

plot(x = seq(10, 37, by = 0.1), lactin2_1995(temp = seq(10, 37, by = 0.1), a = 0.04, b = -1.25, tmax = 40, delta_t =  6))
