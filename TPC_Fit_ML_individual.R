################################################################################
### Fit TPCs for each of the genotypes seperately using rTPC package
################################################################################

### load packages

library(ggplot2); library(cowplot); library(dplyr); library(tidyr); library(rTPC); library(nls.multstart)

### load data

### load growth rate data

growth_data <- read.csv("StartPop_TPC.csv")

### drop blanks

growth_data <- growth_data %>% filter(Genotype != 'blank')

### drop NA's and -Inf

growth_data <- growth_data %>% filter(!is.na(Growth.Rate.Hours) & is.finite(Growth.Rate.Hours))


### make some plots to look at the data

# growth rate in days

growth_data <- growth_data %>% mutate(growth_days = Growth.Rate.Hours*24)

ggplot(data = filter(growth_data, Genotype == '34'), aes(x = Temperature, y = growth_days, color = Genotype)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k = 6), se = FALSE)


### want to fit the Lactin-2 model to each of the tpc's for each genotype and 
### then save the output into a file so that we can compare plasticity in 
### in morphology/movement to how much fitness changes with temperature

### I think we should be able to write a loop that goes through, 
### fits the tpc for each genotype, and then saves the resulting fit parameters

### first, we will just do a single fit to see what we need to do

### we will just fit to the entire dataset 

### fit lactin-2 model using rTPC to get reasonable initial values

start_vals <- get_start_vals(x = growth_data$Temperature, y = growth_data$growth_days, model_name = 'lactin2_1995')

lower_vals <- get_lower_lims(x = growth_data$Temperature, y = growth_data$growth_days, model_name = 'lactin2_1995')

upper_vals <- get_upper_lims(x = growth_data$Temperature, y = growth_data$growth_days, model_name = 'lactin2_1995')

# fit lactin-2 model to all of the data

all_fit <- nls_multstart(growth_days ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
              data = growth_data,
              iter = c(3,3,3,3),
              start_lower = start_vals - 10,
              start_upper = start_vals + 10,
              lower = lower_vals,
              upper = upper_vals,
              supp_errors = 'Y',
              convergence_count = FALSE)

### this should give us the parameter estimates

all_fit$m$getPars()

### need to write a loop that goes through all of the datasets and saves the 
### resulting parameter estimates

### create a dataframe to save the data in

tpc_fit_data <- as.data.frame(matrix(nrow = 20, ncol = 5))

colnames(tpc_fit_data) <- c('Genotype', 'rho', 'theta', 'tmax', 'delta_t')

### convert genotype in the growth data to a factor

growth_data$Genotype <- as.factor(growth_data$Genotype)

### go ahead and fill in the genotype column of the new data

tpc_fit_data$Genotype <- levels(growth_data$Genotype)


### set up the for loop 

for(i in 1:length(unique(growth_data$Genotype))){
  
  data <- growth_data %>% filter(as.numeric(Genotype) == i)
  
  ### get starting values and upper and lower bounds
  
  start_vals <- get_start_vals(x = data$Temperature, y = data$growth_days, model_name = 'lactin2_1995')
  
  lower_vals <- get_lower_lims(x = data$Temperature, y = data$growth_days, model_name = 'lactin2_1995')
  
  upper_vals <- get_upper_lims(x = data$Temperature, y = data$growth_days, model_name = 'lactin2_1995')
  
  ### fit the model
  
  fit <- nls_multstart(growth_days ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                data = data,
                iter = c(3,3,3,3),
                start_lower = start_vals - 10,
                start_upper = start_vals + 10,
                lower = lower_vals,
                upper = upper_vals,
                supp_errors = 'Y',
                convergence_count = FALSE)
  
  ### save the parameter estimates to the ith row in our loop
  
  tpc_fit_data[i, 2:5] <-  fit$m$getPars()
  
}



### just want to see whether relationships exist within the TPC data

pairs(tpc_fit_data[,2:5])

summary(lm(rho ~ log(delta_t), data = tpc_fit_data))

ggplot(data = tpc_fit_data, aes(x = delta_t, y = log(rho))) + 
  geom_point() + geom_smooth(method = 'lm')

summary(lm(theta ~ rho, data = tpc_fit_data))

ggplot(data = tpc_fit_data, aes(x = log(rho), y = theta)) + 
  geom_point() + geom_smooth(method = 'lm')




























