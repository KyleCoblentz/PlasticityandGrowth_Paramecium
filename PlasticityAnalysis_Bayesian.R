################################################################################
### Phenotypic plasticity analysis -- Bayesian version
################################################################################

### load libraries

library(ggplot2); library(dplyr); library(cowplot); library(ggcorrplot); library(ggbiplot); library(brms); library(cmdstanr); library(loo); library(rstan)

### allow loo to use all cores

options(mc.cores = 4)

### ### load data

ind_data <- read.csv('Plast_StartPop_Data.csv')

pca_data <- ind_data %>% select(-c(X, mean_grey, sd_grey, sd_area,
                                   sd_perimeter, sd_major, sd_minor,
                                   sd_ar, duration, N_frames, id))

pca_mvt <- prcomp(pca_data[,-1], center = TRUE, scale = TRUE)

pca_ind_plot <- ggbiplot(pca_mvt, groups = pca_data$file, ellipse = TRUE) + theme(legend.position = 'none')

pca_ind_plot

### look at mean data

mean_pca_data <- pca_data %>% group_by(file) %>% summarize_all(.funs = mean)

mean_pca <- prcomp(mean_pca_data[,-1], center = TRUE, scale = TRUE)

pca_mean_plot <- ggbiplot(mean_pca, labels = mean_pca_data$file, groups = mean_pca_data$file) + theme(legend.position = 'none')

pca_mean_plot

### look at correlation plot for mean data

ggcorrplot(cor(mean_pca_data[,-1]))

################################################################################
### Reduce the number of phenotypes 
################################################################################

### as in our previous work we can probably safely ignore some variables from the beginning. 
### In the copepod foraging paper, we focused on length, width, aspect ratio, mean turning angle,
### gross speed, net displacement, and standard deviation of gross speed

### we can start by looking at PCA's and correlations when we limit our data to just 
### these variables

pca_data <- ind_data %>% select(c(file, mean_major, mean_minor, mean_ar, mean_turning, gross_speed,
                                  net_disp, sd_gross_speed))

pca_mvt <- prcomp(pca_data[,-1], center = TRUE, scale = TRUE)

pca_ind_plot <- ggbiplot(pca_mvt, groups = pca_data$file, ellipse = TRUE) + theme(legend.position = 'none')

pca_ind_plot

### look at mean data

mean_pca_data <- pca_data %>% group_by(file) %>% summarize_all(.funs = mean)

mean_pca <- prcomp(mean_pca_data[,-1], center = TRUE, scale = TRUE)

pca_mean_plot <- ggbiplot(mean_pca, labels = mean_pca_data$file, groups = mean_pca_data$file) + theme(legend.position = 'none')

pca_mean_plot

### look at correlation plot for mean data

ggcorrplot(cor(mean_pca_data[,-1]))

### overall this is pretty similar to our previous findings
### there are some correlations, but the strongest is between 
### aspect ratio and gross speed

### plot of the relationship

ggplot(data = mean_pca_data, aes(x = mean_ar, y = gross_speed)) + geom_point() + geom_smooth(method = 'lm')

### can look at some of the others just to visualize the strengths of the correlations
### for the correlation matrix

ggplot(data = mean_pca_data, aes(x = mean_major, y = mean_minor)) + geom_point() + geom_smooth(method = 'lm')

ggplot(data = mean_pca_data, aes(x = mean_major, y = gross_speed)) + geom_point() + geom_smooth(method = 'lm')

ggplot(data = mean_pca_data, aes(x = mean_turning, y = gross_speed)) + geom_point() + geom_smooth(method = 'lm')

################################################################################
### Plasticity in phenotypes
################################################################################

### for each of the phenotypes of interest, we can examine their relationships with 
### temperature and ask about plasticity using "random regression models" which 
### we can implement in a Bayesian framework to get uncertainty on measures of 
### plasticity for each trait in response to temperature

### modify the data 

### first need to extract temperature and genotype data from the 
### individual level data

ind_data <- ind_data %>% mutate(Temperature = sapply(strsplit(ind_data$file, "_"), function(x) x[3] ),
                                Genotype = sapply(strsplit(ind_data$file, "_"), function(x) x[4]))

### can then get mean data

mean_data <- ind_data %>% group_by(file, Temperature, Genotype) %>% select(-id) %>% summarise_all(.funs = median)

# make temperature numeric

ind_data$Temperature <- as.numeric(ind_data$Temperature)

mean_data$Temperature <- as.numeric(mean_data$Temperature)

# make genotype a factor

ind_data$Genotype <- as.factor(ind_data$Genotype)

mean_data$Genotype <- as.factor(mean_data$Genotype)

# drop 37 since we only have measurements for 2 genotypes

ind_data <- ind_data %>% filter(Temperature < 37)

mean_data <- mean_data %>% filter(Temperature < 37)

# add a variable that is centered temperature

# ind_data <- ind_data %>% mutate(TempCenter = Temperature-mean(unique(Temperature)))

# mean_data <- mean_data %>% ungroup() %>% mutate(TempCenter = Temperature-mean(unique(Temperature)))

################################################################################
### Fitting Bayesian "Random Regression Models"
################################################################################

### start with mean_major (paramecium length)

### first just visualize the relationship with temperature

ggplot(data = ind_data, aes(x = Temperature, y = mean_major, color = Genotype)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k  = 6), se = FALSE)

ggplot(data = mean_data, aes(x = Temperature, y = mean_major, color = Genotype)) + geom_point() + 
  geom_smooth()


### will plan to fit a series of models of increasing complexity.
### for fixed effects we will fit models from linear to cubic
### for random effects will do the intercept and then 
### higher order terms 

### we will then use waic to ask about model support

### simplest model -- linear with random intercept

# for now, we will just use default priors

mod_1_mean_major <- brm(formula = mean_major ~ Temperature + (1|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_1_mean_major)

plot(mod_1_mean_major)

conditional_effects(mod_1_mean_major)

mod_1_mean_major_loo <- loo(mod_1_mean_major)

### linear with random intercept and slope

mod_2_mean_major <- brm(formula = mean_major ~ Temperature + (1 + Temperature|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_2_mean_major)

plot(mod_2_mean_major)

conditional_effects(mod_2_mean_major)

mod_2_mean_major_loo <- loo(mod_2_mean_major)

### quadratic with random intercept

mod_3_mean_major <- brm(formula = mean_major ~ poly(Temperature, 2, raw = TRUE) + (1|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_3_mean_major)

plot(mod_3_mean_major)

conditional_effects(mod_3_mean_major)

mod_3_mean_major_loo <- loo(mod_3_mean_major)

### quadratic with random intercept and slope

mod_4_mean_major <- brm(formula = mean_major ~ poly(Temperature, 2, raw = TRUE) + (1 + poly(Temperature, 1, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_4_mean_major)

plot(mod_4_mean_major)

conditional_effects(mod_4_mean_major)

mod_4_mean_major_loo <- loo(mod_4_mean_major)

### quadratic with random intercept, slope, and quadratic term

mod_5_mean_major <- brm(formula = mean_major ~ poly(Temperature, 2, raw = TRUE) + (1 + poly(Temperature, 2, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_5_mean_major)

plot(mod_5_mean_major)

conditional_effects(mod_5_mean_major)

mod_5_mean_major_loo <- loo(mod_5_mean_major)

### cubic with random intercept

mod_6_mean_major <- brm(formula = mean_major ~ poly(Temperature, 3, raw = TRUE) + (1|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_6_mean_major)

plot(mod_6_mean_major)

conditional_effects(mod_6_mean_major)

mod_6_mean_major_loo <- loo(mod_6_mean_major)

### cubic with random intercept and slope

mod_7_mean_major <- brm(formula = mean_major ~ poly(Temperature, 3, raw = TRUE) + (1 + poly(Temperature, 1, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_7_mean_major)

plot(mod_7_mean_major)

conditional_effects(mod_7_mean_major)

mod_7_mean_major_loo <- loo(mod_7_mean_major)

### cubic with random intercept, slope, and quadratic

mod_8_mean_major <- brm(formula = mean_major ~ poly(Temperature, 3, raw = TRUE) + (1 + poly(Temperature, 2, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_8_mean_major)

plot(mod_8_mean_major)

conditional_effects(mod_8_mean_major)

mod_8_mean_major_loo <- loo(mod_8_mean_major)


### cubic with random intercept, slope, quadratic, and cubic

mod_9_mean_major <- brm(formula = mean_major ~ poly(Temperature, 3, raw = TRUE) + (1 + poly(Temperature, 3, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_9_mean_major)

plot(mod_9_mean_major)

conditional_effects(mod_9_mean_major)

mod_9_mean_major_waic <- waic(mod_9_mean_major)

### this model doesn't fit well ... 

### we can compare all of the models to ask which one is the best at predicting the data

loo_compare(mod_1_mean_major_loo, mod_2_mean_major_loo, mod_3_mean_major_loo,
            mod_4_mean_major_loo, mod_5_mean_major_loo, mod_6_mean_major_loo,
            mod_7_mean_major_loo, mod_8_mean_major_loo)


### next we will want to think about the interpretation of this model and visualize
### the predictions of the model

### plot predictions of the most complicated model.

### will need "newdata" data frame with TempCenter and Genotypes

range(ind_data$Temperature)

### TempCenter ranges from -15.16 to 8.84

### for predictions we can go by 0.2

seq(10, 34, by = 0.2)

# 121 rows for each genotype

Genotype <- rep(unique(ind_data$Genotype), each = 121)

# vectors for temperatures

Temperature <- rep(seq(10, 34, by = 0.2), times = 20)

# put together into a new data frame

new_data <- data.frame(Genotype, Temperature)

predict_mean_major <- posterior_epred(mod_4_mean_major,
                                      newdata = new_data)

# get means for each of the predictions

predict_mean_major <- apply(predict_mean_major, 2, median)

predict_mean_major <- data.frame(Genotype, Temperature, Prediction = predict_mean_major)

# now make the plot

ggplot(data = mean_data, aes(x = Temperature, y = mean_major, color = Genotype)) + 
  geom_point() + geom_line(data = predict_mean_major, aes(x = Temperature, y = Prediction, color = Genotype))

# plot for a single genotype

ggplot(data = filter(ind_data, Genotype == "G33"), aes(x = Temperature, y = mean_major)) + 
  geom_point() + geom_line(data = filter(predict_mean_major, Genotype == "G33"), aes(x = Temperature, y = Prediction))


plot(x = ranef(mod_4_mean_major)$Genotype[,1,2], y = ranef(mod_8_mean_major)$Genotype[,1,2])

### I think maybe fixing it quadratic is fine ... super correlated random effects between 
### 

### now we can move through each of the different phenotypes. Or maybe just focus on size and movement

################################################################################
### Minor axis -- width
################################################################################

### visualize relationship with temperature

ggplot(data = ind_data, aes(x = Temperature, y = mean_minor, color = Genotype)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k  = 6), se = FALSE)

ggplot(data = mean_data, aes(x = Temperature, y = mean_minor, color = Genotype)) + geom_point() + 
  geom_smooth()

### simplest model -- linear with random intercept

# for now, we will just use default priors

mod_1_mean_minor <- brm(formula = mean_minor ~ Temperature + (1|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_1_mean_minor)

plot(mod_1_mean_minor)

conditional_effects(mod_1_mean_minor)

mod_1_mean_minor_loo <- loo(mod_1_mean_minor)

### linear with random intercept and slope

mod_2_mean_minor <- brm(formula = mean_minor ~ Temperature + (1 + Temperature|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_2_mean_minor)

plot(mod_2_mean_minor)

conditional_effects(mod_2_mean_minor)

mod_2_mean_minor_loo <- loo(mod_2_mean_minor)


### quadratic with random intercept

mod_3_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 2, raw = TRUE) + (1|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_3_mean_minor)

plot(mod_3_mean_minor)

conditional_effects(mod_3_mean_minor)

mod_3_mean_minor_loo <- loo(mod_3_mean_minor)


### quadratic with random intercept and slope

mod_4_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 2, raw = TRUE) + (1 + poly(Temperature, 1, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_4_mean_minor)

plot(mod_4_mean_minor)

conditional_effects(mod_4_mean_minor)

mod_4_mean_minor_loo <- loo(mod_4_mean_minor)

### quadratic with random intercept, slope, and quadratic term

mod_5_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 2, raw = TRUE) + (1 + poly(Temperature, 2, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_5_mean_minor)

plot(mod_5_mean_minor)

conditional_effects(mod_5_mean_minor)

mod_5_mean_minor_loo <- loo(mod_5_mean_minor)

### cubic with random intercept

mod_6_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 3, raw = TRUE) + (1|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_6_mean_minor)

plot(mod_6_mean_minor)

conditional_effects(mod_6_mean_minor)

mod_6_mean_minor_loo <- loo(mod_6_mean_minor)

### cubic with random intercept and slope

mod_7_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 3, raw = TRUE) + (1 + poly(Temperature, 1, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_7_mean_minor)

plot(mod_7_mean_minor)

conditional_effects(mod_7_mean_minor)

mod_7_mean_minor_loo <- loo(mod_7_mean_minor)

### cubic with random intercept, slope, and quadratic

mod_8_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 3, raw = TRUE) + (1 + poly(Temperature, 2, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr', cores = getOption("mc.cores", 1))

summary(mod_8_mean_minor)

plot(mod_8_mean_minor)

conditional_effects(mod_8_mean_minor)

mod_8_mean_minor_loo <- loo(mod_8_mean_minor)

### cubic with random intercept, slope, quadratic, and cubic

mod_9_mean_minor <- brm(formula = mean_minor ~ poly(Temperature, 3, raw = TRUE) + (1 + poly(Temperature, 3, raw = TRUE)|Genotype), data = ind_data,
                        backend = 'cmdstanr')

summary(mod_9_mean_minor)

plot(mod_9_mean_minor)

conditional_effects(mod_9_mean_minor)

mod_9_mean_minor_loo <- loo(mod_9_mean_minor)

### compare models 

loo_compare(mod_1_mean_minor_loo, mod_2_mean_minor_loo, mod_3_mean_minor_loo,
            mod_4_mean_minor_loo, mod_5_mean_minor_loo, mod_6_mean_minor_loo,
            mod_7_mean_minor_loo, mod_8_mean_minor_loo)

# put together into a new data frame

new_data <- data.frame(Genotype, Temperature)

predict_mean_minor <- posterior_epred(mod_8_mean_minor,
                                      newdata = new_data)

# get means for each of the predictions

predict_mean_minor <- apply(predict_mean_minor, 2, median)

predict_mean_minor <- data.frame(Genotype, Temperature, Prediction = predict_mean_minor)

# now make the plot

ggplot(data = mean_data, aes(x = Temperature, y = mean_minor, color = Genotype)) + 
  geom_point() + geom_line(data = predict_mean_minor, aes(x = Temperature, y = Prediction, color = Genotype))

# plot for a single genotype

ggplot(data = filter(ind_data, Genotype == "G44"), aes(x = Temperature, y = mean_minor)) + 
  geom_point() + geom_line(data = filter(predict_mean_minor, Genotype == "G44"), aes(x = Temperature, y = Prediction))


plot(x = ranef(mod_7_mean_major)$Genotype[,1,2], y = ranef(mod_7_mean_minor)$Genotype[,1,2])

summary(lm(ranef(mod_8_mean_minor)$Genotype[,1,2] ~ ranef(mod_8_mean_major)$Genotype[,1,2]))










