################################################################################
### Phenotypic plasticity analysis -- Bayesian version
################################################################################

### load libraries

library(ggplot2); library(dplyr); library(cowplot); library(ggcorrplot); library(ggbiplot); library(brms); library(cmdstanr); library(loo)

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

mod_1_mean_major <- brm(formula = mean_major ~ Temperature + (1|Genotype), data = ind_data)

summary(mod_1_mean_major)

plot(mod_1_mean_major)

# trace plot looks good 

conditional_effects(mod_1_mean_major)

ranef(mod_1_mean_major)

### linear with random intercept and slope

mod_2_mean_major <- brm(formula = mean_major ~ Temperature + (1 + Temperature|Genotype), data = ind_data)

summary(mod_2_mean_major)

plot(mod_2_mean_major)

### quadratic with random intercept

mod_3_mean_major <- brm(formula = mean_major ~ poly(Temperature, 2, raw = TRUE) + (1|Genotype), data = ind_data)

# there is a warning, but things look fine. Could run chains for longer and it will go away

summary(mod_2_mean_major)

plot(mod_3_mean_major)

### quadratic with random intercept and slope

mod_4_mean_major <- brm(formula = mean_major ~ poly(Temperature, 2, raw = TRUE) + (1 + Temperature|Genotype), data = ind_data)

summary(mod_4_mean_major)

plot(mod_4_mean_major)

### quadratic with random intercept, slope, and quadratic term

mod_5_mean_major <- brm(formula = mean_major ~ poly(Temperature, 2, raw = TRUE) + (1 + Temperature + I(Temperature^2)|Genotype), data = ind_data)

### very slow ... but seems fine 

summary(mod_5_mean_major)

plot(mod_5_mean_major)

### cubic with random intercept





























