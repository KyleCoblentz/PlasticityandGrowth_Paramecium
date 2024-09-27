################################################################################
### Phenotypic plasticity analysis -- Bayesian version
################################################################################

### load libraries

library(ggplot2); library(dplyr); library(cowplot); library(ggcorrplot); library(ggbiplot); library(brms); library(cmdstanr)

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













































