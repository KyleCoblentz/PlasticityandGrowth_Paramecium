################################################################################
### analysis of phenotypic plasticity data
################################################################################

### load packages

library(ggplot2); library(dplyr); library(cowplot); library(ggcorrplot); library(ggbiplot); library(mgcv); library(gratia)

### load data

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

### looks fairly similar to previous data

### want to look at relationships with temperature 

### first need to extract temperature and genotype data from the 
### individual level data

ind_data <- ind_data %>% mutate(Temperature = sapply(strsplit(ind_data$file, "_"), function(x) x[3] ),
                                Genotype = sapply(strsplit(ind_data$file, "_"), function(x) x[4]))

### can then get mean data

mean_data <- ind_data %>% group_by(file, Temperature, Genotype) %>% select(-id) %>% summarise_all(.funs = median)

### look at some visualizations of how phenotypes might be changing with genotypes across temperature

# start with some of the size related variables

# mean_major

# first make temperature numeric

ind_data$Temperature <- as.numeric(ind_data$Temperature)

mean_data$Temperature <- as.numeric(mean_data$Temperature)

# make genotype a factor

ind_data$Genotype <- as.factor(ind_data$Genotype)

# drop 37 since we only have measurements for 2 genotypes

ind_data <- ind_data %>% filter(Temperature < 37)

mean_data <- mean_data %>% filter(Temperature < 37)


ggplot(data = ind_data, aes(x = Temperature, y = mean_major, color = Genotype, group = Genotype)) + geom_point() +
  geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k = 5), se = FALSE)

ggplot(data = mean_data, aes(x = Temperature, y = mean_major)) + geom_point() + geom_smooth(method = 'lm', formula = y ~ poly(x, 4, raw = T))

gam_major <- gam(log(mean_major) ~ s(Temperature, bs = 'tp', k = 5, by = Genotype) + Genotype, data = ind_data, 
                 method = 'REML')

summary(gam_major)

gam_plot_major <- draw(gam_major, residuals = TRUE)

gam_plot_alt <- ggplot(data = ind_data, aes(x = Temperature, y = mean_major, group = Genotype)) + geom_point(alpha = 0.25) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = 'tp', k = 5)) + 
  facet_wrap(facets = "Genotype") + theme_cowplot()

# look at area

area_plot <- ggplot(data = ind_data, aes(x = Temperature, y = mean_area, group = Genotype)) + 
  geom_line(stat = 'smooth', method = 'gam', formula = y ~ s(x, bs = 'tp', k = 5), color = 'black', alpha = 0.5) + 
  theme_cowplot() + ylab("Cell Area")

# look at mean_minor

ggplot(data = ind_data, aes(x = Temperature, y = mean_minor, color = Genotype)) + geom_point()

ggplot(data = mean_data, aes(x = Temperature, y = mean_minor, color = Genotype)) + geom_point() + geom_line()

gam_minor <- gam(log(mean_minor) ~ s(Temperature, bs = 'tp', k = 5, by = Genotype) + Genotype, data = ind_data,
                 method = 'REML')

summary(gam_minor)

draw(gam_minor, residuals = TRUE)

# look at mean_ar

ggplot(data = ind_data, aes(x = Temperature, y = mean_ar, color = Genotype)) + geom_point()

ggplot(data = mean_data, aes(x = Temperature, y = mean_ar)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, degree = 3, raw = T))

gam_ar <-  gamm(mean_ar ~ s(Temperature, bs = 'tp', k = 5) + s(Genotype, bs = 're'), data = ind_data,
               method = 'REML')

gam_ar

draw(gam_ar, residuals = TRUE)



### check out speed then

ggplot(data = ind_data, aes(x = Temperature, y = gross_speed, color = Genotype)) + geom_point()

ggplot(data = mean_data, aes(x = Temperature, y = gross_speed, color = Genotype)) + geom_point() + geom_line()

gam_gross_speed <- gam(log(gross_speed) ~ s(Temperature, bs = 'tp', k = 5, by = Genotype) + Genotype, data = ind_data,
                       method = 'REML')

summary(gam_gross_speed)

draw(gam_gross_speed, residuals = TRUE)

### turning 

ggplot(data = ind_data, aes(x = Temperature, y = mean_turning, color = Genotype)) + geom_point()

ggplot(data = mean_data, aes(x = Temperature, y = mean_turning, color = Genotype)) + geom_point() + geom_line()

gam_mean_turning <- gam(mean_turning ~ s(Temperature, bs = 'tp', k = 5, by = Genotype) + Genotype, data = ind_data,
                       method = 'REML')

summary(gam_mean_turning)

draw(gam_mean_turning, residuals = TRUE)

################################################################################
### use random regression models to look at gxe and estimate plasticity
################################################################################

################################################################################
### mean_major
################################################################################

library(lme4)

mod1_mean_major <- lmer(mean_major ~ Temperature + (1|Genotype), data = ind_data)

summary(mod_mean_major)

ggplot(data = ind_data, aes(x = Temperature, y = mean_major, color = Genotype)) + geom_point() 
  
mod2_mean_major <- lmer(mean_major ~ poly(Temperature,2, raw = T) + (1|Genotype), data = ind_data)

summary(mod2_mean_major)

AIC(mod1_mean_major, mod2_mean_major)

mod3_mean_major <- lmer(mean_major ~ poly(Temperature, 3, raw = T) + (1|Genotype), data = ind_data)

summary(mod3_mean_major)

AIC(mod1_mean_major, mod2_mean_major, mod3_mean_major)

mod4_mean_major <- lmer(mean_major ~ poly(Temperature, 4, raw = T) + (1|Genotype), data = ind_data)

AIC(mod1_mean_major, mod2_mean_major, mod3_mean_major, mod4_mean_major)

mod3_mean_major_rslope <- lmer(mean_major ~ poly(Temperature, 3, raw = T) + (1 + Temperature|Genotype), data = ind_data)

summary(mod3_mean_major_rslope2)

AIC(mod1_mean_major, mod2_mean_major, mod3_mean_major, mod4_mean_major, mod3_mean_major_rslope)

Temp_pred <- rep(unique(ind_data$Temperature), times = length(unique(ind_data$Genotype)))

Genotype_pred <- rep(unique(ind_data$Genotype), each = length(unique(ind_data$Temperature)))

new_data <- as.data.frame(cbind(Temp_pred, as.character(Genotype_pred)), stringsAsFactors = TRUE)

colnames(new_data) <- c('Temperature', 'Genotype')

new_data$Temperature <- as.numeric(paste(new_data$Temperature))

mod3_mean_major_pred <- data.frame(Genotype = new_data$Genotype,
                                   Temperature = new_data$Temperature,
                                   mean_major = predict(mod3_mean_major_rslope, newdata = new_data))

ggplot(data = mod3_mean_major_pred, aes(x = Temperature, y = mean_major, color = Genotype)) + geom_line() +
  geom_point(data = mean_data, aes(x = Temperature, y = mean_major, color = Genotype))






