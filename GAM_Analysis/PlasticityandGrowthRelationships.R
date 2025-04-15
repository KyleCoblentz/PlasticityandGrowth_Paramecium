################################################################################
### Morphological plasticity across temperatures and variation in TPCs 
################################################################################

### load libraries

library(dplyr); library(ggplot2); library(cowplot); library(ggbiplot); library(RColorBrewer)

### load data 

# morphology data

morph_data <- read.csv('Plast_StartPop_Data.csv')

# modify data 

morph_data <- morph_data %>% mutate(Line = as.factor(sapply(strsplit(morph_data$file, split = "_"),'[[', 4)),
                                    Temperature = as.numeric(sapply(strsplit(morph_data$file, split = "_"),'[[', 3)),
                                    LineTemp = paste(Line, Temperature, sep = "_"))

morph_data <- filter(morph_data, Temperature < 37)

morph_summ_data <- morph_data %>% group_by(Line, LineTemp, Temperature) %>% 
  select(-c(X, file, id)) %>% summarise_all(.funs = median)

# growth data

growth_data <- read.csv('StartPop_TPC.csv')

# modify growth rate data

growth_data <- growth_data %>% filter(Genotype != 'blank')

growth_data$Genotype <- as.factor(paste0('G', growth_data$Genotype))

colnames(growth_data)[5] <- 'Line'

growth_data <- growth_data %>% filter(!is.na(Growth.Rate.Hours))

growth_data <- growth_data %>% filter(is.finite(Growth.Rate.Hours))

### Try to reduce the dimensionality of trait variation using a PCA

# select only the traits that we 'care' about

morph_summ_data <- morph_summ_data %>% select(c(mean_area, mean_major, mean_minor,
                                                mean_ar, mean_turning, sd_turning,
                                                gross_speed, sd_gross_speed, Line,
                                                Temperature, LineTemp)) %>% ungroup()

morph_pca <- prcomp(x = morph_summ_data %>% select(-c(Line, Temperature, LineTemp)),
       center = TRUE, scale = TRUE)

morph_pca

morph_pca$sdev/sum(morph_pca$sdev)

ggbiplot(morph_pca)

# axis 1 is largely size and axis 2 largely varies with aspect ratio and speed

# this variation matches what we also see across outcrossed lines

pca_scores <- predict(morph_pca, newdata = morph_summ_data %>% select(-c(Line, Temperature, LineTemp)))

morph_summ_data <- cbind(morph_summ_data, pca_scores[,c('PC1', 'PC2')])

# to make PC1 a little more easily interpretable, we can its negative so that 
# larger numbers mean a larger size

morph_summ_data$PC1 <- -morph_summ_data$PC1

# make some graphs of how different phenotypes are changing with temperature

ggplot(data = morph_summ_data, aes(x = Temperature, y = PC1)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, degree = 3, raw = TRUE), se = FALSE)

ggplot(data = morph_summ_data, aes(x = Temperature, y = PC1, group = Line, color = Line)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, degree = 3, raw = TRUE), se = FALSE)

# pca 2 

ggplot(data = morph_summ_data, aes(x = Temperature, y = PC2)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, degree = 3, raw = TRUE), se = FALSE)

ggplot(data = morph_summ_data, aes(x = Temperature, y = PC2, group = Line, color = Line)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ poly(x, degree = 3, raw = TRUE), se = FALSE)

### generally see that size is changing nonlinearly with a minimum at middle 
### temperature. In speed, we see a nonlinear but generally increasing 
### relationship with temperature.

### Now we want to look at whether there are relationships between how 
### morphology is changing with temperature and growth rates

# put together growth and morphology data

growth_morph_data <- left_join(growth_data, morph_summ_data, by = c('Line', 'Temperature'))

# plot pca1 versus growth rate with points color coded by temperature

ggplot(data = growth_morph_data, aes(x = PC1, y = Growth.Rate.Hours, color = Temperature, group = Line)) + geom_point() + 
  scale_color_viridis_c() + geom_smooth(method = 'lm', se = FALSE)

# plot pca2

ggplot(data = growth_morph_data, aes(x = PC2, y = Growth.Rate.Hours, color = Temperature, group = Line)) + geom_point() + 
  scale_color_viridis_c() + geom_smooth(method = 'lm', se = FALSE)

### we generally see that outcrossed lines show relationships between size and 
### growth rate in which increasing temperatures generally lead to increaseing growth rates
### and sizes of individuals within outcrossed lines are typically decreasing. 
### However we have some weird line that does the opposite (gets bigger and grows faster).

### in terms of speed, we see a positive relationship. That is, as temperature
### increases, outcrossed lines get faster, and, as they do, they have a higher
### population growth rate. 

### Our hypothesis here was that the more plastic an outcrossed line is in its traits
### the less variation it should show in its population growth rate. That is, by chaning
### their traits, they are able to buffer the effects of temperature on their growth rates.

### To look at this hypothesis we will begin by asking whether outcrossed lines with greater 
### ranges in their traits also have shallower slopes. In other words, we ask whether lines 
### with greater changes in their traits shower a weaker relationship between variation in the 
### the trait and variation in growth rates across temperatures.

### Start with PC1. Fit a model and pull slopes. Also calculate ranges in size.
### put these in a data frame and see if there is a relationship

PC1_fit <- lm(formula = Growth.Rate.Hours ~ PC1*Line, data = growth_morph_data)

summary(PC1_fit)

# create a data frame with Line, ranges in PC1, and the slope of the fitted relationship
# between PC1 and growth rate across temperatures

PC1_plast_data <- data.frame(Line = unique(growth_morph_data$Line), 
                    Slope = abs(c(coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG120'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG33'],
                                coef(PC1_fit)[2],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG110'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG35'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG41'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG44'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG69'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG16'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG38'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG72'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG19'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG77'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG104'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG89'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG20'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG34'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG30'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG90'],
                                coef(PC1_fit)[2] + coef(PC1_fit)['PC1:LineG36'])))

### get ranges and add

Range_PC1 <- morph_summ_data %>% group_by(Line) %>% summarise(Range = max(PC1) - min(PC1))

PC1_plast_data <- left_join(PC1_plast_data, Range_PC1, by = 'Line')

ggplot(data = PC1_plast_data, aes(x = Range, y = Slope)) + geom_point() + 
  geom_smooth(method = 'lm')

summary(lm(Slope ~ Range, data = PC1_plast_data))

### do the same thing for PC2

PC2_fit <- lm(formula = Growth.Rate.Hours ~ PC2*Line, data = growth_morph_data)

summary(PC2_fit)

# create a data frame with Line, ranges in PC2, and the slope of the fitted relationship
# between PC2 and growth rate across temperatures

PC2_plast_data <- data.frame(Line = unique(growth_morph_data$Line), 
                             Slope = abs(c(coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG120'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG33'],
                                           coef(PC2_fit)[2],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG110'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG35'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG41'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG44'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG69'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG16'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG38'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG72'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG19'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG77'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG104'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG89'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG20'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG34'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG30'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG90'],
                                           coef(PC2_fit)[2] + coef(PC2_fit)['PC2:LineG36'])))


### get ranges and add

Range_PC2 <- morph_summ_data %>% group_by(Line) %>% summarise(Range = max(PC2) - min(PC2))

PC2_plast_data <- left_join(PC2_plast_data, Range_PC2, by = 'Line')

ggplot(data = PC2_plast_data, aes(x = Range, y = Slope)) + geom_point() + 
  geom_smooth(method = 'lm')

summary(lm(Slope ~ Range, data = PC2_plast_data))

### weak evidence for a relationship between slopes and growth rates with 
### size and a bit stronger evidence for a relationship with speed. 

### now we want to be able to connect this with the TPC data more generally.
### specifically, we know that there is some evidence among the outcrossed lines 
### that there is a specialist-generalist tradeoff. What we want to know is whether
### these plastic trait response differences among outcrossed lines can potentially 
### explain that trade off.










































