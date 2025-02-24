################################################################################
### GAM analysis of morphological and fitness plasticity with tempeature
### in Paramecium
################################################################################

### load packages 

library(dplyr); library(ggplot2); library(gratia); library(cowplot); library(mgcv)

### load data

morph_data <- read.csv(file = 'Plast_StartPop_Data.csv')

### create columns for line and temperature

morph_data <- morph_data %>% mutate(Line = as.factor(sapply(strsplit(morph_data$file, split = "_"),'[[', 4)),
                                      Temperature = as.numeric(sapply(strsplit(morph_data$file, split = "_"),'[[', 3)),
                                    LineTemp = paste(Line, Temperature, sep = "_"))

morph_data <- filter(morph_data, Temperature < 37)


### start with cell area 

ggplot(data = morph_data, aes(x = Temperature, y = log(mean_area), color = Line)) + 
         geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

### fit gam to area data

gam_area <- gam(formula = mean_area ~ Line + s(Temperature, k = 5, bs = 'tp', by = Line), data = morph_data, method = 'REML')

summary(gam_area)

draw(gam_area, residuals = TRUE)

ggplot(data = morph_data, aes(x = Temperature, y = mean_area)) + geom_point()

newdata_area <- data.frame(Temperature = rep(seq(10, 34, by = 0.1),20),
                           Line = rep(unique(morph_data$Line), each = 241))

area_predictions <- predict.gam(gam_area, newdata = newdata_area)

newdata_area <- cbind(newdata_area, area_predictions)

area_predictions <- newdata_area %>% group_by(Line) %>% 
  summarise(mean_area = mean(area_predictions),
            min_area = min(area_predictions),
            max_area = max(area_predictions),
            range_area = max_area-min_area,
            stand_range_area = range_area/max_area)

### length

ggplot(data = morph_data, aes(x = Temperature, y = mean_major, color = Line)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

### fit gam to area data

gam_length <- gam(formula = mean_major ~ Line + s(Temperature, k = 5, bs = 'tp', by = Line), data = morph_data, method = 'REML')

summary(gam_length)

draw(gam_length, residuals = TRUE)

newdata_length <- data.frame(Temperature = rep(seq(10, 34, by = 0.1),20),
                           Line = rep(unique(morph_data$Line), each = 241))

length_predictions <- predict.gam(gam_length, newdata = newdata_length)

newdata_length <- cbind(newdata_length, length_predictions)

length_predictions <- newdata_length %>% group_by(Line) %>% 
  summarise(mean_length = mean(length_predictions),
            min_length = min(length_predictions),
            max_length = max(length_predictions),
            range_length = max_length-min_length, 
            stand_range_length = range_length/max_length)

### width

ggplot(data = morph_data, aes(x = Temperature, y = mean_minor, color = Line)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

### fit gam to area data

gam_width <- gam(formula = mean_minor ~ Line + s(Temperature, k = 5, bs = 'tp', by = Line), data = morph_data, method = 'REML')

summary(gam_width)

draw(gam_width, residuals = TRUE)

newdata_width <- data.frame(Temperature = rep(seq(10, 34, by = 0.1),20),
                            Line = rep(unique(morph_data$Line), each = 241))

width_predictions <- predict.gam(gam_width, newdata = newdata_width)

newdata_width <- cbind(newdata_width, width_predictions)

width_predictions <- newdata_width %>% group_by(Line) %>% 
  summarise(mean_width = mean(width_predictions),
            min_width = min(width_predictions),
            max_width = max(width_predictions),
            range_width = max_width-min_width, 
            stand_range_width = range_width/max_width)

### aspect ratio 

ggplot(data = morph_data, aes(x = Temperature, y = mean_ar, color = Line)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

### fit gam to area data

gam_ar <- gam(formula = mean_ar ~ Line + s(Temperature, k = 5, bs = 'tp', by = Line), data = morph_data, method = 'REML')

summary(gam_ar)

draw(gam_ar, residuals = TRUE)

newdata_ar <- data.frame(Temperature = rep(seq(10, 34, by = 0.1),20),
                         Line = rep(unique(morph_data$Line), each = 241))

ar_predictions <- predict.gam(gam_ar, newdata = newdata_ar)

newdata_ar <- cbind(newdata_ar, ar_predictions)

ar_predictions <- newdata_ar %>% group_by(Line) %>% 
  summarise(mean_ar = mean(ar_predictions),
            min_ar = min(ar_predictions),
            max_ar = max(ar_predictions),
            range_ar = max_ar-min_ar, 
            stand_range_ar = range_ar/max_ar)

### speed 


ggplot(data = morph_data, aes(x = Temperature, y = gross_speed, color = Line)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

### fit gam to area data

gam_speed <- gam(formula = gross_speed ~ Line + s(Temperature, k = 5, bs = 'tp', by = Line), data = morph_data, method = 'REML')

summary(gam_speed)

draw(gam_speed, residuals = TRUE)

ggplot(data = morph_data, aes(x = Temperature, y = gross_speed)) + geom_point()

newdata_speed <- data.frame(Temperature = rep(seq(10, 34, by = 0.1),20),
                         Line = rep(unique(morph_data$Line), each = 241))

speed_predictions <- predict.gam(gam_speed, newdata = newdata_speed)

newdata_speed <- cbind(newdata_speed, speed_predictions)

speed_predictions <- newdata_speed %>% group_by(Line) %>% 
  summarise(mean_speed = mean(speed_predictions),
            min_speed = min(speed_predictions),
            max_speed = max(speed_predictions),
            range_speed = max_speed-min_speed, 
            stand_range_speed = range_speed/max_speed)





### load TPC data

tpc_data <- read.csv(file = 'StartPop_TPC.csv')

tpc_data <- tpc_data %>% filter(Genotype != 'blank')

tpc_data$Genotype <- as.factor(paste0('G', tpc_data$Genotype))

colnames(tpc_data)[5] <- 'Line'

# plot tpc data

ggplot(data = tpc_data, aes(x = Temperature, y = Growth.Rate.Hours, color = Line)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

# fit gam to growth data

tpc_data <- tpc_data %>% filter(!is.na(Growth.Rate.Hours))

tpc_data <- tpc_data %>% filter(!is.infinite(Growth.Rate.Hours))

gam_growth <- gam(formula = Growth.Rate.Hours ~ Line + s(Temperature, k = 5, bs = 'tp', by = Line), data = tpc_data, method = 'REML')

summary(gam_growth)

draw(gam_growth, residuals = TRUE)

newdata_tpc <- data.frame(Temperature = rep(seq(10, 34, by = 0.1),20),
                          Line = rep(unique(morph_data$Line), each = 241))

predict_tpc <- predict.gam(gam_growth, newdata = newdata_tpc)

newdata_tpc <- cbind(newdata_tpc, predict_tpc)

tpc_predictions <- newdata_tpc %>% group_by(Line) %>% 
  summarise(mean_growth = mean(predict_tpc),
            min_growth = min(predict_tpc),
            max_growth = max(predict_tpc),
            range_growth = max_growth-min_growth,
            stand_range_growth = range_growth/max_growth)

tpc_descriptions <- newdata_tpc %>% group_by(Line) %>% 
  summarise(topt = Temperature[which.max(predict_tpc)],
            ctmin = Temperature[which.min(abs(predict_tpc[which(Temperature < 20)] - 0))],
            ctmax = Temperature[which(Temperature>32)][which.min(abs(predict_tpc[which(Temperature > 32)] - 0))],
            thermal_range = ctmax - ctmin)



col_data <- full_join(area_predictions, tpc_predictions, by = "Line")

col_data <- full_join(col_data, length_predictions, by = 'Line')

col_data <- full_join(col_data, width_predictions, by = "Line")

col_data <- full_join(col_data, ar_predictions, by = "Line")

col_data <- full_join(col_data, speed_predictions, by = "Line")

col_data <- full_join(col_data, tpc_descriptions, by = "Line")

ggplot(data = col_data, aes(x = stand_range_area, y = stand_range_growth)) + geom_point() +
  geom_smooth(method = 'lm')

plot(x = seq(10, 34, by = 0.1), 
     y = select(filter(newdata_length, Line == 'G102'), length_predictions)$length_predictions/(max(select(filter(newdata_length, Line == 'G102'), length_predictions)$length_predictions)))

points(select(filter(newdata_ar, Line == 'G102'), ar_predictions)$ar_predictions/max(select(filter(newdata_ar, Line == 'G102'), ar_predictions)$ar_predictions))

plot(x = select(filter(newdata_length, Line == 'G69'), length_predictions)$length_predictions, 
     y = select(filter(newdata_tpc, Line == 'G69'), predict_tpc)$predict_tpc)

points(select(filter(newdata_length, Line == 'G102'), length_predictions)$length_predictions/(max(select(filter(newdata_length, Line == 'G102'), length_predictions)$length_predictions)))

### want to create a data set with normalized data for all of the variables to look at correlations

### can potentially also just extract the predictions for each varaible at particular temperatures
### and then look at correlations.























