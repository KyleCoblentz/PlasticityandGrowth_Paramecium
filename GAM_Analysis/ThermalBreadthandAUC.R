################################################################################
### Code to examine thermal breadth and the potential relationships to 
### morphological plasticity
################################################################################

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(mgcv); library(gratia); library(ggbiplot)

### load phenotype data

phenotype_data <- read.csv('Plast_StartPop_Data.csv')

phenotype_data <- phenotype_data %>% mutate(Line = as.factor(sapply(strsplit(phenotype_data$file, split = "_"),'[[', 4)),
                      Temperature = as.numeric(sapply(strsplit(phenotype_data$file, split = "_"),'[[', 3)),
                      LineTemp = paste(Line, Temperature, sep = "_"))

phenotype_data <- filter(phenotype_data, Temperature < 37)

### load tpc data

tpc_data <- read.csv('StartPop_TPC.csv')

tpc_data <- tpc_data %>% filter(!is.na(Growth.Rate.Hours) & !is.infinite(Growth.Rate.Hours))

tpc_data <- tpc_data %>% filter(Genotype != 'blank')

tpc_data$Genotype <- as.factor(paste0('G', tpc_data$Genotype))

colnames(tpc_data)[which(colnames(tpc_data) == 'Genotype')] <- 'Line'

tpc_data$Date.paramecium.introduced <- as.factor(tpc_data$Date.paramecium.introduced)

### plot of tpc data

ggplot(data = tpc_data, aes(x = Temperature, y = Growth.Rate.Hours, color = Line)) + 
  geom_point() + geom_smooth(method = 'gam', formula = y ~ s(x, k = 6, bs = 'tp'), se = FALSE)

### fit a GAM to TPC data

tpc_gam <- gam(formula = Growth.Rate.Hours ~ Line + s(Temperature, k = 6, bs = 'tp', by = Line) + s(Date.paramecium.introduced, bs = 're'), data = tpc_data,
               method = 'REML')

summary(tpc_gam)

draw(tpc_gam)

### now we want to write a function that uses the GAM predictions to calculate
### the thermal breadth (the range of temperatures for which the growth rate 
### is at least x% of the maximum growth rate). We might as well also pull 
### the information about what that maximum growth rate is and the the 
### temperature at which it occurs (the topt).

tpc_summary <- function(gam_fit, factor, temp_low, temp_high, breadth_proportion){
  
  out_data <- matrix(nrow = 20, ncol = 8) # columns are Line, max_growth, t_opt, t_min, t_max, thermal_range_min, thermal_range_max, thermal_range
  
  for(i in 1:length(unique(factor))){
    
    line_prediction_data <- data.frame(Line = rep(unique(factor)[i], times = length(seq(temp_low, temp_high, by = 0.1))),
                                       Date.paramecium.introduced = rep('2024-02-16', times = length(seq(temp_low, temp_high, by = 0.1))),
                                       Temperature = seq(temp_low, temp_high, by = 0.1))
    
    gam_prediction <- predict.gam(gam_fit, newdata = line_prediction_data, exclude = "s(Date.paramecium.introduced)")
    
    prediction_together <- cbind(line_prediction_data, gam_prediction)
    
    max_growth <- max(prediction_together$gam_prediction)
    
    t_opt <- prediction_together$Temperature[which(prediction_together$gam_prediction == max_growth)]
    
    below_opt <- prediction_together[which(prediction_together$Temperature < t_opt),]
    
    t_min <- below_opt$Temperature[which(abs(below_opt$gam_prediction) == min(abs(below_opt$gam_prediction)))]
    
    above_opt <- prediction_together[which(prediction_together$Temperature > t_opt),]
    
    t_max <- above_opt$Temperature[which(abs(above_opt$gam_prediction) == min(abs(above_opt$gam_prediction)))]
    
    in_range <- filter(prediction_together, gam_prediction >= breadth_proportion*max_growth)
    
    thermal_range_min <- in_range$Temperature[1]
    
    thermal_range_max <- in_range$Temperature[nrow(in_range)]
    
    thermal_range <- thermal_range_max - thermal_range_min
    
    out_data[i,1] <- as.character(unique(factor)[i])
    
    out_data[i,2] <- max_growth
    
    out_data[i,3] <- t_opt
    
    out_data[i,4] <- t_min
    
    out_data[i,5] <- t_max
    
    out_data[i,6] <- thermal_range_min
    
    out_data[i,7] <- thermal_range_max
    
    out_data[i,8] <- thermal_range
    
    print(i)
    
  }
  
  out_data <- as.data.frame(out_data)
  
  colnames(out_data) <- c("Line", 'max_growth', 't_opt', 't_min',
                          't_max', 'thermal_range_min', 'thermal_range_max', 'thermal_range')
  
  return(out_data)
  
  
}

### see if function works

tpc_summary_data <- tpc_summary(gam_fit = tpc_gam, factor = tpc_data$Line, temp_low = 0, temp_high = 45, breadth_proportion = 0.5)

# need to fix how line is reported ... but everything else seems fine

ggplot(data = filter(tpc_summary_data, max_growth > 0.06), aes(x = thermal_range, y = max_growth)) + geom_point() + geom_smooth(method = 'lm')

ggplot(data = tpc_summary_data, aes(x = thermal_range, y = max_growth)) + geom_point() + geom_smooth(method = 'lm')

summary(lm(max_growth ~ thermal_range, data = filter(tpc_summary_data, max_growth > 0.06)))

### load morphology data and calculate plasticity in morphology and movement to compare
### to the thermal breadth

morph_data <- read.csv('Plast_StartPop_Data.csv')

morph_data <- morph_data %>% mutate(Line = as.factor(sapply(strsplit(morph_data$file, split = "_"),'[[', 4)),
                                                      Temperature = as.numeric(sapply(strsplit(morph_data$file, split = "_"),'[[', 3)),
                                                      LineTemp = paste(Line, Temperature, sep = "_"))

morph_data <- filter(morph_data, Temperature < 37)

### summarise morphological data to get medians at each temperature

med_morph_data <- morph_data %>% group_by(Line, Temperature, LineTemp) %>%
  select(-c(X, file, mean_grey, sd_grey, sd_area, sd_perimeter, sd_major,
            sd_minor, sd_ar, duration, N_frames, id)) %>%
  summarise_all(.funs = median)

### do the same thing for the growth data

med_morph_growth <- full_join(med_morph_data, tpc_data, by = c('Temperature', 'Line'))

med_morph_growth <- filter(med_morph_growth, Temperature < 37)

### now also make a growth_cell_area variable

med_morph_growth <- mutate(med_morph_growth, growth_cell_area = Growth.Rate.Hours*mean_area)

### plot tpcs in terms of cell area

ggplot(data = med_morph_growth, aes(x = Temperature, y = growth_cell_area, group = Line, color = Line)) + geom_point() +
  geom_smooth(method  = 'gam', formula = y ~ s(x, k = 5, bs = 'tp'), se = FALSE)

area_fit <- gam(formula = growth_cell_area ~ Line + s(Temperature, bs = 'tp', k = 6, by = Line), data = med_morph_growth, method = 'REML')

summary(area_fit)

area_fit_no_int <- gam(formula = growth_cell_area ~ Line + s(Temperature, bs = 'tp', k = 6), data = med_morph_growth, method = 'REML')

summary(area_fit_no_int)

area_fit_no_Line <- gam(formula = growth_cell_area ~ s(Temperature, bs = 'tp', k = 6), data = med_morph_growth, method = 'REML')

summary(area_fit_no_Line)

AIC(area_fit, area_fit_no_int, area_fit_no_Line)

### look for relationship between sizes, ar, and speed and growth rates for each of the temperatures

ggplot(data = filter(med_morph_growth, Temperature == 10), aes(x = gross_speed, y = Growth.Rate.Hours)) + geom_point() + 
  geom_smooth(method = 'lm')

### for each of the morphological variables, calculate how much they change across temperature and other summary statistics

ggplot(data = filter(morph_data, Line == 'G102'), aes(x = Temperature, y = mean_area)) + geom_point() + 
  geom_smooth()

### get data for area

summ_area <- morph_data %>% group_by(Line, Temperature) %>% summarise(area = median(mean_area)) %>% ungroup() %>% group_by(Line) %>% 
  summarise(mean_area = median(area),
            max_area = max(area),
            min_area = min(area),
            range_area = max_area - min_area,
            stand_range_area = range_area/max_area)

### get data for length

summ_length <- morph_data %>% group_by(Line, Temperature) %>% summarise(length = median(mean_major)) %>% ungroup() %>% group_by(Line) %>% 
  summarise(mean_length = median(length),
            max_length = max(length),
            min_length = min(length),
            range_length = max_length - min_length,
            stand_range_length = range_length/max_length)

### predictions for width

summ_width <- morph_data %>% group_by(Line, Temperature) %>% summarise(width = median(mean_minor)) %>% ungroup() %>% group_by(Line) %>% 
  summarise(mean_width = mean(width),
            max_width = max(width),
            min_width = min(width),
            range_width = max_width - min_width,
            stand_range_width = range_width/max_width)

### predictions for aspect ratio

summ_ar <- morph_data %>% group_by(Line, Temperature) %>% summarise(ar = median(mean_ar)) %>% ungroup() %>% group_by(Line) %>% 
  summarise(mean_ar = mean(ar),
            max_ar = max(ar),
            min_ar = min(ar),
            range_ar = max_ar - min_ar,
            stand_range_ar = range_ar/max_ar)

### predictions for speed

summ_speed <- morph_data %>% group_by(Line, Temperature) %>% summarise(speed = median(gross_speed)) %>% ungroup() %>% group_by(Line) %>% 
  summarise(mean_speed = mean(speed),
            max_speed = max(speed),
            min_speed = min(speed),
            range_speed = max_speed - min_speed,
            stand_range_speed = range_speed/max_speed)

### put these together into a single data frame

plast_morph_data <- full_join(summ_area, summ_length, by = "Line")

plast_morph_data <- full_join(plast_morph_data, summ_width, by = "Line")

plast_morph_data <- full_join(plast_morph_data, summ_ar, by = "Line")

plast_morph_data <- full_join(plast_morph_data, summ_speed, by = "Line")

### put this data together with the tpc summary data

plast_morph_tpc <- full_join(plast_morph_data, tpc_summary_data, by = 'Line')

ggplot(data= plast_morph_tpc, aes(x = stand_range_speed, y = as.numeric(thermal_range))) + geom_point() + 
  geom_smooth(method = 'lm')

### try with area of hull in pca

pca_1 <- prcomp(med_morph_data[,c(4,6,7,8,9, 10, 11, 15, 19)], center = TRUE, scale = TRUE)

ggbiplot(pca_1)

pca2 <- pca(med_morph_data[,c(4,6,7,8,9, 10, 11, 15, 19)], scale = TRUE)

summary(pca2)

biplot(pca2)

ordihull(pca2, groups = med_morph_data$Line)

hull_data <- t(summary(ordihull(pca2, groups = med_morph_data$Line)))

hull_data <- cbind(rownames(hull_data), data.frame(hull_data, row.names=NULL))

colnames(hull_data)[1] <- 'Line'

plast_growth_data <- full_join(plast_morph_tpc, hull_data, by = 'Line')

ggplot(data = filter(plast_growth_data, Line != 'G20'), aes(x = range_speed, y = as.numeric(t_max) - as.numeric(t_min))) + geom_point() + geom_smooth(method = 'lm')

summary(lm(range_growth ~ Area, data = plast_growth_data))










