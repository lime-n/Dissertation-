library(tidyverse)
t_barren <- read.csv("t.barren.csv")
t_open_shrubland <- read.csv("t.open_shrubland.csv")
t_closed_shrubland <- read.csv("t.closed_shrubland.csv")
t_cropland <- read.csv("t.cropland.csv")
t_deciduous_needleleaf <- read.csv("t.deciduous_needleleaf.csv")
t_deciduous_broadleaf <- read.csv("t.deciduous_broadleaf.csv")
t_evergreen_broadleaf <- read.csv("t.evergreen_broadleaf.csv")
t_evergreen_needleleaf <- read.csv("t.evergreen_needleleaf.csv")
t_grassland <- read.csv("t.grassland.csv")
t_mixed_forest <- read.csv("t.mixed_forest.csv")
t_mosaic <- read.csv("t.mosaic.csv")
t_urban <- read.csv("t.urban.csv")
t_wetland <- read.csv("t.wetland.csv")
t_savanna <- read.csv("t.savanna.csv")
t_woody_savanna <- read.csv("t.woody_savanna.csv")
t_water <- read.csv("t.water.csv")


names(t_barren)[5] <- "barren"
names(t_open_shrubland)[5] <- "open_shrubland"
names(t_closed_shrubland)[5] <- "closed_shrubland"
names(t_cropland)[5] <- "cropland"
names(t_deciduous_needleleaf)[5] <- "deciduous_needleleaf"
names(t_deciduous_broadleaf)[5] <- "deciduous_broadleaf"
names(t_evergreen_broadleaf)[5] <- "evergreen_broadleaf"
names(t_evergreen_needleleaf)[5] <- "evergreen_needleleaf"
names(t_grassland)[5] <- "grassland"
names(t_mixed_forest)[5] <- "mixed_forest"
names(t_mosaic)[5] <- "mosaic"
names(t_urban)[5] <- "urban"
names(t_wetland)[5] <- "wetland"
names(t_savanna)[5] <- "savanna"
names(t_woody_savanna)[5] <- "woody_savanna"


t_barren <-t_barren[, -1]
t_open_shrubland <- t_open_shrubland[, -1]
t_closed_shrubland <- t_closed_shrubland[, -1]
t_cropland <- t_cropland[, -1]
t_deciduous_needleleaf <- t_deciduous_needleleaf[, -1]
t_deciduous_broadleaf <- t_deciduous_broadleaf[, -1]
t_evergreen_broadleaf <- t_evergreen_broadleaf[, -1]
t_evergreen_needleleaf <- t_evergreen_needleleaf[, -1]
t_grassland <- t_grassland[, -1]
t_mixed_forest <- t_mixed_forest[, -1]
t_mosaic <- t_mosaic[, -1]
t_urban <- t_urban[, -1]
t_wetland <-t_wetland [, -1]
t_savanna <- t_savanna[, -1]
t_woody_savanna <- t_woody_savanna[, -1]
t_water <- t_water[, -1]




t_barren <- t_barren %>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_closed_shrubland <- t_closed_shrubland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_cropland <- t_cropland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_deciduous_needleleaf <- t_deciduous_needleleaf%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_deciduous_broadleaf <- t_deciduous_broadleaf%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_evergreen_broadleaf <- t_evergreen_broadleaf%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_evergreen_needleleaf <- t_evergreen_needleleaf%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_grassland <- t_grassland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_mixed_forest <- t_mixed_forest%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_mosaic <- t_mosaic%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_open_shrubland <- t_open_shrubland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_urban <- t_urban%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_wetland <- t_wetland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_savanna <- t_savanna%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_woody_savanna <-t_woody_savanna%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_water <- t_water %>% pivot_longer(-c(1, 2, 3, 5, 6, 7))


test.all <- rbind(t_barren,t_open_shrubland,t_closed_shrubland,t_cropland,t_deciduous_needleleaf,t_deciduous_broadleaf,t_evergreen_broadleaf,t_evergreen_needleleaf,t_grassland,t_mixed_forest,t_mosaic,t_open_shrubland,t_urban,t_wetland,t_savanna,t_woody_savanna,t_water)
test.all <- test.all[, -2]
test.all <- test.all[, c(4, 5, 3, 2, 6, 7, 1)]

colnames(test.all) <- c("X", "Y", "year", "layer", "landcover", "pland_change", "encounter")

#library(missRanger)
#test.all.1<- missRanger(test.all, pmm.k = 3, num.trees = 100)

ly_names <- tibble(layer = 1:17, 
                   ly_name = c("Prairie_Potholes"     ,
                               "Boreal_Hardwood_Transition"   ,  
                               "Lower_Great_Lakes_St.Lawrence_Plain",
                               "Badlands_And_Prairies"    ,    
                               "Shortgrass_Prairie"    ,      
                               "Central_Mixed_Grass_Prairie"   ,      
                               "Oaks_And_Prairies"      ,   
                               "Eastern_Tallgrass_Prairie"  ,
                               "Prairie_Hardwood_Transition" ,       
                               "Central_Hardwoods"     ,  
                               "West_Gulf_Coastal_Plain_Ouachitas",
                               "Mississippi_Alluvial_Valley" ,        
                               "Southeastern_Coastal_Plain" ,   
                               "Appalachian_Mountains" ,        
                               "Piedmont"      ,                      
                               "New_England_Mid_Atlantic_Coast" , 
                               "Peninsular_Florida" ))

test.all <- test.all %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

test.all <- test.all[complete.cases(test.all),]
test.all.2 <- test.all[test.all$pland_change <=0,]
test.all.3 <- test.all[test.all$pland_change >=0,]

library(stepfor)

tbest <- stepfor(test.all$encounter, test.all[, -1], alpha=0.2)

shape.pivot <- do.call(cbind, unname(Map(function(x, z) {
  tmp <- as.data.frame(model.matrix(~x - 1))
  if (ncol(tmp) == 1 & class(tmp[[1]]) == "numeric") {
    names(tmp) <- paste0(names(tmp), z)
  }
  tmp
}, test.all, names(test.all))))
names(shape.pivot) <- sub("^x", "", names(shape.pivot))

names(tbest$model)[-c(1, 29, 28, 27, 26)]
my.formula <- reformulate(names(tbest$model)[-c(1, 29, 28, 27, 26)], 'sqrt(encounter)')
bbest <- rlm(my.formula, shape.pivot)


#Estimate the fitted values, residuals, Cook?s distances and the leverages
rjack<-rstudent(bbest)
yhat<-fitted(bbest)
rstud<-rstandard(bbest)
h<-hatvalues(bbest)
d<-cooks.distance(bbest)


# checking for normality
# statistical test


shapiro.test(rjack[sample(5000)])
# graphical methods
par(mfrow=c(1,2))
qqnorm(rjack[sample(5000)])
qqline(rjack[sample(5000)])
hist(rjack,xlab="Jackknife residuals",main="Jackknife residuals")
graphics.off()
# checking for constant error variance
plot(yhat,rjack,xlab="Predicted values", ylab="Jackknife residuals")
abline(h=0)
abline(h=2,lty=2)
abline(h=-2,lty=2)
rjack[sample(5000)]identify(yhat,rjack,id)

# 5. for detecting and reporting outliers and influential observations
# see the lecture's notes.
test.hat <- data.frame(yhat = yhat, rjack = rjack)
test.hat$id <- 1:nrow(test.hat)
test.hat <- test.hat[test.hat$rjack < 3 & test.hat$rjack >-3,]
shape.pivot$id <- 1:nrow(shape.pivot)
shape.pivot <- shape.pivot[shape.pivot$id %in% test.hat$id,]
