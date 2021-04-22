
  
agg.2010 <- aggregate(encounter_rate ~ layer, lm_l_2010, mean)
agg.2010$sd <- aggregate(encounter_rate ~ layer, lm_l_2010, stats::sd)$encounter_rate
agg.2010$se <- aggregate(encounter_rate ~ layer, lm_l_2010, std.error)$encounter_rate   
agg.2010 <- agg.2010 %>% mutate(year = 2010)

agg.2011 <- aggregate(encounter_rate ~ layer, lm_l_2011, mean)
agg.2011$sd <- aggregate(encounter_rate ~ layer, lm_l_2011, stats::sd)$encounter_rate
agg.2011$se <- aggregate(encounter_rate ~ layer, lm_l_2011, std.error)$encounter_rate 
agg.2011 <- agg.2011 %>% mutate(year = 2011)

agg.2012 <- aggregate(encounter_rate ~ layer, lm_l_2012, mean)
agg.2012$sd <- aggregate(encounter_rate ~ layer, lm_l_2012, stats::sd)$encounter_rate
agg.2012$se <- aggregate(encounter_rate ~ layer, lm_l_2012, std.error)$encounter_rate 
agg.2012 <- agg.2012 %>% mutate(year = 2012)

agg.2013 <- aggregate(encounter_rate ~ layer, lm_l_2013, mean)
agg.2013$sd <- aggregate(encounter_rate ~ layer, lm_l_2013, stats::sd)$encounter_rate
agg.2013$se <- aggregate(encounter_rate ~ layer, lm_l_2013, std.error)$encounter_rate 
agg.2013 <- agg.2013 %>%mutate(year = 2013)

agg.2014 <- aggregate(encounter_rate ~ layer, lm_l_2014, mean)
agg.2014$sd <- aggregate(encounter_rate ~ layer, lm_l_2014, stats::sd)$encounter_rate
agg.2014$se <- aggregate(encounter_rate ~ layer, lm_l_2014, std.error)$encounter_rate  
agg.2014 <- agg.2014 %>% mutate(year = 2014)

agg.2015 <- aggregate(encounter_rate ~ layer, lm_l_2015, mean)
agg.2015$sd <- aggregate(encounter_rate ~ layer, lm_l_2015, stats::sd)$encounter_rate
agg.2015$se <- aggregate(encounter_rate ~ layer, lm_l_2015, std.error)$encounter_rate   
agg.2015 <- agg.2015 %>%mutate(year = 2015)

agg.2016 <- aggregate(encounter_rate ~ layer, lm_l_2016, mean)
agg.2016$sd <- aggregate(encounter_rate ~ layer, lm_l_2016, stats::sd)$encounter_rate
agg.2016$se <- aggregate(encounter_rate ~ layer, lm_l_2016, std.error)$encounter_rate  
agg.2016 <- agg.2016 %>% mutate(year = 2016)

agg.2017 <- aggregate(encounter_rate ~ layer, lm_l_2017, mean)
agg.2017$sd <- aggregate(encounter_rate ~ layer, lm_l_2017, stats::sd)$encounter_rate
agg.2017$se <- aggregate(encounter_rate ~ layer, lm_l_2017, std.error)$encounter_rate  
agg.2017 <- agg.2017 %>% mutate(year = 2017)

agg.2018<- aggregate(encounter_rate ~ layer, lm_l_2018, mean)
agg.2018$sd <- aggregate(encounter_rate ~ layer, lm_l_2018, stats::sd)$encounter_rate
agg.2018$se <- aggregate(encounter_rate ~ layer, lm_l_2018, std.error)$encounter_rate  
agg.2018 <- agg.2018 %>% mutate(year = 2018)

agg.2019 <- aggregate(encounter_rate ~ layer, lm_l_2019, mean)
agg.2019$sd <- aggregate(encounter_rate ~ layer, lm_l_2019, stats::sd)$encounter_rate
agg.2019$se <- aggregate(encounter_rate ~ layer, lm_l_2019, std.error)$encounter_rate 
agg.2019 <- agg.2019 %>% mutate(year = 2019)



agg.all <- rbind(agg.2010,agg.2011,agg.2012,agg.2013,agg.2014,agg.2015,agg.2016,agg.2017,agg.2018,agg.2019)

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

agg.all <- agg.all %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

agg.all <- agg.all %>% pivot_wider(names_from = ly_name, values_from = encounter_rate)
agg.all <- agg.all[, ]

#split the pland dataframe into layers
test.1.l1 <- test.1[test.1$layer %in% 1,]
test.1.l1$year <- as.integer(test.1.l1$year)

#layer1
layer.1.enc <- agg.all[1:4][complete.cases(agg.all$Prairie_Potholes),]

test.2 <-layer.1.enc %>% inner_join(test.1.l1, by="year")



t1 <- lapply(
  c("pland_00_water", "pland_01_evergreen_needleleaf", "pland_02_evergreen_broadleaf","pland_03_deciduous_needleleaf", "pland_04_deciduous_broadleaf" , "pland_05_mixed_forest",        
    "pland_06_closed_shrubland" ,    "pland_07_open_shrubland"   ,    "pland_08_woody_savanna"   ,     "pland_09_savanna"  ,           
    "pland_10_grassland"     ,       "pland_11_wetland" ,             "pland_12_cropland" ,            "pland_13_urban" ,              
    "pland_14_mosiac"  ,             "pland_15_barren"),
  function(x) summary(lm(as.formula(paste0("layer ~ year + ", x)), test.1.l1))
)


lm1 <- lm(layer ~ year + pland_00_water, test.1.l1)
lm2 <- lm(layer ~ year + pland_01_evergreen_needleleaf, test.1.l1)
lm3 <- lm(layer ~ year + pland_02_evergreen_broadleaf, test.1.l1)
lm4 <- lm(layer ~ year + pland_03_deciduous_needleleaf, test.1.l1)
lm5 <- lm(layer ~ year + pland_04_deciduous_broadleaf, test.1.l1)
lm6 <- lm(layer ~ year + pland_05_mixed_forest, test.1.l1)
lm7 <- lm(layer ~ year + pland_06_closed_shrubland, test.1.l1)
lm8 <- lm(layer ~ year + pland_07_open_shrubland, test.1.l1)
lm9 <- lm(layer ~ year + pland_08_woody_savanna, test.1.l1)
lm10 <- lm(layer ~ year + pland_09_savanna, test.1.l1)
lm11 <- lm(layer ~ year + pland_10_grassland, test.1.l1)
lm12 <- lm(layer ~ year + pland_11_wetland, test.1.l1)
lm13 <- lm(layer ~ year + pland_12_cropland, test.1.l1)
lm14 <- lm(layer ~ year + pland_13_urban, test.1.l1)
lm15 <- lm(layer ~ year + pland_14_mosiac, test.1.l1)
lm16 <- lm(layer ~ year + pland_15_barren, test.1.l1)

mtable(lm1,lm2,lm3,lm4,lm5,lm6,lm7,lm8,lm9,lm10,lm11,lm12,lm13,lm14,lm15,lm16, summary.stats=c("Likelihood-ratio", "N","sigma","R-squared","F","p","N"))


x <- paste0('lm', 1:16)
mtable.test <- as.data.frame(sapply(mtable1[x], `[[`, 'sumstat')) %>% round(3)
setDT(mtable.test, keep.rownames = TRUE)[]
mtable.test <- mtable.test %>% pivot_longer(-c(1))
mtable.test <- mtable.test %>% pivot_wider(names_from = rn, values_from = value)
mtable.test$name <- c(0:15)
colnames(mtable.test)[1] <- "landcover"
mtable.test <- mtable.test %>%
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)
mtable.test <- mtable.test[, c(13, 1:12)]








my.formula <- y ~ x
ggplot(test.2, aes(year,Prairie_Potholes)) +
  geom_line(colour="darkblue") +
  labs(y="Encounter rate in BCR:Prairie Potholes", x = "Years (2010-2019)") +
  geom_smooth(method="lm", se=F, formula=my.formula, level = 0.95) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula
               = my.formula, parse = TRUE) +  geom_errorbar(aes(ymin = Prairie_Potholes -
                                                                  se, ymax = Prairie_Potholes + se), size = 0.5, width=0.2)
