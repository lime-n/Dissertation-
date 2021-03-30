setwd("C:/Users/b284718/OneDrive - University of East Anglia/Emiljan")
library(raster)
library(sf)


habs <- stack("l.forest_cover.urban.tif")
rf <- stack("output/m.urban.tif")

#To make a colour ramp centered on zero:
rf <- rf[[c(1, 3)]]
map_proj <- st_crs("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")
habs <- projectRaster(habs, crs = map_proj$proj4string, method = "ngb")

habs <- habs[[3:13]]
minus <- function(a, b){
  a-b
}
habs_change <- overlay(habs[[3:11]], habs[[2]], fun=minus)
habs_change <- stack(habs[[1]], habs_change)
names(habs_change) <- c("layer", "y2011","y2012", "y2013", "y2014", "y2015","y2016", "y2017", "y2018", "y2019")

############
# Below is a pipeline to extract random points from within species range and compare these to habitat change
df<-getValues(habs_change)
rpoints <- data.frame(xyFromCell(habs_change, 1:ncell(habs_change)),getValues(habs_change))
rpoints <- subset(rpoints,!is.na(layer)) #limit to non-na points
y2011 <- rpoints[rpoints$y2011 != 0, c(1, 2)]
y2012 <- rpoints[rpoints$y2012 != 0, c(1, 2)]
y2013 <- rpoints[rpoints$y2013 != 0, c(1, 2)]
y2014 <- rpoints[rpoints$y2014 != 0, c(1, 2)]
y2015 <- rpoints[rpoints$y2015 != 0, c(1, 2)]
y2016 <- rpoints[rpoints$y2016 != 0, c(1, 2)]
y2017 <- rpoints[rpoints$y2017 != 0, c(1, 2)]
y2018 <- rpoints[rpoints$y2018 != 0, c(1, 2)]
y2019 <- rpoints[rpoints$y2019 != 0, c(1, 2)]


#randomly select 1000 points:
y2011 <- y2011[sample(1:nrow(y2011),1000,replace = F),] 
y2012 <- y2011[sample(1:nrow(y2012),1000,replace = F),] 
y2013 <- y2011[sample(1:nrow(y2013),1000,replace = F),] 
y2014 <- y2011[sample(1:nrow(y2014),1000,replace = F),] 
y2015 <- y2011[sample(1:nrow(y2015),1000,replace = F),] 
y2016 <- y2011[sample(1:nrow(y2016),1000,replace = F),] 
y2017 <- y2011[sample(1:nrow(y2017),1000,replace = F),] 
y2018 <- y2011[sample(1:nrow(y2018),1000,replace = F),] 
y2019 <- y2011[sample(1:nrow(y2019),1000,replace = F),] 


#now extract data for bird encounter rate change and habitat change
# Example habitat change raster:

e_change.1 <- extract(rf, y2011) #extract cell values for bird change
h_change.1 <- data.frame(extract(habs_change[[c(1, 2)]],y2011))#extract cell values for habitat change
h_change.1$year <- 2011
colnames(h_change.1)[2]<- "water"
y2011$id <- 1:nrow(y2011)
h_change.1$id <- 1:nrow(h_change.1)
h_change.1 <- merge(h_change.1, y2011)

e_change.2 <- extract(rf, y2012) #extract cell values for bird change
h_change.2 <- data.frame(extract(habs_change[[c(1, 3)]],y2012))#extract cell values for habitat change
h_change.2$year <- 2012
colnames(h_change.2)[2]<- "water"
y2012$id <- 1:nrow(y2012)
h_change.2$id <- 1:nrow(h_change.2)
h_change.2 <- merge(y2012, h_change.2)

e_change.3 <- extract(rf, y2013) #extract cell values for bird change
h_change.3 <- data.frame(extract(habs_change[[c(1, 4)]],y2013))#extract cell values for habitat change
h_change.3$year <- 2013
colnames(h_change.3)[2]<- "water"
y2013$id <- 1:nrow(y2013)
h_change.3$id <- 1:nrow(h_change.3)
h_change.3 <- merge(y2013, h_change.3)

e_change.4 <- extract(rf, y2014) #extract cell values for bird change
h_change.4 <- data.frame(extract(habs_change[[c(1, 5)]],y2014))#extract cell values for habitat change
h_change.4$year <- 2014
colnames(h_change.4)[2]<- "water"
y2014$id <- 1:nrow(y2014)
h_change.4$id <- 1:nrow(h_change.4)
h_change.4 <- merge(y2014, h_change.4)

e_change.5 <- extract(rf, y2015) #extract cell values for bird change
h_change.5 <- data.frame(extract(habs_change[[c(1, 6)]],y2015))#extract cell values for habitat change
h_change.5$year <- 2015
colnames(h_change.5)[2]<- "water"
y2015$id <- 1:nrow(y2015)
h_change.5$id <- 1:nrow(h_change.5)
h_change.5 <- merge(y2015, h_change.5)

e_change.6 <- extract(rf, y2016) #extract cell values for bird change
h_change.6 <- data.frame(extract(habs_change[[c(1, 7)]],y2016))#extract cell values for habitat change
h_change.6$year <- 2016
colnames(h_change.6)[2]<- "water"
y2016$id <- 1:nrow(y2016)
h_change.6$id <- 1:nrow(h_change.6)
h_change.6 <- merge(y2016, h_change.6)

e_change.7 <- extract(rf, y2017) #extract cell values for bird change
h_change.7 <- data.frame(extract(habs_change[[c(1, 8)]],y2017))#extract cell values for habitat change
h_change.7$year <- 2017
colnames(h_change.7)[2]<- "water"
y2017$id <- 1:nrow(y2017)
h_change.7$id <- 1:nrow(h_change.7)
h_change.7 <- merge(y2017, h_change.7)

e_change.8 <- extract(rf, y2018) #extract cell values for bird change
h_change.8 <- data.frame(extract(habs_change[[c(1, 9)]],y2018))#extract cell values for habitat change
h_change.8$year <- 2018
colnames(h_change.8)[2]<- "water"
y2018$id <- 1:nrow(y2018)
h_change.8$id <- 1:nrow(h_change.8)
h_change.8 <- merge(y2018, h_change.8)

e_change.9 <- extract(rf, y2019) #extract cell values for bird change
h_change.9 <- data.frame(extract(habs_change[[c(1, 10)]],y2019))#extract cell values for habitat change
h_change.9$year <- 2019
colnames(h_change.9)[2]<- "water"
y2019$id <- 1:nrow(y2019)
h_change.9$id <- 1:nrow(h_change.9)
h_change.9 <- merge(y2019, h_change.9)

r.all.change <- rbind(e_change.1,e_change.2,e_change.3,e_change.4,e_change.5,e_change.6,e_change.7,e_change.8,e_change.9)
h.all.change <- rbind(h_change.1,h_change.2,h_change.3,h_change.4,h_change.5,h_change.6,h_change.7,h_change.8,h_change.9)
c_change <- data.frame(r.all.change, h.all.change)
c_change <- c_change[, -2]
names(c_change) <- c("encounter","id", "layer", "water", "year", "X", "Y")


write.csv(c_change, "t.urban.csv")


#NB Worth considering whether to look at overall trend in habitat and encounter rate across the period, and not just year-to-year change!

