setwd("C:/Users/b284718/OneDrive - University of East Anglia/Emiljan")
library(raster)

habs <- stack("landcover_cropland.tif")
rf <- stack("rf_2019.tif")

#To make a colour ramp centered on zero:

breaks <- seq(-0.305,0.305,0.01) # these are the break points for colours - make sure they span the full range of values across all your encounter rate change rasters, and are symetical around zero! 


cols <- colorRampPalette(c("red4","orange","gold","gray95","lightskyblue","dodgerblue","blue3"))(length(breaks)-1) #this is the colour ramp - change the colours as you like, bu make sure they are symetrical around the "gray95", which represents the zero point!


############
# Below is a pipeline to extract random points from within species range and compare these to habitat change
df<-getValues(rf)
rpoints <- data.frame(xyFromCell(rf, 1:ncell(rf)),getValues(rf))
rpoints <- subset(rpoints,!is.na(layer)) #limit to non-na points
rpoints$layer <- NULL #get rid of column to leave just points
dim(rpoints)

#randomly select 1000 points:
rpoints <- rpoints[sample(1:nrow(rpoints),1000,replace = F),] 

plot(rf)
points(rpoints,cex=0.2)

#now extract data for bird encounter rate change and habitat change
# Example habitat change raster:
hab_change <- habs[[10]]-habs[[9]] #
e_change <- extract(rf, rpoints) #extract cell values for bird change
h_change <- extract(hab_change,rpoints)#extract cell values for habitat change

plot(h_change,e_change)
abline(lm(e_change~h_change))
summary(lm(e_change~h_change))

#NB Worth considering whether to look at overall trend in habitat and encounter rate across the period, and not just year-to-year change!