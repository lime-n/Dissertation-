library(tidyverse)
library(raster)
library(sf)
library(ranger)

#load the necessary data
altitude <- raster("altitudes.tif")
bird <- read.csv("Bird_Dataset_2019.csv")
bird <- bird[, -c(1, 2)]

#use this to rasterise
uk_map <- raster('Uk_map.gri')
uk_map <- crop(uk_map, altitude)

#use this for plotting
uk.map <- getData('GADM', country='GBR', level=3)
uk.map <- crop(uk.map, uk_map)


#Presence/absence
bird$Species_observed <- ifelse(bird$Pop_Index > 0, 1, 0)

#create a neighbourhood buffer radius
neighborhood_radius <- 0.04 * ceiling(max(res(altitude))) / 2

#Select the years of study, presence and absence data along with coordinates, then concatenate into a Polygon then list.
bird_buff <- bird %>% 
  distinct(Year, Lat, Long,Species_observed) %>% 
  # for 2019 use 2018 landcover data
  mutate(year_lc = if_else(as.integer(Year) > 2013, 
                           as.character(2013), as.character(Year)),
         year_lc = paste0("y", year_lc)) %>% 
  # convert to spatial features
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>% 
  # transform to tif projection
  st_transform(crs = projection(altitude)) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(Year, Species_observed, geometry))

#aggregate the lists into a single dataframe
bird_buff = bird_buff$data %>%
  dplyr::bind_rows() %>% #bind the list into a data.frame
  sf::st_cast(to = "POINT") %>% #convert polygon to a list of points
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry() #drop the geometry column


#rasterize the dataset
bird.raster <- bird_buff %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% st_transform(crs = projection(altitude)) %>% rasterize(altitude)

#plot the data
plot(uk.map)
plot(bird.raster$Species_observed, col = c("Red", "Yellow"), legend=TRUE, add = TRUE)

#Or calculate for each individual year:
bird_buff_year <- bird_buff %>% group_by(Year) %>% mutate(id = row_number()) %>% pivot_wider(names_from = Year, values_from = Species_observed, values_fill = list(Species_observed = 0)) %>% select(-id)
bird.raster_year <- bird_buff_year %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% st_transform(crs = projection(altitude)) %>% rasterize(altitude)

#Then plot their changes; use a function or simply minus them
minus_function <- function(a, b){
  a - b
}
#overlay allows to apply a function to the raster
bird_changes <- overlay(bird.raster_year$X2013, bird.raster_year$X1994, fun=minus_function)

plot(uk.map)
plot(bird_changes,add=TRUE, col=c("yellow","red", "blue") )


#The rise of the Machines!

#Prepare the dataset for training
ebird_split <- bird %>% 
  # select only the columns to be used in the model
  select(Species_observed, Year, Woodland, Scrubland, Grassland, Heath, Farmland, Urban, Waterbodies, Coastal,
         starts_with("Temp"),
         starts_with("Precip"))%>% 
  drop_na()
# split 80/20
ebird_split <- ebird_split %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(ebird_split, nrow)

detection_freq <- mean(ebird_split$train$Species_observed)

ebird_split$train$Species_observed <- factor(ebird_split$train$Species_observed)
# grow random forest

#Check 1. and 2. before the random forest model. Decide which option to take, depending on the option, before proceding you can either:



rf <- ranger(formula =  Species_observed ~ ., 
             data = ebird_split$train,
             importance = "impurity",
             probability = TRUE,
             replace = TRUE, 
             sample.fraction = c(detection_freq, detection_freq),)

#create a dataframe containing the environmental variables contained in the random forest algorithim

environment <- bird[, c(2,3, 4, 7:38)] %>% drop_na()

#create a prediction model; This is the next stage of machine learning relative to the ranger package

predict_rf <- predict(rf, data = environment, type = "response")

#bind the results with the environment dataset
ML_dataset <- bind_cols(environment, encounter_rate = predict_rf$predictions) %>% 
  select(Lat, Long, encounter_rate ) %>% 
  mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))


#create a raster
ML_raster <- ML_dataset %>% 
  # convert to spatial features
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>% 
  st_transform(crs = projection(uk_map)) %>% 
  # rasterize
  rasterize(uk_map)

ML_raster <- ML_raster[[2]]

plot(uk.map)
plot(ML_raster$encounter_rate, add = TRUE, legend = FALSE)
