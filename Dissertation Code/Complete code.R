#Creating the BCR map used in my methodology:
library(tidyverse)
library(broom)
library(rgdal)
#load the geopackage
p <- readOGR(dsn = "bcr.gpkg")
#turn the geopackage into a tibble
tidy_bcr <- tidy(p)
#create a plot to see what it looks like
ggplot(p, aes(x = long, y = lat, group = group)) +
  geom_polygon(color = "black", size = 0.1, fill = "lightgrey") +
  coord_equal() +
  theme_minimal()

p$id <- row.names(p)
#add bcr_name to the dataframe
tidy_bcr <- left_join(tidy_bcr, p@data)
#select bcr regions
bcrcop <- data.frame(bcr_name = sort(p@data$bcr_name),
                     bcr_code = c(28,17,12,24,19,22,13,26,30,21,31,29,23,11,18,27,25))
tidy_bcr <- left_join(tidy_bcr, bcrcop)

#label the regions so they centre on the map
bcrLabel <- tidy_bcr %>%
  group_by(bcr_name) %>%
  summarise(label_long = mean(range(long)), label_lat = mean(range(lat)), bcr_code = mean(bcr_code)) %>%
mutate(label_long = replace(label_long, bcr_name %in% c("Badlands And Prairies", "Shortgrass Prairie ", " Central Mixed Grass Prairie ", "Oaks And Prairies ", "West Gulf Coastal Plain/Ouachitas ", "Peninsular Florida  ", "Central Hardwoods","Eastern Tallgrass Prairie","Prairie Hardwood Transition","Boreal Hardwood Transition","Lower Great Lakes/St. Lawrence Plain","New England/Mid-Atlantic Coast"), 
                            values=c(-105.1,-103.1,-100.4,-97.1,-94.1,-82,-90,-89.7,-89.4,-83.8,-77,-73.7)),
       label_lat = replace(label_lat, bcr_name %in% c("Prairie Potholes", "Southeastern Coastal Plain ", "Piedmont  ", "Appalachian Mountains ","Mississippi Alluvial Valley "), 
                           values=c(48.5,32.5 ,36.3 ,38.3,33.9 )))
#create plot
map <- ggplot(tidy_bcr, aes(x = long, y = lat, group = group, fill = bcr_name)) +
  geom_polygon(color = "black", size = 0.1) +
  coord_equal() +
  theme_void()  +
  theme(plot.title = element_text(margin = margin(t = 40, b = -40)))
#plot wit textx
c <- map + geom_text(data = bcrLabel, mapping = aes(x=label_long, y = label_lat, label = bcr_code, group = NA), cex = 4, col = "white") + theme(legend.box.background = element_rect(colour = "black"), legend.position = "bottom")

####

#Begin the filtering of the eBird data
library(auk)
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select

# setup data directory
dir.create("data", showWarnings = FALSE)
#select the red-headed woodpecker txt file and the sampling text file
ebd <- auk_ebd("rehwoo.txt", 
               file_sampling = "sampling.txt")
#Filter the data
ebd_filters <- ebd %>% 
  auk_species("Red-headed Woodpecker") %>% 
  # june, use * to get data from any year
  auk_date(date = c("*-05-01", "*-08-31")) %>% 
  # restrict to the standard traveling and stationary count protocols
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
  auk_complete()

#create a file to store the values in
data_dir <- "data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}
#name the files
f_ebd <- file.path(data_dir, "rehwoo.txt")
f_sampling <- file.path(data_dir, "esampling.txt")

# only run if the files don't already exist
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
}

#fill the data with zeros
ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebd_zf <- ebd_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# additional filtering
ebd_zf_filtered <- ebd_zf %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last 10 years of data
    year >= 2010,
    # 10 or fewer observers
    number_observers <= 10)

ebird <- ebd_zf_filtered %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
write_csv(ebird, "data/rehwoo.csv", na = "")

####
#begin covariate analysis
library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

bcr <- read_sf("data/gis-data.gpkg", "bcr") %>% filter(bcr_code %in% c(11, 12, 13, 17, 18, 19,21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)) %>% st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                                                                                                                                                                        "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
ebird <- read_csv("data/rehwoo.csv")
begin_year <- format(min(ebird$observation_date), "%Y.01.01")
# end date for ebird data
end_year <- format(max(ebird$observation_date), "%Y.12.31")
tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                extent = bcr %>% st_buffer(dist = 10000), 
                begin = begin_year, end = end_year, 
                outDirPath = "data", job = "modis",
                MODISserverOrder = "LPDAAC") %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)

landcover <- list.files("data/modis", "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% 
  stack()
# label layers with year
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
landcover

max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()

neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
ebird_buff <- ebird %>% 
  distinct(year = format(observation_date, "%Y"),
           locality_id, latitude, longitude) %>% 
  # for 2019 use 2018 landcover data
  mutate(year_lc = if_else(as.integer(year) > max_lc_year, 
                           as.character(max_lc_year), year),
         year_lc = paste0("y", year_lc)) %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to modis projection
  st_transform(crs = projection(landcover)) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(year, locality_id, geometry))

calculate_pland <- function(yr, regions, lc) {
  locs <- st_set_geometry(regions, NULL)
  exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
    map(~ count(., landcover = value)) %>% 
    tibble(locs, data = .) %>% 
    unnest(data)
}
# iterate over all years extracting landcover for all checklists in each
lc_extract <- ebird_buff %>% 
  mutate(pland = map2(year_lc, data, calculate_pland, lc = landcover)) %>% 
  select(pland) %>% 
  unnest(cols = pland)

pland <- lc_extract %>% 
  # calculate proporiton
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

lc_names <- tibble(landcover = 0:15,
                   lc_name = c("pland_00_water", 
                               "pland_01_evergreen_needleleaf", 
                               "pland_02_evergreen_broadleaf", 
                               "pland_03_deciduous_needleleaf", 
                               "pland_04_deciduous_broadleaf", 
                               "pland_05_mixed_forest",
                               "pland_06_closed_shrubland", 
                               "pland_07_open_shrubland", 
                               "pland_08_woody_savanna", 
                               "pland_09_savanna", 
                               "pland_10_grassland", 
                               "pland_11_wetland", 
                               "pland_12_cropland", 
                               "pland_13_urban", 
                               "pland_14_mosiac", 
                               "pland_15_barren"))
pland <- pland %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s%>% 
pland <- pland %>% group_by(year) %>% mutate(nrow = row_number()) %>%
  pivot_wider(names_from = year, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% select(-nrow)

# save
write_csv(pland, "data/modis_pland_location-year.csv")
#aggregate the resolution to get 5x5 modis cells
agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- bcr %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r) %>% 
  # remove any empty cells at edges
  trim()
r <- writeRaster(r, filename = "data/prediction-surface.tif", overwrite = TRUE)
#extract points
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% 
  mutate(id = row_number())
#convert to polygon with neighbourhood radius
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
#extract values at every radius
lc_extract_pred <- landcover %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019)) %>% 
  tibble(id = r_cells$id, layer = r_cells$layer, data = .) %>% 
  unnest(data)
#remove nas
lc_extract_pred <- lc_extract_pred[complete.cases(lc_extract_pred),]
#calculat epland
pland_pred <- lc_extract_pred %>% 
  count(id,layer, y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019, landcover) %>% 
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(.))

lc_names <- tibble(y2010 = 0:15,y2011 = 0:15,y2012 = 0:15,y2013 = 0:15,y2014 = 0:15,y2015 = 0:15,y2016 = 0:15,y2017 = 0:15,y2018 = 0:15,y2019 = 0:15,
                   lc_name = c("pland_00_water", 
                               "pland_01_evergreen_needleleaf", 
                               "pland_02_evergreen_broadleaf", 
                               "pland_03_deciduous_needleleaf", 
                               "pland_04_deciduous_broadleaf", 
                               "pland_05_mixed_forest",
                               "pland_06_closed_shrubland", 
                               "pland_07_open_shrubland", 
                               "pland_08_woody_savanna", 
                               "pland_09_savanna", 
                               "pland_10_grassland", 
                               "pland_11_wetland", 
                               "pland_12_cropland", 
                               "pland_13_urban", 
                               "pland_14_mosiac", 
                               "pland_15_barren"))
#take the names of years
cols <- grep('y\\d+', names(pland_pred), value = TRUE)
#replace landcover values by names in lc_names
pland_pred[cols] <- Map(function(x, y) lc_names$lc_name[match(y, x)],
                        lc_names[cols], pland_pred[cols])


#time series
pland_pred.water <- pland_pred
pland_pred.deciduous_broadleaf <- pland_pred
pland_pred.mixed_forest <- pland_pred
pland_pred.woody_savanna <- pland_pred
pland_pred.savanna <- pland_pred
pland_pred.grassland <- pland_pred
pland_pred.wetland  <- pland_pred
pland_pred.cropland <- pland_pred
pland_pred.urban  <- pland_pred

#Replace those beloning to the landcover name by pland, and those not belonging to it the value 0
pland_pred.water[, 2:11] <- ifelse(pland_pred[, 2:11] == "pland_00_water", pland_pred$pland, 0)
pland_pred.deciduous_broadleaf[, 2:11] <-  ifelse(pland_pred[, 2:11] == "pland_04_deciduous_broadleaf", pland_pred$pland, 0)
pland_pred.mixed_forest[, 2:11] <-  ifelse(pland_pred[, 2:11] == "pland_05_mixed_forest", pland_pred$pland, 0)
pland_pred.woody_savanna[, 2:11] <-  ifelse(pland_pred[, 2:11] == "pland_08_woody_savanna", pland_pred$pland, 0)
pland_pred.savanna[, 2:11] <-  ifelse(pland_pred[, 2:11] == "pland_09_savanna", pland_pred$pland, 0)
pland_pred.grassland[, 2:11] <-  ifelse(pland_pred[, 2:11] == "pland_10_grassland", pland_pred$pland, 0)
pland_pred.wetland[, 2:11]  <-  ifelse(pland_pred[, 2:11] == "pland_11_wetland", pland_pred$pland, 0)
pland_pred.cropland[, 2:11] <- ifelse(pland_pred[, 2:11] == "pland_12_cropland", pland_pred$pland, 0)
pland_pred.urban[, 2:11]  <- ifelse(pland_pred[, 2:11] == "pland_13_urban", pland_pred$pland, 0)

#filter
pland_pred.water <- pland_pred.water[, -12]
pland_pred.deciduous_broadleaf <- pland_pred.deciduous_broadleaf[, -12]
pland_pred.mixed_forest <- pland_pred.mixed_forest[, -12]
pland_pred.woody_savanna <- pland_pred.woody_savanna[, -12]
pland_pred.savanna <- pland_pred.savanna[, -12]
pland_pred.grassland <- pland_pred.grassland[,-12]
pland_pred.wetland <- pland_pred.wetland[, -12]
pland_pred.cropland <- pland_pred.cropland[, -12]
pland_pred.urban <- pland_pred.urban[, -12]

#join the landcover dataframe with the coordinates from the bcr range
pland_coords.water <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.water, by = "id")

pland_coords.deciduous_broadleaf <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.deciduous_broadleaf, by = "id")

pland_coords.mixed_forest <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.mixed_forest, by = "id")

pland_coords.woody_savanna <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.woody_savanna, by = "id")

pland_coords.savanna <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.savanna, by = "id")

pland_coords.grassland <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.grassland, by = "id")

pland_coords.wetland <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.wetland, by = "id")

pland_coords.cropland <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.cropland, by = "id")

pland_coords.urban <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred.urban, by = "id")



#create landcover maps
forest_cover.water <- pland_coords.water %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.deciduous_broadleaf <- pland_coords.deciduous_broadleaf %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.mixed_forest <- pland_coords.mixed_forest %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()


forest_cover.woody_savanna <- pland_coords.woody_savanna %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.savanna <- pland_coords.savanna %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.grassland <- pland_coords.grassland %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.wetland <- pland_coords.wetland %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.cropland <- pland_coords.cropland %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

forest_cover.urban <- pland_coords.urban %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()



#write landcover maps
writeRaster(forest_cover.water, "forest_cover.water.tif")
writeRaster(forest_cover.deciduous_broadleaf, "forest_cover.deciduous_broadleaf.tif")
writeRaster(forest_cover.mixed_forest, "forest_cover.mixed_forest.tif")
writeRaster(forest_cover.woody_savanna, "forest_cover.woody_savanna.tif")
writeRaster(forest_cover.savanna, "forest_cover.savanna.tif")
writeRaster(forest_cover.grassland, "forest_cover.grassland.tif")
writeRaster(forest_cover.wetland, "forest_cover.wetland.tif")
writeRaster(forest_cover.cropland, "forest_cover.cropland.tif")
writeRaster(forest_cover.urban, "forest_cover.urban.tif")


#Create elevation and soil maps
#download the file; replace elevation with soil to get soil and repeat the code for soil
f_dem <- "elevation_01_05_1km_uint16.tif"
if (!file.exists(file.path("data", f_dem))) {
  download.file(paste0("https://data.earthenv.org/texture/", f_dem),
                file.path("data", f_dem), mode = "wb")
}

elev <- raster("data/elevation_01_05_1km_uint16.tif")
# crop, buffer bcr by 10 km to provide a little wiggly room
elev <- bcr %>% 
  st_buffer(dist = 10000) %>% 
  st_transform(crs = projection(elev)) %>% 
  crop(elev, .) %>% 
  projectRaster(crs = projection(landcover))

ebird_buff_noyear <- ebird %>% 
  distinct(locality_id, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(elev)) %>% 
  st_buffer(dist = neighborhood_radius)
# extract elevation values and calculate median and sd
locs <- st_set_geometry(ebird_buff_noyear, NULL) %>% 
  mutate(id = row_number())
elev_checklists <- exact_extract(elev, ebird_buff_noyear, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = exp(mean(log(.$value),na.rm = TRUE)),
                   elevation_sd = exp(sd(log(.$value),na.rm = TRUE )))) %>% 
  # join to lookup table to get locality_id
  bind_cols(locs, .)

elev_pred <- exact_extract(elev, r_cells, progress = FALSE) %>% 
  map_dfr(~ tibble(elevation_median = exp(mean(log(.$value), na.rm = TRUE)),
                   elevation_sd = exp(sd(log(.$value), na.rm = TRUE)))) %>% 
  # join to lookup table to get locality_id
  bind_cols(st_drop_geometry(r_cells), .)


pland_elev_checklist <- inner_join(pland, elev_checklists, by = "locality_id")
write_csv(pland_elev_checklist, "data/year.pland-elevation_location-year.csv")
# prediction surface covariates
pland_elev_pred <- inner_join(pland_coords.water, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.water-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.deciduous_broadleaf, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.decidous_broadleaf-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.mixed_fored, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.mixed_forest-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.woody_savanna, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.woody_savanna-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.savanna, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.savanna-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.cropland, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.cropland-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.grassland, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.grassland-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.wetland, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.wetland-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords.urban, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/l.urban-elev_prediction-surface.csv")
pland_elev_pred <- inner_join(pland_coords10, elev_pred, by = "id")
write_csv(pland_elev_pred, "data/pland2019-elev_prediction-surface.csv")

#create encounter rate models; repeat this for each different landcover class
library(sf)
library(raster)
library(dggridR)
library(lubridate)
library(ranger)
library(scam)
library(PresenceAbsence)
library(verification)
library(ebirdst)
library(fields)
library(gridExtra)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection
# set random number seed to insure fully repeatable results
set.seed(1)
# setup output directory for saved results
if (!dir.exists("output")) {
  dir.create("output")
}
######################################################################################################################################
######################################################################################################################################
ebird <- read_csv("data/rehwoo.csv") %>% 
  # year required to join to habitat data
  mutate(year = year(observation_date))
# modis habitat covariates


k <- read_csv("data/year.pland-elev_location-year.csv") %>% arrange(id)
l <- read_csv("data/year.pland-slope_location-year.csv") %>% arrange(id)

habitat <- cbind(l, k[,15:16])

cols <- grep('y\\d+', names(habitat))
habitat$elevation_median[rowSums(habitat[cols] == 0) == length(cols)] <- 0
habitat$elevation_sd[rowSums(habitat[cols] == 0) == length(cols)] <- 0
habitat$slope_median[rowSums(habitat[cols] == 0) == length(cols)] <- 0
habitat$slope_sd[rowSums(habitat[cols] == 0) == length(cols)] <- 0

[1] "pland_00_water"                "pland_04_deciduous_broadleaf"  "pland_05_mixed_forest"         "pland_08_woody_savanna"       
[5] "pland_09_savanna"              "pland_11_wetland"              "pland_12_cropland"             "pland_14_mosiac"              
[9] "pland_13_urban"                "pland_10_grassland"            "pland_15_barren"               "pland_01_evergreen_needleleaf"
[13] "pland_02_evergreen_broadleaf"  "pland_03_deciduous_needleleaf" "pland_07_open_shrubland"       "pland_06_closed_shrubland"   

habitat1 <- habitat[habitat$lc_name %in% "pland_00_water",]


habitat1$year <- as.integer(habitat1$year)
rm(k, l)
# combine ebird and habitat data

ebird_habitat <- inner_join(ebird, habitat1, by = c("locality_id", "year"))
# prediction surface


j <- read_csv("data/l.water-elev_prediction-surface.csv") %>% arrange(id)
h <- read_csv("data/l.water-slope_prediction-surface.csv") %>% arrange(id)

pred_surface_2010 <- cbind(j, h[,15:16])

cols <- grep('y\\d+', names(pred_surface_2010))
pred_surface_2010$elevation_median[rowSums(pred_surface_2010[cols] == 0) == length(cols)] <- 0
pred_surface_2010$elevation_sd[rowSums(pred_surface_2010[cols] == 0) == length(cols)] <- 0
pred_surface_2010$slope_median[rowSums(pred_surface_2010[cols] == 0) == length(cols)] <- 0
pred_surface_2010$slope_sd[rowSums(pred_surface_2010[cols] == 0) == length(cols)] <- 0



rm(j, h)
# latest year of landcover data
max_lc_year <- max(pred_surface$year)
r <- raster("data/layer.prediction-surface.tif")
# load gis data for making maps
map_proj <- st_crs("ESRI:102003")
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
bcr <- read_sf("data/gis-data.gpkg", "bcr") %>% filter(bcr_code %in% c(11, 12, 13, 17, 18, 19,21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)) %>%
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf("data/gis-data.gpkg", "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf("data/gis-data.gpkg", "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

bb <- st_bbox(c(xmin = -0.1, xmax = 0.1, ymin = -0.1, ymax = 0.1), 
              crs = 4326) %>% 
  st_as_sfc() %>% 
  st_sf()
# random points
pts <- st_sample(bb, 500) %>% 
  st_sf(as.data.frame(st_coordinates(.)), geometry = .) %>% 
  rename(lat = Y, lon = X)
# contruct a hexagonal grid with ~ 5 km between cells
dggs <- dgconstruct(spacing = 5)
# for each point, get the grid cell
pts$cell <- dgGEO_to_SEQNUM(dggs, pts$lon, pts$lat)$seqnum
# sample one checklist per grid cell
pts_ss <- pts %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup()
# generate polygons for the grid cells
hexagons <- dgcellstogrid(dggs, unique(pts$cell), frame = FALSE) %>% 
  st_as_sf()
ggplot() +
  geom_sf(data = hexagons) +
  geom_sf(data = pts, size = 0.5) +
  geom_sf(data = pts_ss, col = "red") +
  theme_bw()
set.seed(1)

dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id and week number for each checklist
checklist_cell <- ebird_habitat %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
         year = year(observation_date),
         week = week(observation_date))
# sample one checklist per grid cell per week
# sample detection/non-detection independently 
ebird_ss <- checklist_cell %>% 
  group_by(species_observed, year, week, cell) %>% 
  sample_n(size = 1) %>% 
  ungroup()

nrow(ebird_habitat)
count(ebird_habitat, species_observed) %>% 
  mutate(percent = n / sum(n))
# after sampling
nrow(ebird_ss)
count(ebird_ss, species_observed) %>% 
  mutate(percent = n / sum(n))

pct_before <- count(ebird_habitat, species_observed) %>% 
  mutate(proportion = n / sum(n)) %>% 
  filter(species_observed) %>% 
  pull(proportion) %>% 
  round(3) %>%
  scales::percent()
pct_after <- count(ebird_ss, species_observed) %>% 
  mutate(proportion = n / sum(n)) %>% 
  filter(species_observed) %>% 
  pull(proportion) %>% 
  round(3) %>%
  scales::percent()

ebird_split <- ebird_ss %>% 
  # select only the columns to be used in the model
  select(species_observed, 
         time_observations_started, duration_minutes,
         effort_distance_km, number_observers,
         starts_with("y"),
         starts_with("elevation_"),
         starts_with("slope_"))%>% 
  drop_na()
# split 80/20
ebird_split <- ebird_split %>% 
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))
map_int(ebird_split, nrow)

detection_freq <- mean(ebird_split$train$species_observed)

ebird_split$train$species_observed <- factor(ebird_split$train$species_observed)
# grow random forest

ebird_split$train$year <- NULL

rf <- ranger(formula =  species_observed ~ ., 
             data = ebird_split$train,
             importance = "impurity",
             probability = TRUE,
             replace = TRUE, 
             sample.fraction = c(detection_freq, detection_freq),)

occ_pred <- rf$predictions[, 2]
# convert the observered response back to a numeric value from factor
occ_obs <- ebird_split$train$species_observed %>% 
  as.logical() %>% 
  as.integer()
rf_pred_train <- tibble(obs = occ_obs, pred = occ_pred) %>% 
  drop_na()
# fit calibration model
calibration_model <- scam(obs ~ s(pred, k = 5, bs = "mpi"), 
                          gamma = 1.4,
                          data = rf_pred_train)

search_hours <- ebird_split$train %>% 
  mutate(hour = floor(time_observations_started)) %>%
  count(hour) %>% 
  mutate(pct = n / sum(n)) %>% 
  filter(pct >= 0.01)
# constrained peak time
t_peak <- pd_time %>% 
  filter(floor(time_observations_started) %in% search_hours$hour) %>% 
  top_n(1, wt = desc(time_observations_started)) %>% 
  pull(time_observations_started)
t_peak

human_time <- str_glue("{h}:{m} {ap}", 
                       h = floor(t_peak),
                       m = str_pad(round((t_peak %% 1) * 60), 2, pad = "0"),
                       ap = ifelse(t_peak > 12, "PM", "AM"))
#prepare checklist data for prediction using missRanger to fill missing values
#library(missRanger)
#pred_surface <- missRanger(pred_surface, pmm.k = 3, num.trees = 100)
pred_surface_eff_2010 <- pred_surface_2010 %>% 
  mutate(time_observations_started = t_peak,
         duration_minutes = 60,
         effort_distance_km = 1,
         number_observers = 1)


# predict
pred_rf_2010 <- predict(rf, data = pred_surface_eff_2010, type = "response")
pred_rf_2010 <- pred_rf_2010$predictions[, 2]
# apply calibration models
pred_rf_cal_2010 <- predict(calibration_model, 
                            data.frame(pred = pred_rf_2010), 
                            type = "response")
# add to prediction surface
pred_er_2010 <- bind_cols(pred_surface_2010, encounter_rate = pred_rf_2010_cal_2010) %>% 
  select(latitude, longitude, encounter_rate, id, layer) %>% 
  mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))
#select the names of the columns
cols <- grep('y\\d+', names(pred_surface_2010))
#add those columns to the dataframe with predictions
pred_er_2010 <- pred_er_2010 %>% mutate(pred_surface_2010[, c(cols)])
cols <- grep('y\\d+', names(pred_er_2010))
#select only those predictions with landcover values otherwise turn them to 0
pred_er_2010$encounter_rate[rowSums(pred_er_2010[cols] == 0) == length(cols)] <- 0

pred_er_2010 <- pred_er_2010[, -c(cols)]

r_pred_2010 <- pred_er_2010 %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred_2010 <- r_pred_2010[[-1]]
# save the raster

writeRaster(r_pred_2010, "output/p.barren.tif")

##gather both increasing and decreasing landcover changes relative to regions and encounter rate
library(raster)
library(sf)


habs <- stack("l.forest_cover.water.tif")
rf <- stack("output/p.water.tif")

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


write.csv(c_change, "t.water.csv")


#filter into groups and then select for changes in pland
library(tidyverse)
t_cropland <- read.csv("q.cropland.csv")
t_deciduous_broadleaf <- read.csv("q.deciduous_broadleaf.csv")
t_grassland <- read.csv("q.grassland.csv")
t_mixed_forest <- read.csv("q.mixed_forest.csv")
t_mosaic <- read.csv("q.mosaic.csv")
t_urban <- read.csv("q.urban.csv")
t_wetland <- read.csv("q.wetland.csv")
t_savanna <- read.csv("q.savanna.csv")
t_woody_savanna <- read.csv("q.woody_savanna.csv")
t_water <- read.csv("q.water.csv")

names(t_cropland)[5] <- "cropland"
names(t_deciduous_broadleaf)[5] <- "deciduous_broadleaf"
names(t_grassland)[5] <- "grassland"
names(t_mixed_forest)[5] <- "mixed_forest"
names(t_mosaic)[5] <- "mosaic"
names(t_urban)[5] <- "urban"
names(t_wetland)[5] <- "wetland"
names(t_savanna)[5] <- "savanna"
names(t_woody_savanna)[5] <- "woody_savanna"

t_cropland <- t_cropland[, -1]
t_deciduous_broadleaf <- t_deciduous_broadleaf[, -1]
t_grassland <- t_grassland[, -1]
t_mixed_forest <- t_mixed_forest[, -1]
t_mosaic <- t_mosaic[, -1]
t_urban <- t_urban[, -1]
t_wetland <-t_wetland [, -1]
t_savanna <- t_savanna[, -1]
t_woody_savanna <- t_woody_savanna[, -1]
t_water <- t_water[, -1]


t_cropland <- t_cropland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_deciduous_broadleaf <- t_deciduous_broadleaf%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_grassland <- t_grassland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_mixed_forest <- t_mixed_forest%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_urban <- t_urban%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_wetland <- t_wetland%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_savanna <- t_savanna%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_woody_savanna <-t_woody_savanna%>% pivot_longer(-c(1, 2, 3, 5, 6, 7))
t_water <- t_water %>% pivot_longer(-c(1, 2, 3, 5, 6, 7))


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

Woodland <- rbind(t_deciduous_broadleaf, t_mixed_forest)
Woodland <- Woodland[, -2]
Woodland <- Woodland[, c(4, 5, 3, 2, 6, 7, 1)]
colnames(Woodland) <- c("X", "Y", "year", "layer", "landcover", "pland_change", "encounter")
Woodland <- Woodland %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

Woodland <- Woodland[complete.cases(Woodland),]
Woodland.2 <- Woodland[Woodland$pland_change <=0,]
Woodland.3 <- Woodland[Woodland$pland_change >=0,]


#landcover
Woodland.2$landcover <- as.factor(Woodland.2$landcover)
Woodland.2$landcover <- as.numeric(Woodland.2$landcover)
Woodland.3$landcover <- as.factor(Woodland.3$landcover)
Woodland.3$landcover <- as.numeric(Woodland.3$landcover)

#ly_name
Woodland.2$ly_name <- as.factor(Woodland.2$ly_name)
Woodland.2$ly_name <- as.numeric(Woodland.2$ly_name)
Woodland.3$ly_name <- as.factor(Woodland.3$ly_name)
Woodland.3$ly_name <- as.numeric(Woodland.3$ly_name)

#correlation test
cor.test(Woodland.3$ly_name, Woodland.3$encounter)
cor.test(Woodland.3$landcover, Woodland.3$encounter)
cor.test(Woodland.2$ly_name, Woodland.2$encounter)
cor.test(Woodland.2$landcover, Woodland.2$encounter)

#
savannas <- rbind(t_woody_savanna, t_savanna)
savannas <- savannas[, -2]
savannas <- savannas[, c(4, 5, 3, 2, 6, 7, 1)]
colnames(savannas) <- c("X", "Y", "year", "layer", "landcover", "pland_change", "encounter")
savannas <- savannas %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

savannas <- savannas[complete.cases(savannas),]
savannas.2 <- savannas[savannas$pland_change <=0,]
savannas.3 <- savannas[savannas$pland_change >=0,]

#landcover
savannas.2$landcover <- as.factor(savannas.2$landcover)
savannas.2$landcover <- as.numeric(savannas.2$landcover)
savannas.3$landcover <- as.factor(savannas.3$landcover)
savannas.3$landcover <- as.numeric(savannas.3$landcover)

#ly_name
savannas.2$ly_name <- as.factor(savannas.2$ly_name)
savannas.2$ly_name <- as.numeric(savannas.2$ly_name)
savannas.3$ly_name <- as.factor(savannas.3$ly_name)
savannas.3$ly_name <- as.numeric(savannas.3$ly_name)

#correlation test
cor.test(savannas.3$ly_name, savannas.3$encounter)
cor.test(savannas.3$landcover, savannas.3$encounter)
cor.test(savannas.2$ly_name, savannas.2$encounter)
cor.test(savannas.2$landcover, savannas.2$encounter)

agricultural <- rbind(t_cropland, t_grassland)
agricultural <- agricultural[, -2]
agricultural <- agricultural[, c(4, 5, 3, 2, 6, 7, 1)]
colnames(agricultural) <- c("X", "Y", "year", "layer", "landcover", "pland_change", "encounter")
agricultural <- agricultural %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

agricultural <- agricultural[complete.cases(agricultural),]
agricultural.2 <- agricultural[agricultural$pland_change <=0,]
agricultural.3 <- agricultural[agricultural$pland_change >=0,]

#landcover
agricultural.2$landcover <- as.factor(agricultural.2$landcover)
agricultural.2$landcover <- as.numeric(agricultural.2$landcover)
agricultural.3$landcover <- as.factor(agricultural.3$landcover)
agricultural.3$landcover <- as.numeric(agricultural.3$landcover)

#ly_name
agricultural.2$ly_name <- as.factor(agricultural.2$ly_name)
agricultural.2$ly_name <- as.numeric(agricultural.2$ly_name)
agricultural.3$ly_name <- as.factor(agricultural.3$ly_name)
agricultural.3$ly_name <- as.numeric(agricultural.3$ly_name)

#correlation test
cor.test(agricultural.3$ly_name, agricultural.3$encounter)
cor.test(agricultural.3$landcover, agricultural.3$encounter)
cor.test(agricultural.2$ly_name, agricultural.2$encounter)
cor.test(agricultural.2$landcover, agricultural.2$encounter)

urban <- rbind(t_urban)
urban <- urban[, -2]
urban <- urban[, c(4, 5, 3, 2, 6, 7, 1)]
colnames(urban) <- c("X", "Y", "year", "layer", "landcover", "pland_change", "encounter")
urban <- urban %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

urban <- urban[complete.cases(urban),]
urban.2 <- urban[urban$pland_change <=0,]
urban.3 <- urban[urban$pland_change >=0,]

#landcover
urban.2$landcover <- as.factor(urban.2$landcover)
urban.2$landcover <- as.numeric(urban.2$landcover)
urban.3$landcover <- as.factor(urban.3$landcover)
urban.3$landcover <- as.numeric(urban.3$landcover)

#ly_name
urban.2$ly_name <- as.factor(urban.2$ly_name)
urban.2$ly_name <- as.numeric(urban.2$ly_name)
urban.3$ly_name <- as.factor(urban.3$ly_name)
urban.3$ly_name <- as.numeric(urban.3$ly_name)

#correlation test
cor.test(urban.3$ly_name, urban.3$encounter)
cor.test(urban.3$landcover, urban.3$encounter)
cor.test(urban.2$ly_name, urban.2$encounter)
cor.test(urban.2$landcover, urban.2$encounter)
cor.test(urban.3$pland_change, urban.3$encounter)


water_bodies <- rbind(t_water, t_wetland)
water_bodies <- water_bodies[, -2]
water_bodies <- water_bodies[, c(4, 5, 3, 2, 6, 7, 1)]
colnames(water_bodies) <- c("X", "Y", "year", "layer", "landcover", "pland_change", "encounter")
water_bodies <- water_bodies %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

water_bodies <- water_bodies[complete.cases(water_bodies),]
water_bodies.2 <- water_bodies[water_bodies$pland_change <=0,]
water_bodies.3 <- water_bodies[water_bodies$pland_change >=0,]

#landcover
water_bodies.2$landcover <- as.factor(water_bodies.2$landcover)
water_bodies.2$landcover <- as.numeric(water_bodies.2$landcover)
water_bodies.3$landcover <- as.factor(water_bodies.3$landcover)
water_bodies.3$landcover <- as.numeric(water_bodies.3$landcover)

#ly_name
water_bodies.3$ly_name <- as.factor(water_bodies.3$ly_name)
water_bodies.3$ly_name <- as.numeric(water_bodies.3$ly_name)
water_bodies.3$ly_name <- as.factor(water_bodies.3$ly_name)
water_bodies.3$ly_name <- as.numeric(water_bodies.3$ly_name)

#correlation test
cor.test(water_bodies.3$ly_name, water_bodies.3$encounter)
cor.test(water_bodies.3$landcover, water_bodies.3$encounter)
cor.test(water_bodies.2$ly_name, water_bodies.2$encounter)
cor.test(water_bodies.2$landcover, water_bodies.2$encounter)

#create a big dataframe containing both decreases and increases
test.all <- rbind(t_cropland,t_deciduous_broadleaf,t_grassland,t_mixed_forest,t_urban,t_wetland,t_savanna,t_woody_savanna,t_water)
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
#split this dataframe
test.all <- test.all[complete.cases(test.all),]
test.all.2 <- test.all[test.all$pland_change <=0,]
test.all.3 <- test.all[test.all$pland_change >=0,]

round(max(test.all.3$encounter), 2)/10

#calculate ranges of probability
myIntervals <- c("0 - 0.1", "0.1 - 0.2 ", "0.2  - 0.3", "0.3  - 0.4","0.4  - 0.5","0.5  - 0.6","0.6  - 0.6","0.7  - 0.8","0.8 - 0.8","0.9  - 1")
test.all.3$encounter_range<- myIntervals[findInterval(test.all.3$encounter, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))]

tt3.l <- test.all.3 %>% count(encounter_range, landcover, ly_name) %>% group_by(ly_name)

#plot these ranges but first make sure the colours match the number of bcr regions
library(RColorBrewer)
brks <- seq(0, 1, 0.057)
cols <- colorRampPalette(c("red4","gray95","orange","gold","lightskyblue","dodgerblue","blue3"))(length(brks)-1) #this is the colour ramp - change the colours as you like, bu make sure they are symetrical around the "gray95", which represents the zero point!

forestlossline <- ggplot(tt3.l, aes(x=factor(encounter_range), y=n, group=ly_name, label=round(n, 1))) + 
  geom_bar(aes(fill=ly_name), stat='identity') +
  scale_colour_manual(values=cols) +
  scale_fill_manual(values=cols)+
  ylab("Counts of Encounters in BCR regions for Landcover classes") +
  xlab("Encounter rate Ranges") +
  facet_wrap(~landcover, scales="free") +
  theme_classic() +
  scale_x_discrete(label = function(x) stringr::str_trunc(x, 12))+
  theme(
    axis.title.x=element_text(face="bold", size=13, family="TT Times New Roman"),
    axis.title.y=element_text(face="bold", size=13, family="TT Times New Roman"),
    axis.text.x=element_text(face="bold", size =10, family="TT Times New Roman", angle=90), 
    axis.text.y=element_text(face="bold",  size =10, family="TT Times New Roman"),
    legend.title=element_blank(),
    legend.text=element_text(color="black", size =10, face="bold", family="TT Times New Roman"),
    legend.justification=c(1.2,1),
    plot.title=element_text(face="bold", size = 18, hjust=0.5, colour = "black"),
    axis.line=element_blank(),
    legend.key.height=unit(.1, "cm")
  )


#calculate the stepforward regression using the stepfor function built
library(stepfor)

tbest1 <- stepfor(savannas.2$encounter, savannas.2[, -c(1, 4, 2, 6)], alpha=0.2)



#get points from landcover changes
#calculate a buffer radius around each point of observation
neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
ebird_buff <- Woodland.3[Woodland.3$landcover %in% "deciduous_broadleaf",] %>% 
  distinct(year, X, Y,encounter) %>% 
  # for 2019 use 2018 landcover data
  mutate(year_lc = if_else(as.integer(year) > max_lc_year, 
                           as.character(max_lc_year), as.character(year)),
         year_lc = paste0("y", year_lc)) %>% 
  # convert to spatial features
  st_as_sf(coords = c("X", "Y")) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(year, encounter, geometry))

#raster
agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
########################################################

#conver to dataframe with geometry points
ebird_buff <- ebird_buff %>% arrange(year_lc)
#2010
bird_buff <- NULL
for(i in 1:9){
  bird_buff[[i]] <- ebird_buff$data[i] %>%
    bind_cols %>%
    st_cast(to = "POINT") %>%
    dplyr::mutate(
      X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
      Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
    ) %>%
    sf::st_drop_geometry()}

#rbind the lists
bird_buff <- do.call(rbind.data.frame, bird_buff)


bird_buff$year <- as.integer(bird_buff$year)
bird_buff$X <- as.numeric(bird_buff$X)
bird_buff$Y <- as.numeric(bird_buff$Y)

#convert into raster
bird.raster <- bird_buff %>% st_as_sf(coords = c("X", "Y")) %>% rasterize(r) 
map_proj <- st_crs("ESRI:102003")
r_bird<- projectRaster(bird.raster, crs = map_proj$proj4string, method = "ngb")

#extract points and create a spatial points
rpoints <- rasterToPoints(r_bird) %>% data.frame()
xy <- SpatialPointsDataFrame(coords = rpoints[, 1:2], rpoints, proj4string = CRS("+proj=aea +lat_0=37.5 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))


#Load the maps into R and filter the coordinates points as you load them also

#
map_proj <- st_crs("ESRI:102003")
##########################
p.water <- stack("output/p.water.tif")
p.water <- projectRaster(p.water[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.evergreen_needleleaf <- stack("output/p.evergreen_needleleaf.tif")
p.evergreen_needleleaf <- projectRaster(p.evergreen_needleleaf[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.evergreen_broadleaf <- stack("output/p.evergreen_broadleaf.tif")
p.evergreen_broadleaf <- projectRaster(p.evergreen_broadleaf[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.deciduous_broadleaf <- stack("output/p.deciduous_broadleaf.tif")
p.deciduous_broadleaf <- projectRaster(p.deciduous_broadleaf[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.deciduous_needleleaf <- stack("output/p.deciduous_needleleaf.tif")
p.deciduous_needleleaf <- projectRaster(p.deciduous_needleleaf[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.mixed_forest <- stack("output/p.mixed_forest.tif")
p.mixed_forest <- projectRaster(p.mixed_forest[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.closed_shrubland <- stack("output/p.closed_shrubland.tif")
p.closed_shrubland <- projectRaster(p.closed_shrubland[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.open_shrubland <- stack("output/p.open_shrubland.tif")
p.open_shrubland <- projectRaster(p.open_shrubland[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.woody_savanna <- stack("output/p.woody_savanna.tif")
p.woody_savanna <- projectRaster(p.woody_savanna[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.savanna <- stack("output/p.savanna.tif")
p.savanna <- projectRaster(p.savanna[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.grassland <- stack("output/p.grassland.tif")
p.grassland <- projectRaster(p.grassland[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.cropland <- stack("output/p.cropland.tif")
p.cropland <- projectRaster(p.cropland[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.urban <- stack("output/p.urban.tif")
p.urban <- projectRaster(p.urban[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.wetland <- stack("output/p.wetland.tif")
p.wetland <- projectRaster(p.wetland[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.mosaic <- stack("output/p.mosaic.tif")
p.mosaic <- projectRaster(p.mosaic[[1]], crs = map_proj$proj4string, method = "ngb")
##########################
p.barren <- stack("output/p.barren.tif")
p.barren <- projectRaster(p.barren[[1]], crs = map_proj$proj4string, method = "ngb")


######################################################################################
xy_barren <- xy_barren[,-c(1, 7:9)]
xy_water <- xy_water[,-c(1, 7:9)]
xy_evergreen_needleleaf <- xy_evergreen_needleleaf[,-c(1, 7:9)]
xy_evergreen_broadleaf <- xy_evergreen_broadleaf[,-c(1, 7:9)]
xy_deciduous_broadleaf <- xy_deciduous_broadleaf[,-c(1, 7:9)]
xy_mixed_forest <- xy_mixed_forest[,-c(1, 7:9)]
xy_open_shrubland <- xy_open_shrubland[,-c(1, 7:9)]
xy_savanna <- xy_savanna[,-c(1, 7:9)]
xy_grassland <- xy_grassland[,-c(1, 7:9)]
xy_wetland <- xy_wetland[,-c(1, 7:9)]
xy_cropland <- xy_cropland[,-c(1, 7:9)]
xy_mosaic <- xy_mosaic[,-c(1, 7:9)]
xy_woody_savanna <- xy_woody_savanna[,-c(1, 7:9)]
ixy_barren <- ixy_barren[,-c(1, 7:9)]
ixy_water <- ixy_water[,-c(1, 7:9)]
ixy_evergreen_needleleaf <- ixy_evergreen_needleleaf[,-c(1, 7:9)]
ixy_evergreen_broadleaf <- ixy_evergreen_broadleaf[,-c(1, 7:9)]
ixy_deciduous_needleleaf <- ixy_deciduous_needleleaf[,-c(1, 7:9)]
ixy_deciduous_broadleaf <- ixy_deciduous_broadleaf[,-c(1, 7:9)]
ixy_mixed_forest <- ixy_mixed_forest[,-c(1, 7:9)]
ixy_closed_shrubland <- ixy_closed_shrubland[,-c(1, 7:9)]
ixy_open_shrubland <- ixy_open_shrubland[,-c(1, 7:9)]
ixy_woody_savanna <- ixy_woody_savanna[,-c(1, 7:9)]
ixy_savanna <- ixy_savanna[,-c(1, 7:9)]
ixy_grassland <- ixy_grassland[,-c(1, 7:9)]
ixy_wetland <- ixy_wetland[,-c(1, 7:9)]
ixy_cropland <- ixy_cropland[,-c(1, 7:9)]
ixy_mosaic <- ixy_mosaic[,-c(1, 7:9)]
ixy_urban <- ixy_urban[,-c(1, 7:9)]

#create maps
windows(10, 8)

par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(bcr, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)

# encounter rate
r_max <- ceiling(10 * cellStats(p.water, max)) / 10
brks <- seq(0, r_max, by = 0.006)

cols <- colorRampPalette(c("red4","gray95","orange","gold","lightskyblue","dodgerblue","blue3"))(length(brks)-1) #this is the colour ramp - change the colours as you like, bu make sure they are symetrical around the "gray95", which represents the zero point!

plot(p.water, 
     col = cols, breaks = brks, 
     maxpixels = ncell(p.water),
     legend = FALSE, add = TRUE)
points(xy_water[, 1:2], add=TRUE, pch=".", cex=2, col="darkred")
points(ixy_water[, 1:2], add=TRUE, pch=1, cex=1, col="blue")

# borders
plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)

box()

brks <- seq(0, 1, by = 0.05)
cols <- colorRampPalette(c("red4","gray95","orange","gold","lightskyblue","dodgerblue","blue3"))(length(brks)-1) #this is the colour ramp - change the colours as you like, bu make sure they are symetrical around the "gray95", which represents the zero point!

lbl_brks <- seq(0, 1, by = 0.05)
# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Red-headed Woodpecker encounter rate"
image.plot(zlim = range(brks), legend.only = TRUE, 
           col = cols, breaks = brks,
           smallplot = c(0.25, 0.75, 0.06, 0.09),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))


