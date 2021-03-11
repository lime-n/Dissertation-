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

ebird$species_observed <- ifelse(ebird$species_observed == TRUE, 1,0)

#calculate a buffer radius around each point of observation
neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
ebird_buff <- ebird %>% 
  distinct(year = format(observation_date, "%Y"), latitude, longitude,species_observed) %>% 
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
  nest(data = c(year, species_observed, geometry))

#raster
agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
########################################################

ebird_buff <- ebird_buff %>% arrange(year_lc)
#2010
bird_buff <- ebird_buff$data[1] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff$year <- as.integer(bird_buff$year)
bird_buff$X <- as.numeric(bird_buff$X)
bird_buff$Y <- as.numeric(bird_buff$Y)


bird.raster <- bird_buff %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% rasterize(r) 

rpoints <- rasterToPoints(bird.raster) %>% data.frame()
#####################################################################
#2011
bird_buff.2011 <- ebird_buff$data[2] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2011$year <- as.integer(bird_buff.2011$year)
bird_buff.2011$X <- as.numeric(bird_buff.2011$X)
bird_buff.2011$Y <- as.numeric(bird_buff.2011$Y)

bird.raster.2011 <- bird_buff.2011 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>%  rasterize(r)

rpoints.2011 <- rasterToPoints(bird.raster.2011) %>% data.frame()

rm(bird_buff, bird_buff.2011, bird.raster, bird.raster.2011)
#####################################################################
#2012
bird_buff.2012 <- ebird_buff$data[3] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2012$year <- as.integer(bird_buff.2012$year)
bird_buff.2012$X <- as.numeric(bird_buff.2012$X)
bird_buff.2012$Y <- as.numeric(bird_buff.2012$Y)

bird.raster.2012 <- bird_buff.2012 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>%  rasterize(r)

rpoints.2012 <- rasterToPoints(bird.raster.2012) %>% data.frame()

rm(bird_buff.2012, bird.raster.2012)


rpoints.all <- rbind(rpoints, rpoints.2011, rpoints.2012, rpoints.2013, rpoints.2014, rpoints.2015, rpoints.2016, rpoints.2017, rpoints.2018, rpoints.2019)


####################################################################
#2013
bird_buff.2013 <- ebird_buff$data[4] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2013$year <- as.integer(bird_buff.2013$year)
bird_buff.2013$X <- as.numeric(bird_buff.2013$X)
bird_buff.2013$Y <- as.numeric(bird_buff.2013$Y)

bird.raster.2013 <- bird_buff.2013 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>%  rasterize(r)

rpoints.2013 <- rasterToPoints(bird.raster.2013) %>% data.frame()

rm(bird_buff.2013, bird.raster.2013)
###################################################################
#2014
bird_buff.2014 <- ebird_buff$data[5] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2014$year <- as.integer(bird_buff.2014$year)
bird_buff.2014$X <- as.numeric(bird_buff.2014$X)
bird_buff.2014$Y <- as.numeric(bird_buff.2014$Y)

bird.raster.2014 <- bird_buff.2014 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% rasterize(r)

rpoints.2014 <- rasterToPoints(bird.raster.2014) %>% data.frame()

rm(bird_buff.2014, bird.raster.2014)
##################################################################
#2015
bird_buff.2015 <- ebird_buff$data[6] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2015$year <- as.integer(bird_buff.2015$year)
bird_buff.2015$X <- as.numeric(bird_buff.2015$X)
bird_buff.2015$Y <- as.numeric(bird_buff.2015$Y)

bird.raster.2015 <- bird_buff.2015 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>%  rasterize(r)

rpoints.2015 <- rasterToPoints(bird.raster.2015) %>% data.frame()

rm(bird_buff.2015, bird.raster.2015)
#################################################################
#2016
bird_buff.2016 <- ebird_buff$data[7] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2016$year <- as.integer(bird_buff.2016$year)
bird_buff.2016$X <- as.numeric(bird_buff.2016$X)
bird_buff.2016$Y <- as.numeric(bird_buff.2016$Y)

bird.raster.2016 <- bird_buff.2016 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>%  rasterize(r)

rpoints.2016 <- rasterToPoints(bird.raster.2016) %>% data.frame()

rm(bird_buff.2016, bird.raster.2016)
#################################################################
#2017
bird_buff.2017 <- ebird_buff$data[8] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2017$year <- as.integer(bird_buff.2017$year)
bird_buff.2017$X <- as.numeric(bird_buff.2017$X)
bird_buff.2017$Y <- as.numeric(bird_buff.2017$Y)

bird.raster.2017 <- bird_buff.2017 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% rasterize(r)

rpoints.2017 <- rasterToPoints(bird.raster.2017) %>% data.frame()

rm(bird_buff.2017, bird.raster.2017)
################################################################
#2018
bird_buff.2018 <- ebird_buff$data[9] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2018$year <- as.integer(bird_buff.2018$year)
bird_buff.2018$X <- as.numeric(bird_buff.2018$X)
bird_buff.2018$Y <- as.numeric(bird_buff.2018$Y)

bird.raster.2018 <- bird_buff.2018 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% rasterize(r)

rpoints.2018 <- rasterToPoints(bird.raster.2018) %>% data.frame()

rm(bird_buff.2018, bird.raster.2018)
###############################################################
#2019
bird_buff.2019 <- ebird_buff$data[10] %>%
  bind_cols %>%
  st_cast(to = "POINT") %>%
  dplyr::mutate(
    X =  sf::st_coordinates(geometry)[,1], #retrieve X coord
    Y =  sf::st_coordinates(geometry)[,2]  #retrieve Y coord
  ) %>%
  sf::st_drop_geometry()

bird_buff.2019$year <- as.integer(bird_buff.2019$year)
bird_buff.2019$X <- as.numeric(bird_buff.2019$X)
bird_buff.2019$Y <- as.numeric(bird_buff.2019$Y)

bird.raster.2019 <- bird_buff.2019 %>% st_as_sf(coords = c("X", "Y"), crs = 4326) %>% rasterize(r)

rpoints.2019 <- rasterToPoints(bird.raster.2019) %>% data.frame()

rm(bird_buff.2019, bird.raster.2019)




ebird_buff <- ebird_buff$data %>% map(~ mutate(., ID = row_number()))
ebird_buff <- ebird_buff %>% bind_rows(.id = 'grp')

ebird_buff$grp <- NULL
ebird_buff$ID <- NULL

ebird_buff <- ebird_buff %>% mutate(id = row_number())


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
pland <- pland %>% distinct()%>%
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0))

write_csv(pland, "data/pland_location-year.csv")

ebird_buff <- ebird_buff %>% arrange(year)

#year selection
ebird_buff.2010 <- ebird_buff[ebird_buff$year %in% 2010,]
ebird_buff.2011 <- ebird_buff[ebird_buff$year %in% 2011,]
ebird_buff.2012 <- ebird_buff[ebird_buff$year %in% 2012,]
ebird_buff.2013 <- ebird_buff[ebird_buff$year %in% 2013,]
ebird_buff.2014 <- ebird_buff[ebird_buff$year %in% 2014,]
ebird_buff.2015 <- ebird_buff[ebird_buff$year %in% 2015,]
ebird_buff.2016 <- ebird_buff[ebird_buff$year %in% 2016,]
ebird_buff.2017 <- ebird_buff[ebird_buff$year %in% 2017,]
ebird_buff.2018 <- ebird_buff[ebird_buff$year %in% 2018,]
ebird_buff.2019 <- ebird_buff[ebird_buff$year %in% 2019,]




agg_factor <- round(2 * neighborhood_radius / res(landcover))
r.2010 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2010.1 <- ebird_buff.2010 %>% 
  st_transform(crs = projection(r.2010)) %>% 
  fasterize(r.2010) %>% 
  # remove any empty cells at edges
  trim()

rr2010 <- r.2010

r.2011 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2011.1 <- ebird_buff.2011 %>% 
  st_transform(crs = projection(r.2011)) %>% 
  fasterize(r.2011,) %>% 
  # remove any empty cells at edges
  trim()

r.2012 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2012 <- ebird_buff.2012 %>% 
  st_transform(crs = projection(r.2012)) %>% 
  fasterize(r.2012, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2013 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2013 <- ebird_buff.2013 %>% 
  st_transform(crs = projection(r.2013)) %>% 
  fasterize(r.2013, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2014 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2014 <- ebird_buff.2014 %>% 
  st_transform(crs = projection(r.2014)) %>% 
  fasterize(r.2014, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2015 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2015 <- ebird_buff.2015 %>% 
  st_transform(crs = projection(r.2015)) %>% 
  fasterize(r.2015, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2016 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2016 <- ebird_buff.2016 %>% 
  st_transform(crs = projection(r.2016)) %>% 
  fasterize(r.2016, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2017 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2017 <- ebird_buff.2017 %>% 
  st_transform(crs = projection(r.2017)) %>% 
  fasterize(r.2017, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2018 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2018 <- ebird_buff.2018 %>% 
  st_transform(crs = projection(r.2018)) %>% 
  fasterize(r.2018, field = c("species_observed")) %>% 
  # remove any empty cells at edges
  trim()

r.2019 <- raster(landcover) %>% 
  aggregate(agg_factor) 
r.2019.1 <- ebird_buff.2019 %>% 
  st_transform(crs = projection(r.2019)) %>% 
  fasterize(r.2019, field = c("species_observed", "year")) %>% 
  # remove any empty cells at edges
  trim()


agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 


rat <- levels(r)[[1]]
rat[["landcover"]] <- c("oak","pine", "scrub", "grass")
levels(r_fact) <- rat


r.all <- ebird_buff %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r) %>% 
  # remove any empty cells at edges
  trim()




r.ratify <- deratify(r.all, c("species_observed", "year"))

rpointsall <-  rasterToPoints(r.ratify,  spatial="TRUE") %>% as.data.frame()
rpointsall <- rpointsall[,c(2, 3, 1)]
r.all <- rasterFromXYZ(rpointsall)

#covert to raserpoints then dataframe
rpoints2010 <-  rasterToPoints(r.2010,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2010)
rpoints2011 <-  rasterToPoints(r.2011,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2011)
rpoints2012 <-  rasterToPoints(r.2012,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2012)
rpoints2013 <-  rasterToPoints(r.2013,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2013)
rpoints2014 <-  rasterToPoints(r.2014,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2014)
rpoints2015 <-  rasterToPoints(r.2015,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2015)
rpoints2016 <-  rasterToPoints(r.2016,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2016)
rpoints2017 <-  rasterToPoints(r.2017,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2017)
rpoints2018 <-  rasterToPoints(r.2018,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2018)
rpoints2019 <-  rasterToPoints(r.2019,  spatial="TRUE") %>% as.data.frame() %>% mutate(year = 2019)

rpoints.all <- rbind(rpoints2010,rpoints2011,rpoints2012,rpoints2013,rpoints2014,rpoints2015,rpoints2016,rpoints2017,rpoints2018,rpoints2019)

#filter
rpoints2010 <- rpoints2010[, c(2, 3, 1)]
rpoints2011 <- rpoints2011[, c(2, 3, 1)]
rpoints2012 <- rpoints2012[, c(2, 3, 1)]
rpoints2013 <- rpoints2013[, c(2, 3, 1)]
rpoints2014 <- rpoints2014[, c(2, 3, 1)]
rpoints2015 <- rpoints2015[, c(2, 3, 1)]
rpoints2016 <- rpoints2016[, c(2, 3, 1)]
rpoints2017 <- rpoints2017[, c(2, 3, 1)]
rpoints2018 <- rpoints2018[, c(2, 3, 1)]
rpoints2019 <- rpoints2019[, c(2, 3, 1)]


r.1<-rasterFromXYZ(rpoints2010[,c(1:3)])
r.2<-rasterFromXYZ(rpoints2011[,c(1:3)])
r.3<-rasterFromXYZ(rpoints2012[,c(1:3)])
r.4<-rasterFromXYZ(rpoints2013[,c(1:3)])
r.5<-rasterFromXYZ(rpoints2014[,c(1:3)])
r.6<-rasterFromXYZ(rpoints2015[,c(1:3)])
r.7<-rasterFromXYZ(rpoints2016[,c(1:3)])
r.8<-rasterFromXYZ(rpoints2017[,c(1:3)])
r.9<-rasterFromXYZ(rpoints2018[,c(1:3)])
r.10<-rasterFromXYZ(rpoints2019[,c(1:3)])

#get the dataframe into lists
pinch <- do.call(rbind, st_geometry(ebird_buff.2010$geometry)) %>% 
  as_tibble() %>% setNames(c("X","Y"))





















r <- writeRaster(r, filename = "data/prediction-surface.tif", overwrite = TRUE)

r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% transmute(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)
# extract landcover values within neighborhoods, only needed most recent year
lc_extract_pred <- landcover %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019)) %>% 
  tibble(id = r_cells$id, data = .) %>% 
  unnest(data)
# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>% 
  count(id, y2010, y2011, y2012, y2013, y2014, y2015, y2016, y2017, y2018, y2019) %>% 
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(.))
#take means
#pivot longer for landcover names
pland_pred <- pland_pred %>% pivot_longer(-c(1,12))
#change column names and remove 'y' from year
colnames(pland_pred) <- c("id","pland", "year", "landcover")
pland_pred$year <- sub('.', '', pland_pred$year)
# convert names to be more descriptive
pland_pred <- pland_pred %>%
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)
# tranform to wide format, filling in implicit missing values with 0s
pland_pred <- pland_pred %>% distinct()%>%
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% 
  select(id, year, everything())
# join in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred, by = "id")


forest_cover <- pland_coords %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

writeRaster(forest_cover, "output/single_layer_forestcover.tif")