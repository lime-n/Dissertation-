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
#This is to calculate the percentage of pland cover for each BCR region.

bcr <- read_sf("data/gis-data.gpkg", "bcr") %>% filter(bcr_code %in% c(11, 12, 13, 17, 18, 19,21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)) %>% st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                                                                                                                                                                        "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))

landcover <- list.files("data/modis", "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% 
  stack()
# label layers with year
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- bcr %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r) 
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% mutate(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

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

#read lc_extract_pred and extract individual years
lc_extract_pred <- read.csv("data/layer_lc_extract_pred.csv")
lc_extract_pred <- lc_extract_pred[, -1]
lc_extract_pred <- lc_extract_pred[complete.cases(lc_extract_pred),]

########################################
### Linear models

lay.all <- lc_extract_pred %>% pivot_longer(-c(1, 2, 13))

#test_begin
test.t <- lay.all %>% 
  count(layer, value, name) %>% 
  group_by(layer) %>%
  mutate(pland = (n / sum(n)*1000)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(.))

colnames(test.t) <- c("layer", "landcover","year", "pland")
test.t$year <- sub('.', '', test.t$year)

test.t.land <- test.t %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)


test.t.land <- test.t.land %>% 
  inner_join(ly_names, by = "layer") %>% 
  arrange(layer) %>% 
  select(-layer)

test.t.land <- test.t.land %>% pivot_wider(names_from = ly_name, values_from = pland, values_fill = list(pland = 0))


ly_names <- tibble(layer = 1:17, 
                   ly_name = c("Prairie Potholes"     ,
                               "Boreal Hardwood Transition"   ,  
                               "Lower Great Lakes/St. Lawrence Plain",
                               "Badlands And Prairies"    ,    
                               "Shortgrass Prairie"    ,      
                               "Central Mixed Grass Prairie"   ,      
                               "Oaks And Prairies"      ,   
                               "Eastern Tallgrass Prairie"  ,
                               "Prairie Hardwood Transition" ,       
                               "Central Hardwoods"     ,  
                               "West Gulf Coastal Plain/Ouachitas",
                               "Mississippi Alluvial Valley" ,        
                               "Southeastern Coastal Plain" ,   
                               "Appalachian Mountains" ,        
                               "Piedmont"      ,                      
                               "New England/Mid-Atlantic Coast" , 
                               "Peninsular Florida" ))

 

#####################################
#Linear models end


##############################################################
#Test begin; Select those values relative to habitat
water.2010.0<- test_pland.pivot[test_pland.pivot$lc_name %in% "pland_00_water", ]
colnames(water.2010.0) <- c("id", "layer", "pland_00_water")

#2019
water.2019.0<- test_pland.pivot[test_pland.pivot$lc_name %in% "pland_00_water", -4]
colnames(water.2019.0) <- c("id", "layer", "pland_00_water")

#############################################################
#Join coordinates
#Select the layer
r.water <- r_centers[r_centers$layer %in% 1, ]
#

water.t.2019.0 <- st_transform(r.water, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r.water$id, .) %>%
  rename(longitude = X, latitude = Y) %>% 
  full_join(water.2019.0, by = c("id"))

#Create zeros for area with no values
water.t.2019.0$pland_00_water[is.na(water.t.2019.0$pland_00_water)] <- 0
water.t.2019.0$layer[is.na(water.t.2019.0$layer)] <- 1


#############################################################


forest.2019 <- water.t.2019.0%>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs("ESRI:102003")$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()
plot(forest.2019)
