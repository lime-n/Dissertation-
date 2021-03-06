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
k <- read_csv("data/layer.single/single.pland-elev_location-year.csv")
z <- read_csv("data/layer.single/single.pland-slope_location-year.csv")


habitat <- inner_join(k, z, by = c("id", "locality_id", "year", "pland_00_water", "pland_01_evergreen_needleleaf", "pland_02_evergreen_broadleaf",
                                  "pland_03_deciduous_needleleaf", "pland_04_deciduous_broadleaf", "pland_05_mixed_forest", 
                                  "pland_06_closed_shrubland", "pland_07_open_shrubland", "pland_08_woody_savanna", 
                                  "pland_09_savanna", "pland_10_grassland", "pland_11_wetland", "pland_12_cropland",
                                  "pland_13_urban", "pland_14_mosiac", "pland_15_barren"))


habitat$year <- as.integer(habitat$year)
rm(k, z)
# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))
# prediction surface
j <- read_csv("data/layer.single/single.pland2010-elev_prediction-surface.csv")
l  <-  read_csv("data/layer.single/single.pland2010-slope_prediction-surface.csv")
pred_surface_2010 <- inner_join(j, l, by = c("id","layer.x", "layer.y", "longitude", "latitude", "year", "pland_00_water", "pland_01_evergreen_needleleaf", "pland_02_evergreen_broadleaf",
                                        "pland_03_deciduous_needleleaf", "pland_04_deciduous_broadleaf", "pland_05_mixed_forest", 
                                        "pland_06_closed_shrubland", "pland_07_open_shrubland", "pland_08_woody_savanna", 
                                        "pland_09_savanna", "pland_10_grassland", "pland_11_wetland", "pland_12_cropland",
                                        "pland_13_urban", "pland_14_mosiac", "pland_15_barren"))



rm(j, l)
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
  select(species_observed, year, 
         time_observations_started, duration_minutes,
         effort_distance_km, number_observers,
         starts_with("pland"),
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
# calculate the average observed encounter rates for different 
# categories of estimated encounter rates 
average_encounter <- rf_pred_train %>%
  mutate(pred_cat = cut(rf_pred_train$pred, breaks = seq(0, 1, by=0.02))) %>%
  group_by(pred_cat) %>%
  summarise(pred = mean(pred), obs = mean(obs), checklist_count = n()) %>%
  ungroup()
# plot
cal_pred <- tibble(pred = seq(0, 1, length.out = 100))
cal_pred <- predict(calibration_model, cal_pred, type = "response") %>% 
  bind_cols(cal_pred, calibrated = .)
ggplot(cal_pred) +
  aes(x = pred, y = calibrated) +
  geom_line() +
  geom_point(data = average_encounter, 
             aes(x = pred, y = obs, size = sqrt(checklist_count)),
             show.legend = FALSE, shape = 1) +
  labs(x = "Estimated encounter rate",
       y = "Observed encounter rate",
       title = "Calibration model")

p_fitted <- predict(rf, data = ebird_split$test, type = "response")
# extract probability of detection
p_fitted <- p_fitted$predictions[, 2]
# calibrate
p_calibrated <- predict(calibration_model, 
                        newdata = tibble(pred = p_fitted), 
                        type = "response")
rf_pred_test <- data.frame(id = seq_along(p_calibrated),
                           # actual detection/non-detection
                           obs = ebird_split$test$species_observed,
                           # uncalibrated prediction
                           fit = p_fitted,
                           # calibrated prediction
                           cal = p_calibrated) %>%
  # constrain probabilities to 0-1
  mutate(cal = pmin(pmax(cal, 0), 1)) %>% 
  drop_na()
# mean squared error (mse)
mse_fit <- mean((rf_pred_test$obs - rf_pred_test$fit)^2, na.rm = TRUE)
mse_cal <- mean((rf_pred_test$obs - rf_pred_test$cal)^2, na.rm = TRUE)
# pick threshold to maximize kappa
opt_thresh <- optimal.thresholds(rf_pred_test, opt.methods = "MaxKappa")
# calculate accuracy metrics: auc, kappa, sensitivity, specificity,
metrics_fit <- rf_pred_test %>% 
  select(id, obs, fit) %>% 
  presence.absence.accuracy(threshold = opt_thresh$fit, 
                            na.rm = TRUE, 
                            st.dev = FALSE)
metrics_cal <- rf_pred_test %>% 
  select(id, obs, cal) %>% 
  presence.absence.accuracy(threshold = opt_thresh$cal, 
                            na.rm = TRUE, 
                            st.dev = FALSE)
rf_assessment <- tibble(
  model = c("RF", "Calibrated RF"),
  mse = c(mse_fit, mse_cal),
  sensitivity = c(metrics_fit$sensitivity, metrics_cal$sensitivity),
  specificity = c(metrics_fit$specificity, metrics_cal$specificity),
  auc = c(metrics_fit$AUC, metrics_cal$AUC),
  kappa = c(metrics_fit$Kappa, metrics_cal$Kappa)
)
knitr::kable(rf_assessment, digits = 3)

pi <- enframe(rf$variable.importance, "predictor", "importance")
# plot
ggplot(pi) + 
  aes(x = fct_reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, size = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL, 
       y = "Predictor Importance (Gini Index)") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", size = 0.5))

top_pred <- pi %>% 
  filter(!predictor %in% c("year", "day_of_year")) %>% 
  top_n(n = 15, wt = importance) %>% 
  arrange(desc(importance))

calculate_pd <- function(predictor, model, data, 
                         x_res = 25, n = 1000) {
  # create prediction grid
  rng <- range(data[[predictor]], na.rm = TRUE)
  x_grid <- seq(rng[1], rng[2], length.out = x_res)
  grid <- data.frame(covariate = predictor, x = x_grid, 
                     stringsAsFactors = FALSE)
  names(grid) <- c("covariate", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  s <- sample(seq.int(nrow(data)), size = n, replace = FALSE)
  data <- data[s, ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict
  p <- predict(model, data = grid)
  
  # summarize
  pd <- grid[, c("covariate", predictor)]
  names(pd) <- c("covariate", "x")
  pd$pred <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, covariate, x) %>% 
    dplyr::summarise(pred = mean(pred, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pd)
}

pd <- top_pred %>% 
  mutate(pd = map(predictor, calculate_pd, model = rf, 
                  data = ebird_split$train),
         pd = map(pd, ~ .[, c(2, 3)]),
         pd = map(pd, set_names, nm = c("value",  "encounter_rate"))) %>% 
  unnest(cols = pd)
# calibrate predictions
pd$encounter_rate <- predict(calibration_model, 
                             newdata = tibble(pred = pd$encounter_rate), 
                             type = "response") %>% 
  as.numeric()

# plot
ggplot(pd) +
  aes(x = value, y = encounter_rate) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ as_factor(predictor), nrow = 3, scales = "free") +
  labs(x = NULL, y = "Encounter Rate") +
  theme_minimal() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"))

pd_time <- calculate_pd("time_observations_started",
                        model = rf, 
                        data = ebird_split$train,
                        # make estimates at 30 minute intervals
                        # using a subset of the training dataset
                        x_res = 2 * 24, n = 1000) %>% 
  transmute(time_observations_started = x, encounter_rate = pred)
# histogram
g_hist <- ggplot(ebird_split$train) +
  aes(x = time_observations_started) +
  geom_histogram(binwidth = 1, center = 0.5, color = "grey30",
                 fill = "grey50") +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Hours since midnight",
       y = "# checklists",
       title = "Distribution of observation start times")
# gam
g_pd <- ggplot(pd_time) +
  aes(x = time_observations_started, y = encounter_rate) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(x = "Hours since midnight",
       y = "Probability of reporting",
       title = "Observation start time partial dependence")
# combine
grid.arrange(g_hist, g_pd)

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
pred_er_2010 <- bind_cols(pred_surface_eff_2010, encounter_rate = pred_rf_cal_2010) %>% 
  select(latitude, longitude, encounter_rate, id, layer.y) %>% 
  mutate(encounter_rate = pmin(pmax(encounter_rate, 0), 1))

############
test.2010 <- pred_er_2010 %>% 
  count(layer.y, encounter_rate) %>% 
  group_by(layer.y) 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(.))

x.test <- function(x){

  rp <- sum(pred_er_2010[pred_er_2010$layer.y %in% x, c(3, 5)]$encounter_rate)
  return(rp)
  }
############



r_pred_2010 <- pred_er_2010 %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize
  rasterize(r)
r_pred_2010 <- r_pred_2010[[-1]]
# save the raster
tif_dir <- "output"
if (!dir.exists(tif_dir)) {
  dir.create(tif_dir)
}
writeRaster(r_pred, file.path(tif_dir, "rf-rehwoo2010.tif"), 
            overwrite = TRUE)

r_pred_proj <- projectRaster(r_pred, crs = map_proj$proj4string, method = "ngb")
par(mar = c(3.5, 0.25, 0.25, 0.25))
# set up plot area
plot(bcr, col = NA, border = NA)
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
# encounter rate
r_max <- ceiling(10 * cellStats(r_pred_proj, max)) / 10
brks <- seq(0, r_max, by = 0.025)
lbl_brks <- seq(0, r_max, by = 0.1)
# ebird status and trends color palette
pal <- abundance_palette(length(brks) - 1)
plot(r_pred_proj, 
     col = pal, breaks = brks, 
     maxpixels = ncell(r_pred_proj),
     legend = FALSE, add = TRUE)
# borders
plot(bcr, border = "#000000", col = NA, lwd = 1, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
box()
# legend
par(new = TRUE, mar = c(0, 0, 0, 0))
title <- "Red-headed Woodpecker Encounter Rate"
image.plot(zlim = range(brks), legend.only = TRUE, 
           col = pal, breaks = brks,
           smallplot = c(0.25, 0.75, 0.06, 0.09),
           horizontal = TRUE,
           axis.args = list(at = lbl_brks, labels = lbl_brks,
                            fg = "black", col.axis = "black",
                            cex.axis = 0.75, lwd.ticks = 0.5,
                            padj = -1.5),
           legend.args = list(text = title,
                              side = 3, col = "black",
                              cex = 1, line = 0))

######################################################################################################################################
