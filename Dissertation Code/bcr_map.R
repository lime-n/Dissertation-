library(broom)
library(rgdal)

p <- readOGR(dsn = "bcr.gpkg")
tidy_bcr <- tidy(p)

ggplot(p, aes(x = long, y = lat, group = group)) +
  geom_polygon(color = "black", size = 0.1, fill = "lightgrey") +
  coord_equal() +
  theme_minimal()

p$id <- row.names(p)
tidy_bcr <- left_join(tidy_bcr, p@data)
bcrcop <- data.frame(bcr_name = sort(p@data$bcr_name),
                     bcr_code = c(28,17,12,24,19,22,13,26,30,21,31,29,23,11,18,27,25))
tidy_bcr <- left_join(tidy_bcr, bcrcop)

bcrLabel <- tidy_bcr %>%
  group_by(bcr_name) %>%
  summarise(label_long = mean(range(long)), label_lat = mean(range(lat)), bcr_code = mean(bcr_code)) 
  mutate(label_long = replace(label_long, bcr_name %in% c("Badlands And Prairies", "Shortgrass Prairie ", " Central Mixed Grass Prairie ", "Oaks And Prairies ", "West Gulf Coastal Plain/Ouachitas ", "Peninsular Florida  ", "Central Hardwoods","Eastern Tallgrass Prairie","Prairie Hardwood Transition","Boreal Hardwood Transition","Lower Great Lakes/St. Lawrence Plain","New England/Mid-Atlantic Coast"), 
                              values=c(-105.1,-103.1,-100.4,-97.1,-94.1,-82,-90,-89.7,-89.4,-83.8,-77,-73.7)),
         label_lat = replace(label_lat, bcr_name %in% c("Prairie Potholes", "Southeastern Coastal Plain ", "Piedmont  ", "Appalachian Mountains ","Mississippi Alluvial Valley "), 
                             values=c(48.5,32.5 ,36.3 ,38.3,33.9 )))

map <- ggplot(tidy_bcr, aes(x = long, y = lat, group = group, fill = BCR_regions)) +
  geom_polygon(color = "black", size = 0.1) +
  coord_equal() +
  theme_void()  +
  theme(plot.title = element_text(margin = margin(t = 40, b = -40)))

c <- map + geom_text(data = bcrLabel, mapping = aes(x=label_long, y = label_lat, label = bcr_code, group = NA), cex = 4, col = "white") + theme(legend.box.background = element_rect(colour = "black"), legend.position = "bottom")
c + 
