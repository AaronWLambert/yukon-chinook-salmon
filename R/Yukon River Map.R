library (tidyverse)
library (sf)
library (rnaturalearth)
library (rnaturalearthdata)

library (ggmap)
library (mapdata)
library (maps)
library (viridis)
library (stringr)
library (magick)
library (sp)
library (rgdal)
library (ggspatial) #requires that you have installed "sf"
library (patchwork)
library(sf)
library (terra) #for dealing with spat raster objects
library (tidyterra) #
# Define Workflow Paths ============================================

# set to working directory
wd <- "C:/Users/aaron/Desktop/Yukon Kings/Inseason Forcast Model"

setwd(wd)

# Objects used to save/load data, outputs, or stan/R scripts
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"figs")
dir.stan <- file.path(wd,"Stan")
dir.data <- file.path(wd,"Data")
dir.R <- file.path(wd,"R")


## 1) Adding layers with shapefiles
#First, let's start with our base Alaska map using what we learned in lecture 12. We're going to 
#pretend that Alaska stops at the 180th meridian.

world <- ne_countries (scale = "medium", returnclass = "sf")
# world <-st_transform(world, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #Transform into WGS84 coordinate system
# world <- st_make_valid(st_as_sf(world)) 
ggplot (data = world) + geom_sf() + coord_sf(xlim = c(-170, -160), ylim = c(50, 75)) + theme_bw()


# Downloaded from natural earth. Read in from shp file
ne50rivers <- st_read(file.path(dir.data,"Map Data/Rivers Data/ne_50m_rivers_lake_centerlines.shp"))

ne10rivers <- st_read(file.path(dir.data,"Map Data/Rivers Data/ne_10m_rivers_lake_centerlines.shp"))

ne10rivers_NA <- st_read(file.path(dir.data,"Map Data/Rivers Data/ne_10m_rivers_north_america.shp"))

test.yukon <- ne50rivers[ne50rivers$label=="Yukon",]

test.koyukuk <- ne50rivers[ne50rivers$label=="Koyukuk" ,]

test.porcupine <- ne50rivers[ne50rivers$label=="Porcupine",]

test.tanana <- ne50rivers[ne50rivers$label=="Tanana",]

test.teslin <- ne50rivers[ne50rivers$label=="Teslin",]

yukonRiver <- rbind(test.koyukuk,test.porcupine,test.tanana,test.teslin,test.yukon)

plot.sonar.DF <- data.frame("Site" = c("Eagle Sonar", "Pilot Station Sonar"),
                            "lat" = c(64.7860222,61.93605),
                            "lon" = c(-141.199917,-162.883403))
yukonRiver$Site <- "Yukon River"

# png(filename = file.path(dir.figs,"Yukon Map Chap 1 21Sept23.png"),width = 800, height = 800)
# Lets see....
ggplot () + 
  geom_sf(data = world) + 
  geom_sf(data = ne50rivers) + 
  geom_sf(data = yukonRiver, aes(shape = Site, color = Site))+
  # geom_sf(data = ne10rivers) + 
  # geom_sf(data = ne10rivers_NA) + 
  # geom_sf(data = test.yukon, color = "blue")+
  # geom_sf(data = test.koyukuk, color = "blue")+
  # geom_sf(data = test.porcupine, color = "blue")+
  # geom_sf(data = test.join, color = "red")+
  # geom_sf(data = test.tanana, color = "blue")+
  # geom_sf(data = test.teslin, color = "blue")+
  coord_sf(xlim = c(-135, -175), 
           ylim = c(50, 75)) + 
  geom_point(data = plot.sonar.DF, aes( x= lon, y = lat, shape = Site, color = Site), size = 3)+
  scale_shape_manual(values = c(15,17,3), name = "")+
  scale_color_manual(values = c("black","black","blue"), name = "")+
  scale_linetype_manual(values = c(1,1,1), name = "")+
  theme_bw()+
  labs(x = "Longitude", y = "Latitude")+
  theme(
    legend.position = c(.80,.933),
    # plot.background = element_rect(fill = "transparent", colour = NA),
    legend.title=element_blank(),
    legend.text=element_text(size=14),
    legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
    legend.key = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )

# dev.off()
# Try to make prettier?
alaska.box <- c(left = -180, bottom = 50, right = -125, top = 73)

alaska.terrain <- ggmap::get_stamenmap(bbox = alaska.box, zoom = 5, maptype = "terrain-background",where = "cache")

# Plot the selected map
ggmap(alaska.terrain)

# Function to allow plotting with crs 3875
ggmap_bbox <- function(map) {
  if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(map, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
  
  # Convert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
  map
}

alaska.terrain.fixed <- ggmap_bbox(alaska.terrain)

# Try with rivers
ggmap(alaska.terrain.fixed)+
  coord_sf(crs = st_crs(3857))+
  geom_sf(data = st_transform( test.yukon, crs = 3857), inherit.aes = FALSE, color = "blue")+
  geom_sf(data = st_transform( test.koyukuk, crs = 3857), inherit.aes = FALSE, color = "blue")+
  geom_sf(data = st_transform( test.porcupine, crs = 3857), inherit.aes = FALSE, color = "blue")+
  geom_sf(data = st_transform( test.tanana, crs = 3857), inherit.aes = FALSE, color = "blue")+
  geom_sf(data = st_transform( test.teslin, crs = 3857), inherit.aes = FALSE, color = "blue")+
  # coord_sf(xlim = c(-135, -175), 
  #          ylim = c(50, 75)) + 
  geom_point(data = plot.sonar.DF, aes( x= lon, y = lat, shape = Site), size = 2)
  