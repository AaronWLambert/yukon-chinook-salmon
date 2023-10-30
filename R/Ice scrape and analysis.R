# Aaron Lambert
# 2/21/2023
# Yukon Inseason Chinook Forecast
# Sea ICe
# https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.htmlnegative

# NOAA OISST on ERDDAP server
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

all.timing.logistic <- read.csv(file = file.path(dir.output,"logistic curve parameters All Chinook 1995_2022.csv"))
# Following scrape code from Robert W Schlegel and AJ Smit

# The packages we will use
library(dplyr) # A staple for modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(ggpubr)
library(ggmap)
library(maps)
library(usmap)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(ggsn)
library(cowplot)

# The information for the NOAA OISST data
# rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")

# # Note that there is also a version with lon values from 0 yo 360
# rerddap::info(datasetid = "ncdcOisst21Agg", url = "https://coastwatch.pfeg.noaa.gov/erddap/")
# 
# # This function downloads and prepares data based on user provided start and end dates
# OISST_sub_dl <- function(time_df){
#   OISST_dat <- griddap(datasetx = "ncdcOisst21Agg_LonPM180",
#                        url = "https://coastwatch.pfeg.noaa.gov/erddap/",
#                        time = c(time_df$start, time_df$end),
#                        zlev = c(0, 0),
#                        latitude = c(60, 64),
#                        longitude = c(-163,-169),
#                        fields = "ice")$data %>%
#     mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>%
#     dplyr::rename(t = time) %>%
#     select(longitude, latitude, t, ice) %>%
#     na.omit()
# }
# 
# # Date download range by start and end dates per year
# dl_years <- data.frame(date_index = 1:4,
#                        # start = as.Date(c("1982-01-01", "1990-01-01",
#                        #                   "1998-01-01", "2006-01-01", "2014-01-01")),
#                        # end = as.Date(c("1989-12-31", "1997-12-31",
#                        #                 "2005-12-31", "2013-12-31", "2019-12-31")))
#                        start = as.Date(c("1995-01-01", "2003-01-01",
#                                          "2011-01-01", "2020-01-01")),
#                        end = as.Date(c("2002-12-31", "2010-12-31",
#                                        "2019-12-31", "2022-12-31")))
# 
# # Download all of the data with one nested request
# # The time this takes will vary greatly based on connection speed
# system.time(
#   OISST_data <- dl_years %>%
#     group_by(date_index) %>%
#     group_modify(~OISST_sub_dl(.x)) %>%
#     ungroup() %>%
#     select(longitude, latitude, t, ice)
# )

# saveRDS(OISST_data, file = file.path(dir.data,"OISST_PIC_data_163_169_60_64_25Oct23.RDS"))
OISST_data <- readRDS(file = file.path(dir.data,"OISST_PIC_data_163_169_60_64_25Oct23.RDS"))

# Break out day, month , year , and jdate
yukon.ice  <- OISST_data %>% mutate(year = lubridate::year(t),
                                    month = lubridate::month(t),
                                    day = lubridate::day(t),
                                    jday = lubridate::yday(t))

# saveRDS(object = yukon.ice, file = file.path(dir.data,"PIC Data full 29Aug.RDS"))

# Mapping stuff....
# Get Map info for plotting below
ak <- st_read(dsn = file.path(dir.data,"Map Data/Alaska_Coastline")) #load shape file from right directory
ak <- st_transform(ak, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #Transform into WGS84 coordinate system
ak <- st_make_valid(st_as_sf(ak)) # Make sf object valid


# Canada
can <- st_read(dsn = file.path(dir.data, "Map Data/Canadian_Provinces/PROVINCE.shp"))
can <- st_transform(can, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") #Transform into WGS84 coordinate system
can <- st_make_valid(st_as_sf(can)) # Make sf object valid

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

# Create Inset Map 
inset.map <- 
  ggplot()+
  geom_sf(data = ak) +
  geom_sf(data = can)+# call in sf object
  geom_rect(mapping = aes(xmin = -163, xmax = -169, ymin = 62, ymax = 64), 
            color = "red", size = 1) + # geom_rect specifies the boundaries of your rectangle
  coord_sf(xlim = c(-180, -137), ylim = c(50, 71))+ # controlling axis to disregard dateline
  theme_classic() +
  theme(axis.title = element_blank(),  # Removing remnants of axis title, text, tick marks,
        axis.text = element_blank(),  # lines, and changing plot background to transparent,
        axis.ticks = element_blank(), # and also creating black border around AK map
        axis.line = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75))


# Subset for May observations
yukon.ice.may <- yukon.ice[yukon.ice$month==5 & yukon.ice$year != 1996,]

# Caluclate the mean for May
yukon.ice.may.mean <- yukon.ice.may %>% 
  group_by(longitude, latitude, year) %>% 
  summarize(mean.ice = mean(ice))

# Get the unique long and lat for looping through correlations
lon <- unique(yukon.ice.may$longitude)
lat <- unique(yukon.ice.may$latitude)


# Plot out correlations for all years of PSS available ##############################
# Vector of midpoints
midpoint <- all.timing.logistic[,c(2,5)]

# Matrix to hold the values from correlation loop
cor.mat <- matrix(nrow = length(lat)*length(lon), ncol = 3)

# Name the matrix
colnames(cor.mat) <- c("lon", "lat", "cor")

# Set counter
counter <- 1
# i <- 15
# a <- 15
for (i in 1:length(lat)) {
  for (a in 1:length(lon)) {
    
    ice <- yukon.ice.may.mean[yukon.ice.may.mean$longitude == lon[a] & yukon.ice.may.mean$latitude==lat[i],]
    
    temp.df <- left_join(ice,midpoint)
    
    if(!is.na(ice[1,1]) ){
      
      
      cor.ice <- cor(temp.df$mean.ice,temp.df$mid)
      
      
      
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- cor.ice
    }else{
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- NA
    }
    
    counter <- counter +1
    
  }
  
}

# Convert matrix to df
cor.df <- data.frame(cor.mat)

# Look for the max correlation
cor.df[which.max(cor.df$cor),]$lat


df<- left_join(yukon.ice.may.mean[yukon.ice.may.mean$longitude == cor.df[which.max(cor.df$cor),]$lon &
                               yukon.ice.may.mean$latitude== cor.df[which.max(cor.df$cor),]$lat &
                               yukon.ice.may.mean$year >= 1995,], midpoint)

ggplot(df, aes(x = mean.ice, y = mid))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~`,`~")),label.x = 0)+
  labs(x = "SST (Lat = 63.375, Lon = -168.375)",
       y = "PSS Chinook Salmon Run Midpoint")


pic_map_1995_2022 <- ggplot(data = ak) +
  geom_raster(data = cor.df, aes(x = lon, y = lat, fill = cor),interpolate = T) +
  #   # geom_sf(data = temp.dat)+
  # scale_fill_gradientn(colours = rev(rainbow(7)), na.value = NA) +
  scale_fill_viridis_c(direction = -1, na.value = NA,name = "Correlation (r)")+
  # ggtitle("Correlation Heat Map", subtitle = "May Mean PIC vs Midpoint \n(1995-2022)")+# use ak map data
  geom_sf(fill = "darkgray") +
  geom_sf(data = ne50rivers)+# plot sf object
  geom_sf(data = ne10rivers_NA)+
  geom_sf(data = ne10rivers)+
  coord_sf(xlim = c(-170, -162), ylim = c(60, 64)) + # specify region-wide map coordinates (dimensions)
  # guides(colour = guide_colourbar(barwidth  = unit(30,"cm"),
  #                                 title.position = "top", 
  #                                 title.hjust = 0.5))+
  labs(x = "Longitude", y = "Latitude") +
  # theme_bw()+
  theme(
    # axis.text = element_text(size = 16, color = "black"),
    text = element_text(size = 15),
    axis.line = element_blank(),
    panel.grid = element_blank(), # remove elements for clean background
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.title=element_text(size=15),
    legend.text=element_text(size=12),
    # legend.position = c(0.3, 0.3),
    # legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
    legend.key = element_rect(fill = "white"),
    legend.position = "top", 
    # plot.margin = unit(c(0,0,0,0), "cm")
  )+
  # guides(colour = guide_colourbar(barwidth  = unit(30,"cm"),
  #                                 title.position = "top", 
  #                                 title.hjust = 0.5))
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Save as png
png(file.path(dir.figs,"ICE Cor Map 1995_2022 Plot 26Oct23.png"))

ggdraw() + 
  draw_plot(pic_map_1995_2022) +
  draw_plot(inset.map, width = 0.2, height = 0.2, x = 0.715, y = 0.07)

# Turn off 
dev.off()


# Plot out correlations for  years 2000 onward ##############################

# Vector of midpoints
midpoint <- all.timing.logistic[all.timing.logistic$year>=2000,][,c(2,5)]

# Matrix to hold the values from correlation loop
cor.mat <- matrix(nrow = length(lat)*length(lon), ncol = 3)

colnames(cor.mat) <- c("lon", "lat", "cor")

counter <- 1
# i <- 1
# a <- 1
for (i in 1:length(lat)) {
  for (a in 1:length(lon)) {
    
    
    
    ice <- yukon.ice.may.mean[yukon.ice.may.mean$longitude == lon[a] & yukon.ice.may.mean$latitude==lat[i],]
    
    temp.df <- left_join(midpoint,ice) %>% na.omit()
    
    if(!is.na(ice[1,1]) ){
      
      
      cor.ice <- cor(temp.df$mean.ice,temp.df$mid)
      
      
      
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- cor.ice
    }else{
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- NA
    }
    
    counter <- counter +1
    
  }
  
}

cor.df.trunc <- data.frame(cor.mat)

# DF with largest correlation for plotting on map
h.cor.df <-cor.df.trunc[which.max(cor.df.trunc$cor),]
h.cor.df2 <- cor.df.trunc[which(cor.df.trunc$cor>=0.6),]

# Look at years available for lat lon 62.125 -165.625
unique(yukon.ice$year[yukon.ice$longitude == h.cor.df$lon & 
            yukon.ice$latitude == h.cor.df$lat &
            yukon.ice$month == 5])

# Years available for -168.625 63.125 
unique(yukon.ice$year[yukon.ice$longitude == -168.625 & 
                        yukon.ice$latitude == 63.125 &
                        yukon.ice$month == 5])

# png(filename = file.path(dir.figs,"ICE Cor Map 2000_2022 Plot 29Aug23.png"), width = 1000, height = 800)

pic_map_2000_2022 <- ggplot(data = ak) +
  geom_raster(data = cor.df.trunc, aes(x = lon, y = lat, fill = cor),interpolate = T) +
  #   # geom_sf(data = temp.dat)+
  # scale_fill_gradientn(colours = rev(rainbow(7)), na.value = NA) +
  scale_fill_viridis_c(direction = 1, na.value = NA,name = "Correlation (r)", limits = c(0, 0.9))+
  # ggtitle("Correlation Heat Map", subtitle = "May Mean PIC vs Midpoint \n(1995-2022)")+# use ak map data
  geom_sf(fill = "darkgray") +
  geom_sf(data = ne50rivers)+# plot sf object
  geom_sf(data = ne10rivers_NA)+
  geom_sf(data = ne10rivers)+
  coord_sf(xlim = c(-170, -162), ylim = c(60, 64)) + # specify region-wide map coordinates (dimensions)
  # guides(colour = guide_colourbar(barwidth  = unit(30,"cm"),
  #                                 title.position = "top", 
  #                                 title.hjust = 0.5))+
  labs(x = "Longitude", y = "Latitude") +
  # theme_bw()+
  theme(
    # axis.text = element_text(size = 16, color = "black"),
    text = element_text(size = 15),
    axis.line = element_blank(),
    panel.grid = element_blank(), # remove elements for clean background
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.title=element_text(size=15),
    legend.text=element_text(size=12),
    legend.key.width  = unit(1,"cm"),
    legend.key = element_rect(fill = "white"),
    legend.position = "top", 
    # plot.margin = unit(c(0,0,0,0), "cm")
  )+
  # guides(colour = guide_colourbar(barwidth  = unit(30,"cm"),
  #                                 title.position = "top", 
  #                                 title.hjust = 0.5))
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Save as png
png(file.path(dir.figs,"ICE Cor Map 2000_2022 Plot 26Oct23.png"))

ggdraw() + 
  draw_plot(pic_map_2000_2022) +
  draw_plot(inset.map, width = 0.2, height = 0.2, x = 0.715, y = 0.07)

# Turn off 
dev.off()

cor.df.trunc[which(cor.df.trunc$cor>=0.6),]




# Data frame of mean sst at the location with the highest correlation when restricting to 2000-2022
# full.df <- data.frame("mid" = midpoint.trunc,
#                       "sst" = yukon.ice.may.mean[yukon.ice.may.mean$longitude == -167.625 &
#                                                             yukon.ice.may.mean$latitude== 62.375 &
#                                                             yukon.ice.may.mean$year >= 2000,],
#                       "year" = c(2000:2022))
full.df <- left_join(yukon.ice.may.mean[yukon.ice.may.mean$longitude == -168.625 &
                                          yukon.ice.may.mean$latitude== 63.125 &
                                          yukon.ice.may.mean$year >= 2000,], midpoint)
# Save the df for use in models
# saveRDS(object = full.df,
#         file = file.path(dir.data,"ICE Point -168625_63125 2000_2022.RDS"))

ggplot(full.df, aes(x = mean.ice, y = mid))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~`,`~")),label.x = 0)+
  labs(x = "SST (Lat = 63.125, Lon = -168.625)",
       y = "PSS Chinook Salmon Run Midpoint")


