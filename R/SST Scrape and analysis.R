# Aaron Lambert
# 2/21/2023
# Yukon Inseason Chinook Forecast
# SST 
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

# Midpoint from logistic curve fitting 1995-2022
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
require(ggmap)
require(maps)
require(usmap)
require(rnaturalearth)
require(rnaturalearthhires)
require(sf)
require(ggsn)
library(cowplot)
# The information for the NOAA OISST data
# rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")
# 
# # Note that there is also a version with lon values from 0 yo 360
# # rerddap::info(datasetid = "ncdcOisst21Agg", url = "https://coastwatch.pfeg.noaa.gov/erddap/")
# 
# # This function downloads and prepares data based on user provided start and end dates
# OISST_sub_dl <- function(time_df){
#   OISST_dat <- griddap(datasetx = "ncdcOisst21Agg_LonPM180",
#                        url = "https://coastwatch.pfeg.noaa.gov/erddap/",
#                        time = c(time_df$start, time_df$end),
#                        zlev = c(0, 0),
#                        # latitude = c(62,64),
#                        latitude = c(60,64),
#                        longitude = c(-163,-169),
#                        fields = "sst")$data %>%
#     mutate(time = as.Date(stringr::str_remove(time, "T00:00:00Z"))) %>%
#     dplyr::rename(t = time, temp = sst) %>%
#     select(longitude, latitude, t, temp) %>%
#     na.omit()
# }
# 
# # 
# # # Date download range by start and end dates per year
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
# # # Download all of the data with one nested request
# # # The time this takes will vary greatly based on connection speed
# system.time(
#   OISST_data <- dl_years %>%
#     group_by(date_index) %>%
#     group_modify(~OISST_sub_dl(.x)) %>%
#     ungroup() %>%
#     select(longitude, latitude, t, temp)
# ) # 38 seconds, ~8 seconds per batch

# saveRDS(OISST_data, file = file.path(dir.data,"OISST_data_163_169_62_64_21Sept23.RDS"))
# Visualize area selected
OISST_data <- readRDS(file = file.path(dir.data,"OISST_data_163_169_62_64_21Sept23.RDS"))

# Break out day, month , year , and jdate
yukon.sst  <- OISST_data %>% mutate(year = lubridate::year(t),
                              month = lubridate::month(t),
                              day = lubridate::day(t),
                              jday = lubridate::yday(t))

# Get Map info for plotting below
# Alaska
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


# Subset for May observations #########################################
yukon.sst.may <- yukon.sst[yukon.sst$month==5 & 
                             yukon.sst$year != 1996,]

yukon.sst.may.mean <- yukon.sst.may %>% 
  group_by(longitude, latitude, year) %>% 
  summarize(mean.sst = mean(temp)) 

yukon.sst.may.mean <- as.data.frame(yukon.sst.may.mean)

# Get the unique long and lat
lon <- unique(yukon.sst.may$longitude)
lat <- unique(yukon.sst.may$latitude)


# Plot out correlations for all years of PSS available ##############################
# Vector of midpoints
midpoint <- all.timing.logistic$mid

# Matrix to hold the values from correlation loop
cor.mat <- matrix(nrow = length(lat)*length(lon), ncol = 3)

colnames(cor.mat) <- c("lon", "lat", "cor")

counter <- 1
# i <- 15
# a <- 15
for (i in 1:length(lat)) {
  for (a in 1:length(lon)) {
    
    
    sst <- yukon.sst.may.mean$mean.sst[yukon.sst.may.mean$longitude == lon[a] & yukon.sst.may.mean$latitude==lat[i]]
    
    if(!is.na(sst[1]) ){
      
    
    cor.sst <- cor(sst,midpoint)
    
    
    
    cor.mat[counter,"lon"] <- lon[a]
    cor.mat[counter,"lat"] <- lat[i]
    cor.mat[counter,"cor"] <- cor.sst
    }else{
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- NA
    }
    
    counter <- counter +1
    
  }
  
}

cor.df <- data.frame(cor.mat)

map_1995_2022 <- 
  ggplot(data = ak) +
  geom_raster(data = cor.df, aes(x = lon, y = lat, fill = cor),interpolate = T) +
  #   # geom_sf(data = temp.dat)+
  # scale_fill_gradientn(colours = (rainbow(7)), na.value = NA, name = "r",limits = c(-.82,-.4)) +
  scale_fill_viridis_c(direction = -1, na.value = NA,name = "Correlation (r)")+
  geom_sf(fill = "darkgray") +
  geom_sf(data = ne50rivers)+# plot sf object
  geom_sf(data = ne10rivers_NA)+
  geom_sf(data = ne10rivers)+
  coord_sf(xlim = c(-170, -162), ylim = c(60, 64)) + # specify region-wide map coordinates (dimensions)
  # theme_bw() +
  guides(colour = guide_colourbar(barwidth  = unit(30,"cm")))+
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
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))

# Place inset map into main map
# pdf(file = file.path(dir.figs,"testmap.pdf"))
# png( file.path(dir.figs,"SST Cor Map 1995_2022 Plot 25Oct23.png"))

ggdraw() + 
  draw_plot(map_1995_2022) +
  draw_plot(inset.map, width = 0.2, height = 0.2, x = 0.715, y = 0.07)


# dev.off()
# Plot out correlations for subset of PSS years 2000 and onward ##############################

# Vector of midpoints
midpoint.trunc <- all.timing.logistic$mid[all.timing.logistic$year >= 2000]

# Matrix to hold the values from correlation loop
cor.mat <- matrix(nrow = length(lat)*length(lon), ncol = 3)

colnames(cor.mat) <- c("lon", "lat", "cor")

# Loop over sst observations to get correlation coeficients with midpoint
counter <- 1
# i <- 1
# a <- 1
for (i in 1:length(lat)) {
  for (a in 1:length(lon)) {
    
    
    sst <- yukon.sst.may.mean$mean.sst[yukon.sst.may.mean$longitude == lon[a] &
                                         yukon.sst.may.mean$latitude==lat[i] &
                                         yukon.sst.may.mean$year >= 2000]
    
    if(!is.na(sst[1]) ){
      
      
      cor.sst <- cor(sst,midpoint.trunc)
      
      
      
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- cor.sst
    }else{
      cor.mat[counter,"lon"] <- lon[a]
      cor.mat[counter,"lat"] <- lat[i]
      cor.mat[counter,"cor"] <- NA
    }
    
    counter <- counter +1
    
  }
  
}

# Convert matrix to df
cor.df.trunc <- data.frame(cor.mat)

# png(filename = file.path(dir.figs,"SST Cor Map 2000_2022 Plot 29Aug23.png"), width = 1000, height = 800)

# Plot the results
# ggplot(data = ak) +
#   geom_raster(data = cor.df.trunc, aes(x = lon, y = lat, fill = cor),interpolate = T) +
#   #   # geom_sf(data = temp.dat)+
#   scale_fill_gradientn(colours = (rainbow(7)), na.value = NA, limits = c(-.82,-.35)) +
#   ggtitle("Correlation Heat Map", subtitle = "May Mean SST vs Midpoint \n(2000-2022)")+# use ak map data
#   geom_sf() + 
#   geom_sf(data = ne50rivers)+# plot sf object
#   geom_sf(data = ne10rivers_NA)+
#   geom_sf(data = ne10rivers)+# plot sf object
#   coord_sf(xlim = c(-170, -162), ylim = c(60, 64)) + # specify region-wide map coordinates (dimensions)
#   theme_bw() +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(
#     # axis.text = element_text(size = 16, color = "black"),
#     text = element_text(size = 15),
#     axis.line = element_blank(),
#     panel.grid = element_blank(), # remove elements for clean background
#     plot.background = element_rect(fill = "transparent", colour = NA),
#     # legend.title=element_text(size=20),
#     # legend.text=element_text(size=14),
#     # legend.position = c(0.3, 0.3),
#     legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
#     legend.key = element_rect(fill = "white"))
map_2000_2022 <- 
  ggplot(data = ak) +
  geom_raster(data = cor.df.trunc, aes(x = lon, y = lat, fill = cor),interpolate = T) +
  #   # geom_sf(data = temp.dat)+
  # scale_fill_gradientn(colours = (rainbow(7)), na.value = NA, name = "r",limits = c(-.82,-.4)) +
  scale_fill_viridis_c(direction = -1, na.value = NA, name = "Correlation (r)")+
  geom_sf(fill = "darkgray") +
  geom_sf(data = ne50rivers)+# plot sf object
  geom_sf(data = ne10rivers_NA)+
  geom_sf(data = ne10rivers)+
  coord_sf(xlim = c(-170, -162), ylim = c(60, 64)) + # specify region-wide map coordinates (dimensions)
  # theme_bw() +
  # guides(colour = guide_colourbar(barwidth  = unit(50,"cm")))+
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
    # legend.position = c(0.3, 0.3),
    # legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
    legend.key = element_rect(fill = "white"),
    legend.position = "top", 
    # plot.margin = unit(c(0,0,0,0), "cm")
  )+
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5 ))

map_2000_2022

# Plot the subsetted correlations and the inset map
# png( file.path(dir.figs,"SST Cor Map 2000_2022 Plot 25Oct23.png"))

ggdraw() + 
  draw_plot(map_2000_2022) +
  draw_plot(inset.map, width = 0.2, height = 0.2, x = 0.715, y = 0.07)


# dev.off()
# dev.off()

# Look for strongest correlation
cor.df.trunc[which.min(cor.df.trunc$cor),]

# Subset by areas with high correlation
cor.df.trunc_0.7 <- cor.df.trunc[which(cor.df.trunc$cor <= -0.7),]

# Plot the subseted cells to take a look
ggplot(data = ak) +
  geom_raster(data = cor.df.trunc[which(cor.df.trunc$cor <= -0.7),], aes(x = lon, y = lat, fill = cor),interpolate = T) +
  #   # geom_sf(data = temp.dat)+
  scale_fill_gradientn(colours = (rainbow(7)), na.value = NA) +
  ggtitle("Correlation Heat Map", subtitle = "May Mean SST vs Midpoint \n(2000-2022)")+# use ak map data
  geom_sf() + # plot sf object
  geom_sf(data = ne50rivers)+# plot sf object
  geom_sf(data = ne10rivers_NA)+
  geom_sf(data = ne10rivers)+
  coord_sf(xlim = c(-170, -162), ylim = c(60, 64)) + # specify region-wide map coordinates (dimensions)
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    # axis.text = element_text(size = 16, color = "black"),
    text = element_text(size = 15),
    axis.line = element_blank(),
    panel.grid = element_blank(), # remove elements for clean background
    plot.background = element_rect(fill = "transparent", colour = NA),
    # legend.title=element_text(size=20),
    # legend.text=element_text(size=14),
    # legend.position = c(0.3, 0.3),
    legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
    legend.key = element_rect(fill = "white"))


# Data frame of mean sst at the location with the highest correlation when restricting to 2000-2022
full.df <- data.frame("mid" = midpoint.trunc,
                      "sst" = yukon.sst.may.mean$mean.sst[yukon.sst.may.mean$longitude == -164.875 &
                                                                                   yukon.sst.may.mean$latitude== 62.625 &
                                                                                   yukon.sst.may.mean$year >= 2000],
                      "year" = c(2000:2022))



ggplot(full.df, aes(y = mid, x = sst))+
  geom_point(size = 3)+
  geom_smooth(method = "lm",col = "black")+
  geom_text(aes(label = year), nudge_y = 0.3)+
  stat_cor(aes(label = paste(gsub("R", "r", ..r.label..),
                             gsub("R", "r", ..rr.label..),
                             ..p.label.., sep = "~`,`~")),label.x = 2, size = 5)+
  # coord_cartesian(ylim = c(152,185))+
  scale_y_continuous(labels = c("June 15",
                                "June 17",
                                "June 19",
                                "June 21",
                                "June 23",
                                "June 25",
                                "June 27",
                                "June 29"),
                     breaks = c(166,168,170,172,174,176,178,180))+
  labs(x = expression("May Mean SST"~ (degree*C)),
       y = "PSS Chinook Salmon Run Midpoint")+
  theme_cowplot()+
  theme(axis.title.y = element_text(vjust = 6), plot.margin = margin(5.5,5.5,5.5,15))




# Extract the mean sst from the specific lat lon when subsetted to <= -0.7
i=0
for(i in 1:length(cor.df.trunc_0.7[,1])){
  
  
LA <- as.numeric(cor.df.trunc_0.7$lat[i])
LO <- as.numeric(cor.df.trunc_0.7$lon[i])

if(i ==1){
  tempDF <-yukon.sst.may.mean %>% filter(longitude == LO,
                                         latitude == LA,
                                         year>=2000)
  tempDF2 <- tempDF
}else{
tempDF <-yukon.sst.may.mean %>% filter(longitude == LO,
                              latitude == LA,
                              year>=2000)

tempDF2<-rbind(tempDF2,tempDF)}

}

mean.tengrid <- tempDF2 %>% group_by(year) %>% summarize(mean.sst = mean(mean.sst))

# Data frame of mean sst at the location with correlations <= -.7 when restricting to 2000-2022
full.df.tengrid <- cbind(mean.tengrid, data.frame("mid" = midpoint.trunc))

# Save the df for use in models
# saveRDS(object = full.df,
#         file = file.path(dir.data,"SST Point 2000_2022.RDS"))

ggplot(full.df.tengrid, aes(x = mean.sst, y = mid))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(aes(label = paste(..r.label..,..rr.label..,..p.label.., sep = "~`,`~")),label.x = 2)+
  labs(x = "SST (averaged over 10 grids; r <= -0.7)",
       y = "PSS Chinook Salmon Run Midpoint")



# 
# sst <- yukon.sst.may.mean$mean.sst[yukon.sst.may.mean$longitude == lon[1] & yukon.sst.may.mean$latitude==lat[1]]
# 
# 
# 
# cor(sst,midpoint)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# norton.sst <- yukon.sst %>% 
#   group_by(year, month, day) %>% 
#   summarise(day.mean = mean(temp)) %>%
#   group_by(year,month) %>% 
#   summarise(month.mean = mean(day.mean))
# 
# 
# ggplot(norton.sst, aes(x = year, y = month.mean))+
#   geom_point()+
#   facet_wrap(~month)
# 
# 
# mid.point <- all.timing.logistic[,c(2,5)]
# 
# 
# norton.sst.mp <- left_join(norton.sst, mid.point)
# 
# norton.sst.mp <- norton.sst.mp[norton.sst.mp$year != 1996,]
# 
# ggplot(norton.sst.mp, aes(x = month.mean, y = mid))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~month, scales = "free_x")+
#   stat_cor(method = "pearson", label.y = 182)+
#   labs(x = "Mean SST Celsius",
#        y = "Observed Midpoint")+
#   theme(text = element_text(size = 18))
# 
# 
# may.sst <- norton.sst.mp[norton.sst.mp$month == 5,]
# 
# tt<-cor(x = na.omit(may.sst$month.mean), y = na.omit(may.sst$mid))
# 
# tt <- paste("r =",tt)
# 
# ggplot(may.sst, aes(x = month.mean, y = mid, label = (year)))+
#   geom_point(size = 3)+
#   geom_smooth(method = "lm")+
#   stat_cor(method = "pearson", label.y = 185)+
#   labs(x = "Mean SST Celsius",
#        y = "Observed Midpoint")+
#   geom_text(vjust = -1)+
#   scale_color_colorblind()+
#   theme(text = element_text(size = 18))
# 
# 
# 
