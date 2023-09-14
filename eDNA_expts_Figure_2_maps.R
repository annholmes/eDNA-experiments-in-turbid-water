#Maps for turbidity paper

library(sf) #imports shapefiles
library(ggplot2)
library(tmap) #not used?
library(tmaptools) #not used?
library(leaflet) #not used?
library(dplyr) #not used?
library(ggspatial) # for scale bar and north symbol
library(cowplot) #used for inset maps and arranging multiple plots

options(scipen = 999) #so numbers don't display as scientific notation

#read in shapefiles and data
setwd("/Users/aholmes/Desktop/github/eDNA-expts")

#read shapefiles
#download and unzip data for map from https://geodata.lib.berkeley.edu/catalog/stanford-qh320kj0191
SFE <- st_read("/Users/aholmes/Desktop/github/eDNA-expts/data/water_bodies_carto.shp", 
               stringsAsFactors = FALSE) #all other files must remain in the unzipped directory in order to read the shape file

#shape files for inset maps showing where sampling is located within region
#download and unzip data for map from https://geodata.lib.berkeley.edu/catalog/stanford-vv027px1909
CA_coast <- st_read("/Users/aholmes/Desktop/github/eDNA-expts/stanford-vv027px1909-shapefile/vv027px1909.shp", 
                    stringsAsFactors = FALSE) #all other files must remain in the unzipped directory in order to read the shape file
NAmer <- st_read("/Users/aholmes/Desktop/github/eDNA-expts/stanford-ns372xw1938-shapefile/ns372xw1938.shp", 
                 stringsAsFactors = FALSE) #all other files must remain in the unzipped directory in order to read the shape file

#read data
#file with lat and lon for all EDSM detections from 15 Dec 2016 to 31 Mar 2017
#includes turbidity and secchi values from the original EDSM metadata from USFWS
#data has been cleaned up so that each row is 1 tow rather than each row is 1 or more fish 
EDSM <- read.csv("EDSM_KDTR_1617_turbidity.csv")

#a map for eDNA detections only
eDNA_EDSM <- read.csv("eDNA_EDSM_stations_detections.csv")

#map for sites where water was collected for the filtration experiments
expts_water <- read.csv("Expts_in_turbid_water_map_pts.csv")

#read in point to show city of Sacramento
Sac <- read.csv("Sacramento.csv")

#data exploration
#just take a look at the shapefile
plot(st_geometry(SFE))
#this makes a plain ole map using the SFE shapefile
#it's colored blue, cropped to the study area, and has defined axis labels
ggplot() + 
  geom_sf(data = SFE, fill= "#89CFEF", color = "#89CFEF") + #fill and outline same color
  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and lon
               xlim = c(-122.3,-121.4), ylim = c(37.9,38.6)) + #set limits of map
  scale_x_discrete(breaks = c(-122.2, -122.0, -121.8, -121.6, -121.4), 
                   labels = c("-122.2", "-122.0", "-121.8", "-121.6", "-121.4")) + #label x axis
  scale_y_discrete(breaks = c(38.0, 38.2, 38.4, 38.6), 
                   labels = c("38.0", "38.2", "38.4", "38.6")) #label y axis

#to add the EDSM tow sampling points, first define the tow coords to be used as points
#the sampling metadata lists target lat/lon for the station and the actual start and stop lat/lon for each tow
#here I use the actual start lat and start lon recorded for each tow
sites <- st_as_sf(EDSM, 
                  coords = c("StartLong", "StartLat"), 
                  crs = 4326, 
                  agr = "constant")

#define the eDNA sites; these are the 6 (of 32) paired eDNA-EDSM sites where DSM were detected in the trawl (year 2017)
eDNA_sites <- st_as_sf(eDNA_EDSM, 
                       coords = c("lon", "lat"), 
                       crs = 4326, 
                       agr = "constant")

#sites where water was collected for the experiments
water_sites <- st_as_sf(expts_water, 
                        coords = c("lon", "lat"), 
                        crs = 4326, 
                        agr = "constant")

#definte a point for the city of Sacramento
Sac_site <- st_as_sf(Sac, 
                     coords = c("lon", "lat"), 
                     crs = 4326, 
                     agr = "constant")

#add EDSM tow sampling data to the SFE map
map <- ggplot(data = SFE) +
  theme(panel.grid.major = element_blank(),#remove grid lines
        legend.position = c(0.87, 0.53)) + #set position of legend on figure
  geom_sf(fill= "#bbe4ee", color = "#bbe4ee") + #make water features and outline of water features blue
  geom_sf(data = sites, aes(color=Turbidity), size = 3) +
  geom_sf(data = Sac_site, size = 3, shape = 0) +
  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and long
           xlim = c(-122.4,-121.3), ylim = c(37.9,38.65)) + #set limits of map
  scale_x_discrete(breaks = c(-122.2, -122.0, -121.8, -121.6, -121.4), 
                   labels = c("-122.2", "-122.0", "-121.8", "-121.6", "-121.4")) + #label x axis
  scale_y_discrete(breaks = c(38.0, 38.2, 38.4, 38.6), 
                   labels = c("38.0", "38.2", "38.4", "38.6")) #label y axis
  
map

#annotate map for EDSM turbidity values, scale and north symbol
EDSM_map <- map + scale_color_gradient2(low = "#3366CC", #low end of turbidity gradient is blue
                           mid = "#663500", #midpoint of tubidity is brown
                           high = "#CC0000", #high end of turbidity gradient is red
                           midpoint = 70, #set the midpoint
                           space = "Lab",
                           name = bquote("Turbidity (NTU)"), #add the units to the legend
                           na.value = "#CCCCCC") + #tows that didn't record a turbidity value are gray
  annotation_scale(location = 'bl', #add scale bar in bottom left
                   height = unit(0.4, "cm"), #increase height of scale bar
                   text_cex = 1) + #increase relative size of text
  annotation_north_arrow(location = 'tr') + #add north symbol
  annotate("text", x=-121.55, y=38.605, label= "Sacramento", size = 6) +
  #annotate("text", x=-121.4, y=37.95, label= "A", size = 6) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(colour = "black"))

EDSM_map

#make an inset map of the entire SFE showing the sampling sites in black
inset_SFE <- ggplot(data = CA_coast) + 
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_sf(color = "#999DA0") + 
  geom_sf(data = sites, size = 1) +
  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and lon
           xlim = c(-123.5,-121.0), ylim = c(37.0, 39.0)) 

inset_SFE

#another option: land is shaded but there's less detail in water features
#inset_SFE_option <- ggplot(data = NAmer) + 
#  theme(panel.grid.major = element_blank(),
#        axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank(),
#        axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank()) +
#  geom_sf(color = "#999DA0") + 
#  geom_sf(data = sites, size = 1) +
#  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and lon
#           xlim = c(-123.5,-121.0), ylim = c(37.0, 39.0)) 

#make an inset map of the US west coast with a box around the SFE
inset_CA <- ggplot(data = NAmer) + 
  theme(panel.grid.major = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_sf(color = "#999DA0") + 
  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and lon
           xlim = c(-125,-110), ylim = c(30, 50)) + 
  annotate("rect", xmin=c(-123.5), xmax=c(-121.0), 
           ymin=c(37.0) , ymax=c(39.0), 
           alpha=0.2, color="black") #add box

#the west coast map will go in the top left corner, add that first
EDSM_CA <- ggdraw(EDSM_map) +
  draw_plot(
    {
      inset_CA +
        theme(legend.position = "none")
    },
    x = 0.06, #x-axis location for left edge of plot
    y = .62, #y-axis location for bottom edge of plot
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.33, 
    height = 0.33)

EDSM_CA

#then add the SFE regional map
EDSM_CA_SFE <- ggdraw(EDSM_CA) +
  draw_plot(
    {
      inset_SFE +
        theme(legend.position = "none")
    },
    x = 0.26, #x-axis location for left edge of plot
    y = .62, #y-axis location for bottom edge of plot
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.33, 
    height = 0.33)

#final map of DSM detections by EDSM with insets
#x,y location of inset maps adjusted by trail and error- don't spend too much time of this because the placement changes slightly when the plots are arranged at the end
EDSM_CA_SFE

#for eDNA collections paired with EDSM
eDNA_map <- ggplot(data = SFE) +
  theme(panel.grid.major = element_blank(),#remove grid lines
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  geom_sf(fill= "#bbe4ee", color = "#bbe4ee") + #make water features and outline of water features blue
  geom_sf(data = eDNA_sites, aes(color=Turbidity), size = 3) +
  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and long
           xlim = c(-122.2,-121.65), ylim = c(38.0, 38.25)) +  #set limits of map
  scale_color_gradient2(low = "#3366CC", #low end of turbidity gradient is blue
                        mid = "#663500", #midpoint of tubidity is brown
                        high = "#CC0000", #high end of turbidity gradient is red
                        midpoint = 70,
                        space = "Lab",
                        na.value = "#CCCCCC") +
  #annotate("text", x=-121.64, y=37.98, label= "B", size = 6) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

#map of eDNA detections paired with EDSM
eDNA_map

#for water sites used in filtration experiment
water_map <- ggplot(data = SFE) +
  theme(panel.grid.major = element_blank(),#remove grid lines
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  geom_sf(fill= "#bbe4ee", color = "#bbe4ee") + #make water features and outline of water features blue
  geom_sf(data = water_sites, aes(color=Turbidity), size = 3) +
  geom_sf(data = Sac_site, size = 3, shape = 0) +
  coord_sf(default_crs = sf::st_crs(4326), #x and y positions interpreted as lat and long
           xlim = c(-121.8,-121.3), ylim = c(38.45, 38.65)) +  #set limits of map
  scale_color_gradient2(low = "#3366CC", #low end of turbidity gradient is blue
                        mid = "#663500", #midpoint of tubidity is brown
                        high = "#CC0000", #high end of turbidity gradient is red
                        midpoint = 70,
                        space = "Lab",
                        na.value = "#CCCCCC") +
  annotate("text", x=-121.54, y=38.6, label= "Sacramento", size = 6) +
  #annotate("text", x=-121.34, y=38.62, label= "C", size = 6) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

#map of water collection sites
water_map 

#arrange 3 plots and label A, B, C
#first put two smaller plots in one column
#modify the margin on the left side so the plots will sit closer to the main plot
eDNA_map1 <- eDNA_map + theme(plot.margin=margin(l=-0.8,r=0.5, unit="cm"))
water_map1 <- water_map + theme(plot.margin=margin(l=-0.8,r=0.5, unit="cm"))
right_col <- plot_grid(water_map1, eDNA_map1, 
                       labels = c('C', 'B'), 
                       ncol = 1,
                       label_size = 12)
#then add larger map on the left
final_map <- plot_grid(EDSM_CA_SFE, right_col,
          labels = c('A','',''), 
          ncol = 2,
          label_size = 12,
          rel_heights = c(2,1),
          rel_widths = c(2,0.9))

#recommended file type png, size 3000px (min 900 px)
#conversion at 300 dpi is 10in/25.4cm and 3in/7.62cm
save_plot("plot1k.png", 
          final_map, 
          dpi = 300,
          base_height = 6, #height in inches, equal to 15.24cm
          base_width = 12, #height in inches, equal to 30.48cm
          bg = "white")

