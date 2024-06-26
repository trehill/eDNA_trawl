#Map of study site ####
#https://jtr13.github.io/cc19/using-stamen-maps-for-plotting-spatial-data.html

#Set-up ####
#install.packages(ggmap)

library(tidyr)
library(tidyverse)
library(ggplot2)
library(here)
library(maps)
library(dplyr)
library(ggmap)
library(maps)


#read data 
metatrawl <- read.csv(here::here("Processed_data", 
                                 "trawl",
                                 "metadata", 
                                 "trawl_metadata.csv"),
                      head=TRUE)


metaeDNA <- read.csv(here::here("Processed_data", 
                                "eDNA",
                                "metadata", 
                                "eDNA_metadata.csv"),
                     head=TRUE)

map <- merge(metatrawl, metaeDNA, by=c('set_number'))
map <- subset(map, !is.na(depth))
map <- select(map, c('depth', 'lat_door_in_dd','long_door_in_dd', 'set_number'))
map <- distinct(map)
map <- subset(map, depth != 5) #map is our dataset 

#register API to use ggmaps
register_google('AIzaSyBG4StmTBuRubdtKeUZWPt-1E1F08uWc-U') #register API 


#plot in color
map.for.samples <- get_map(location = c(-128,47.5,-122,51.5),
                           source = 'google',
                           maptype='terrain',
                           api_key = 'AIzaSyBG4StmTBuRubdtKeUZWPt-1E1F08uWc-U') # <- REPLACE WITH YOUR KEY


edna_sample_map<-ggmap(map.for.samples) +
  geom_point(data = map,
             aes(x = long_door_in_dd, y = lat_door_in_dd,
                 colour=depth), size=2) +
  scale_colour_gradient(name='depth of sample', low = '#00AFBB', high = '#006a71') #set colours so dark is deeper

edna_sample_map

#plot in b&w
map.for.samples <- get_map(location = c(-128,47.5,-122,51.5),
                           maptype = 'toner-background', #change this to terrain for terrain background
                           source = c("stamen"),
                           color=c('bw'),
                           api_key = 'AIzaSyBG4StmTBuRubdtKeUZWPt-1E1F08uWc-U') # <- REPLACE WITH YOUR KEY

# change opacity of basemap
mapatt <- attributes(map.for.samples)
map_transparent <- matrix(adjustcolor(map.for.samples, alpha.f = 0.2), nrow = nrow(map.for.samples))
attributes(map_transparent) <- mapatt

map.for.samples <- ggmap(map_transparent) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white'))+
  geom_point(data = map,
             aes(x = long_door_in_dd, y = lat_door_in_dd,
                 colour=depth), size=2) +
  scale_colour_gradient(name='depth of sample', low = '#00AFBB', high = '#006a71') + #set colours so dark is deeper
  ylab(c("latitude")) + xlab(c("longitude")) 

map.for.samples


ggsave("./Outputs/metadata/mapbw.png", 
       plot = map.for.samples,
       width = 5, height = 5, units = "in")

#adding labels to the points for each set number
map.for.samples <- get_map(location = c(-128,47.5,-122,51.5),
                           maptype = 'toner-background', #change this to terrain for terrain background
                           source = c("stamen"),
                           color = c('bw'),
                           api_key = 'AIzaSyBG4StmTBuRubdtKeUZWPt-1E1F08uWc-U') # <- REPLACE WITH YOUR KEY

# change opacity of basemap
mapatt <- attributes(map.for.samples)
map_transparent <- matrix(adjustcolor(map.for.samples, alpha.f = 0.2), nrow = nrow(map.for.samples))
attributes(map_transparent) <- mapatt

map.for.samples <- ggmap(map_transparent) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white'))+
  geom_text(data = map,
            aes(x = long_door_in_dd, y = lat_door_in_dd, label = set_number, fontface = "bold"),
            size = 4, position = position_jitter(width = 0.03, height = 0.03)) +
  ylab("latitude") + xlab("longitude") 

map.for.samples

ggsave("./Outputs/metadata/samples.png", 
       plot = map.for.samples,
       width = 5, height = 5, units = "in")



#Map of Canada
map.canada <- map_data("world", region = "Canada")

p1 <- ggplot(map.canada, aes(x = long, y = lat, group = group))
p1 <- p1 + geom_polygon(fill = "grey", color = "white") # Set fill to white, color to grey, and fill.alpha to 1
p1 <- p1 + theme_void() # Use a void theme for a clean background


p1 # Display the map of Canada

ggsave("./Outputs/metadata/canada.png", 
       plot = p1,
       width = 5, height = 5, units = "in")

#Changing colors 
map.for.samples <- get_map(location = c(-128,48,-123,51),
                           maptype = 'toner-background', #change this to terrain for terrain background
                           source = c("stamen"),
                           color = c('bw'),
                           api_key = 'AIzaSyBG4StmTBuRubdtKeUZWPt-1E1F08uWc-U') # <- REPLACE WITH YOUR KEY


ggmap::ggmap(map.for.samples)

# invert colors in raster
invert <- function(x) rgb(t(255-col2rgb(x))/255)    
m_inv <- as.raster(apply(map.for.samples, 2, invert))

# copy attributes from original object
class(m_inv) <- class(map.for.samples)
attr(m_inv, "bb") <- attr(map.for.samples, "bb")

ggmap(m_inv)

map.for.samples <- m_inv

# change opacity of basemap
mapatt <- attributes(map.for.samples)
map_transparent <- matrix(adjustcolor(map.for.samples, alpha.f = 0.2), nrow = nrow(map.for.samples))
attributes(map_transparent) <- mapatt

ggmap(map_transparent)

map.for.samples <- ggmap(map_transparent) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white'))+
  geom_text(data = map,
            aes(x = long_door_in_dd, y = lat_door_in_dd, label = set_number, fontface = "bold"),
            size = 5, position = position_jitter(width = 0.045, height = 0.045)) +
  ylab("latitude") + xlab("longitude") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = 'white'),
        axis.text = element_text(size = 10))  # Adjust the size of axis text

map.for.samples

ggsave("./Outputs/metadata/samples.png", 
       plot = map.for.samples,
       width = 15, height = 10, units = "in")

#### Determining smallest distance between two sites
# Install and load the geosphere package if not already installed
if (!requireNamespace("geosphere", quietly = TRUE)) {
  install.packages("geosphere")
}
library(geosphere)

map <- map %>%
  rename(set_number = set_number,
         lat = lat_door_in_dd,
         lon = long_door_in_dd)

# Assuming you already have your 'map' dataset loaded

# Function to calculate distance between two points using Haversine formula
calculate_distance <- function(lat1, lon1, lat2, lon2) {
  dist <- distHaversine(c(lon1, lat1), c(lon2, lat2))
  return(dist)
}

# Function to find smallest distance between two set numbers
find_smallest_distance <- function(map_data, set_number1, set_number2) {
  # Filter data for the two set numbers
  data1 <- subset(map_data, set_number == set_number1)
  data2 <- subset(map_data, set_number == set_number2)
  
  # Initialize smallest distance
  smallest_distance <- Inf
  
  # Iterate through all combinations of points and find smallest distance
  for (i in 1:nrow(data1)) {
    for (j in 1:nrow(data2)) {
      dist <- calculate_distance(data1[i, "lat"], data1[i, "lon"],
                                 data2[j, "lat"], data2[j, "lon"])
      if (dist < smallest_distance) {
        smallest_distance <- dist
      }
    }
  }
  
  return(smallest_distance)
}

# Example usage:
# Assuming 'map' is your dataset
set_number1 <- 1
set_number2 <- 2
smallest_dist <- find_smallest_distance(map, set_number1, set_number2)
cat("Smallest distance between set number", set_number1, "and", set_number2, "is:", smallest_dist, "meters\n")


