
# R script showing how the MCP 95% were calculated. In this script, it is also shown how we created a circle buffer to 
# select the nestboxes visited by each of our starlings. The resulting list of visited nestboxes for each male starlings was then used 
# to add the variables we were interested in. The script where we analyzed the relationship between the chosen nestbox and 
# the chosen variables is available in the repository.

##--Libraries--------------------------------------------------------------------
library(adehabitatHR)
library(sp)
library(tidyverse)
library(spatialEco)
library(sf)
library(geosphere)

##--Working directory------------------------------------------------------------
setwd("Y:/Iraida intentando instalar ARK/IMPORTANTE_TESIS/TESIS FLOTANTES/CAPITULO 2. PROSPECTING/02_Prospecting_in_the_starling/03_Analysis/03_Prospecting_Area/04_Focal_males")

##--Loading data-----------------------------------------------------------------

# Male dataset 
males <- read.csv("./data/all_males.csv", header = T, sep = ";", dec = ".") #dataset about male visits

# Nest site coordinates
nest_coordinates <- read.csv("./data/GPS.csv", header=T, dec=",", sep=";")  
nest_coordinates <- nest_coordinates %>% rename(nest_dummy = CAJA)
nest_coordinates <- nest_coordinates %>% rename(nest = CAJA_REAL)


# Modifying part of the dataset to avoid errors: replacing names of nest boxes that contain letters.
males <- males %>% mutate("nest_dummy" = nest) # nestbox E1 is now nestbox 1000
males$date <- as.Date(males$date, format = "%d/%m/%Y") # date format


# Changing names containing letters.
males$nest_dummy[which(males$nest_dummy=="E6")] <- 6000
males$nest_dummy[which(males$nest_dummy=="E5")] <- 5000
males$nest_dummy[which(males$nest_dummy=="E4")] <- 4000
males$nest_dummy[which(males$nest_dummy=="E3")] <- 3000
males$nest_dummy[which(males$nest_dummy=="E2")] <- 2000
males$nest_dummy[which(males$nest_dummy=="E1")] <- 1000
males$nest_dummy <- as.numeric(males$nest_dummy)
males <- males %>% rename(nidos_visitados = Nidos.diferentes) 

# Relocation of one column 
males <- males %>% relocate(nest_dummy, .before = nest)

##  Joining nestbox coordinates to dataset ----
nest_coordinates_dummy <- nest_coordinates[,-c(4:6)]
males <- left_join(males, nest_coordinates_dummy, by = "nest_dummy") 

# Latitude and longitud as numeric 
males$LAT <- as.numeric(males$LAT)
males$LON <- as.numeric(males$LON)

# Dummy code to identify unique 'ring-nest'.
males <- males %>% mutate("ring_nest_dummy" = paste0(ring, nest_dummy)) 

# Duplicate of the dataset
starlings <- males


## Computation of the 95% minimum convex polygon (MCP)

nest_coordinates.sp <- nest_coordinates[, c("nest_dummy", "LON", "LAT")] 

# Numeric coordinates
nest_coordinates.sp$LON <- as.numeric(nest_coordinates.sp$LON) 
nest_coordinates.sp$LAT <- as.numeric(nest_coordinates.sp$LAT)

# Making an spatial object and adjusting projection of nest boxes.
nest_coordinates.spdf <- SpatialPointsDataFrame(coords=(cbind(nest_coordinates.sp$LON, nest_coordinates.sp$LAT)), 
                                                data = nest_coordinates.sp["nest_dummy"], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
nest_coordinates.spdf <- spTransform(nest_coordinates.spdf, CRS("+proj=utm +zone=30 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


# Selecting data to do the MCP 95%.
starlings.sp <- starlings[, c("ring", "LON", "LAT")] 

# Making an spatial object and adjusting projection of nest boxes that have been visited.
starlings.spdf <- SpatialPointsDataFrame(coords=(cbind(starlings.sp$LON, starlings.sp$LAT)), 
                                         data = starlings.sp["ring"], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
starlings.spdf <- spTransform(starlings.spdf, CRS("+proj=utm +zone=30 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))



# Calculating the MCP. 
starlings.mcp <- mcp(starlings.spdf, percent = 95) 
areas_MCP <- starlings.mcp@data # sizes of polygons

# We use th median MCP.
mediana_areaMCP <- median(areas_MCP$area)
sd_areMCP <- sd(areas_MCP$area)
radio_areaMCP <-sqrt((mediana_areaMCP*10000)/pi)
sd_en_metros <-sqrt((sd_areMCP*10000)/pi)
radio_en_grados<-(radio_areaMCP*1/111320)

mean(areas_MCP$area)
range(areas_MCP$area)

# Getting the centroid point of each starling
mcp_centroid <- gCentroid(starlings.spdf, byid=T)
centroid_coords <- coordinates(mcp_centroid)
centroides <- as.data.frame(centroid_coords)
centroides$males <- rownames(centroid_coords)


# We calculate the radius (in meters) of a circle of the same area and make a buffer to find the nestboxes that fall within. 

# Load the dataset containing the centroids.
centroids <- read.csv("./data/all_centroids.csv", sep = ";", head = T)

# Centroids -> proyeccion
centroids.spdf <- SpatialPointsDataFrame(coords=(cbind(centroids$x, centroids$y)), 
                                         data = centroids["male_ring"], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
centroids.spdf <- spTransform(centroids.spdf, CRS("+proj=utm +zone=30 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))


# Finding the distance between nestboxes
distancias <- distm(centroids[,2:3], nest_coordinates.sp[, 2:3], fun = distHaversine)

rownames(distancias) <- centroids$male_ring
colnames(distancias) <- nest_coordinates.sp$nest_dummy

distancias_df <- as.data.frame(distancias)
distancias_df$centroids_male <- rownames(distancias_df)
distancias_df <- distancias_df %>% dplyr::relocate(centroids_male, .before="1")
rownames(distancias_df) <- NULL


# Some organization of the dataset. 

data <- gather(data = distancias_df, key = "nest_dummy", value = "distance", 2:252)
data <- data %>% dplyr::filter(distance <= radio_areaMCP)
data <- data %>% rename(ring = "centroids_male")

nests_within <- data %>% mutate("ring_nest_dummy" = paste0(ring, nest_dummy)) 
cajas_visitadas <- nests_within %>% dplyr::inner_join(males[, c("evento", "ring_nest_dummy")], by= "ring_nest_dummy")

# We create a new column that codes if a given starling was detected or not visiting that nestbox in the recording database.
# 
nests_within <- nests_within %>% mutate("visited" = ifelse(nests_within$ring_nest_dummy %in% cajas_visitadas$ring_nest_dummy, 1, 0))

# Adding year of visiting
nests_within$año_recaptura <- males$año_recaptura[match(nests_within$ring,males$ring)]
nests_within$año_recaptura <- males$año_recaptura[match(nests_within$ring,males$ring)]

# New coding to identify each nest per year. ----
nests_within <- nests_within %>% mutate("YEAR_NEST" = paste(año_recaptura, nest_dummy, sep = "_"))