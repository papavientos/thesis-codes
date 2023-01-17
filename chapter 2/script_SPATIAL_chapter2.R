# Spatial autocorrelation

#' One of the methods to calculate spatial autocorrelation is the Moran's I.
#' The Moran's I quantifies how similar an object is with others around it,
#' allowing to check if different variables are 'grouped' in the space.
#' I am not sure if this is important at this scale (although there are interesting studies on this matter)
#' I think next steps would require us to make comparisons among colonies situated farther away.
#' Giving the organisation of our colony and the biology of the starling, I would not expect to 
#' find clusters of higher or lower reproductive succcess. 
#' Long time ago I asked at ResearchGate ((https://www.researchgate.net/post/How_to_compute_Morans_I_for_multiple_years) 
#' which would be the best way to calculate the Moran's I for multiple years. 
#' I was recommended to calculate the mean reproductive success or make correlograms for each year (https://www.researchgate.net/post/How_to_compute_Morans_I_for_multiple_years?_ec=topicPostOverviewAuthoredQuestions&_sg=fth9pz1_AztSvV76apyluZSDQBhEDmXMUDuI8eV06xMhPWkpspCb8TwP4IQBy6ei6BOIoNOUFRdjq5Ms). 
#' In this correlograms, distance can be treated as a discrete variable (making distance bands or as
#' a continous variable.
#' I also contacted Anne Charmantier to ask for the method they used in one of their publications
#' (https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12448).
#' Marrot, coauthor, told me this: "Indeed, the spanning tree procedure is particularly useful to determine the minimum distance lag to keep all the nest boxes linked, 
#' thereby producing the “best” correlogram in term of precision"
  
#--Libraries
library(sp)
library(spdep) # for the correlogram.
library(splm)
library(ape)
library(fields) # for computing distance matrix between lat/lon points.
library(arm)
library(geoR)
library(vegan)  # for computing spantrees.
library(geosphere)
library(tidyverse)

#--Data
nests <- read.csv("./data/GPS.csv", sep = ";", head = T, dec = ".") # coordinates of nest boxes
rs <- read.csv("./data/exito_rep.csv", sep = ",", head = T, dec = ".") #reproductive success data

# Customizing reproductive success dataframe
rs <- rs[,-1]
rs <- rs %>% filter(! YEAR_NEST == "#N/D") # Deleting no valid observations
rs <- rs %>% rename(first_broods = X1, second_broods = X2) 
rs <- rs %>% mutate(total_vols = first_broods + second_broods) # new variable
rs_join <- rs %>% dplyr::select(1, 4, 7) 
rs_join <- rs_join %>% dplyr::distinct(YEAR_NEST, .keep_all = T) # Remove duplicates


# Nest data
nests_data <- array(0, dim = c(nrow(nests), ncol = 10)) # Empty array
nests_data <- as.data.frame(nests_data) 
nests_data$CAJA <- nests$CAJA 
nests_data$LON <- nests$LON 
nests_data$LAT <- nests$LAT 
nests_data <- nests_data %>% relocate(CAJA, .before = everything()) 
colnames(nests_data)[2:11] <- c(2012:2021) 


nest_data_longer <- nests_data %>% pivot_longer(!c(CAJA, LON, LAT), names_to = "year", values_to = "reproductive_success")
nest_data_longer <- nest_data_longer %>% mutate(year_nest = paste(year, CAJA, sep = "_"))
nest_data_longer <- nest_data_longer %>% mutate(bred = ifelse(nest_data_longer$year_nest %in% rs $YEAR_NEST, 1, 0))

# Together
nestdata_joined_rs <- nest_data_longer %>% left_join(rs_join, by = c("year_nest" = "YEAR_NEST"))
# Nest boxes where there was not breeding event are assigned a 0
nestdata_joined_rs$total_vols[nestdata_joined_rs$bred==0] <- 0 
nestdata_joined_rs <- nestdata_joined_rs %>% group_by(CAJA) %>% mutate(mean = mean(total_vols))
mean.rs <- nestdata_joined_rs %>% dplyr::select(1:3, 10)
mean.rs <- mean.rs %>% distinct(CAJA, .keep_all = T)


# Renaming
geo = mean.rs               

# Making it a spatial objetc
coordinates(geo) <- ~LON+LAT                  

XYstarling = coordinates(geo)
# Calculating distances
xy.dl = rdist.earth(XYstarling, miles = FALSE)
diag(xy.dl) <- 0                               
# Calculating the spanning tree
spanning = spantree(xy.dl)   

# From this, we obtain the lag distance
dmin = max(spanning$dist)
# And the maximun distance detected
dmax = max(xy.dl)         
# Checking coordinates
plot(coordinates(geo))    


# Calculation (Marrot et al. 2017)
nb1 = dnearneigh(as.matrix(XYstarling), 0, dmin, longlat = TRUE) 
correloresm = sp.correlogram(nb1, mean.rs$mean, order = (dmax/dmin), method = "I", zero.policy = TRUE) 
plot(correloresm, ylim = c(-0.5, 0.5), main = "2012-2021 Mean Reproductive Success")

# Dataframe for ggplot
neigh <- data.frame(cor = correloresm$res[,1], lag = c(1:8), upper = 2*(sqrt(correloresm$res[,3])), 
                    lower = (-1)*(2*(sqrt(correloresm$res[,3]))))

# Representation
neigh1 <- ggplot(neigh, aes(x = lag, y = cor))+
  geom_linerange(aes(ymin = (cor+lower), ymax = (cor+upper)), lwd = 0.8) +
  geom_point(size = 4, col = "purple") +
  theme(axis.title.x= element_text(face="plain", size=16,  margin = margin(t=20, r=0, l=0, b=0)),
        axis.title.y= element_text(face="plain", size=16, margin = margin(t=0, r=20, l=0, b=0)),
        axis.text.x = element_text(face="plain", size =12, vjust = 1, hjust = 1),
        axis.text.y = element_text(face="plain", size =12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right", 
        legend.title = element_text(colour="black", size=14, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.3, linetype = "solid"),
        panel.spacing.y = unit(2, "lines"))+
  geom_hline(yintercept = 0, lwd = 0.5, linetype = "dashed")+
  scale_x_continuous(breaks=seq(0,8,1))+
  labs(y="Moran's I", x="Lag") +
  #annotate("text", x = 1.5, y = 0.48, label = diference, parse = T) +
  #annotate("segment", x = 1, xend = 2, y = 0.51, yend = 0.51) + 
  scale_y_continuous(limits = c(-0.20, 0.20), breaks = seq(-0.20, 0.20, by = 0.05)) 

plot(neigh1)
