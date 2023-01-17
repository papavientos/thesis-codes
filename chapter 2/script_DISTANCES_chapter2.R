##--Libraries--------------------------------------------------------------------
library(adehabitatHR)
library(sp)
library(tidyverse)
library(spatialEco)
library(sf)
library(geosphere)


# Nest boxes where focal individuals ended up breeding 
males <- read.csv("./data/reproduction_first_nest.csv", header = T, sep = ";", dec = ".") #dataset de machos de un a?o e
males <- males %>% dplyr::select(2, 5, 7, 6)
males <- males %>% arrange(male_ring)

# Centroids of the prospecting area of focal individuals.
centroids <- read.csv("./data/coords_centroids.csv", sep = ";", head = T)
centroids <- centroids[,-4]
centroids <- centroids %>% relocate(LON_centroid, .before = LAT_centroid)
centroids <- centroids %>% arrange(male_ring)

# Coordinates 

all_nest <- read.csv("./data/GPS.csv", sep = ";", head = T)
all_nest <- all_nest %>% dplyr::select(1, 2, 3)


# Distances between nestboxes and centroids

real_distances <- read.csv("./data/distancia_centroid_caja_rep.csv", sep = ",", head = T)
real_distances <- real_distances %>% arrange(male)
mean(real_distances$distancias_diag)
sd(real_distances$distancias_diag)
boxplot(real_distances$distancias_diag)
range(real_distances$distancias_diag)


# Distance between centroids and any nest of the colony.

shuffle <- all_nest
p <- array(0, dim = c(10000, 1))

for(i in 1:nrow(p)){
  shuffle_loop <- slice_sample(shuffle, n = 139, replace = T)
  distancias_loop <- distm(centroids[,1:2], shuffle_loop[, 3:2], fun = distHaversine)
  rownames(distancias_loop) <- centroids$male_ring
  colnames(distancias_loop) <- shuffle_loop$CAJA
  distancias_loop_diag <- diag(distancias_loop)
  distancias_loop_diag_df <- as.data.frame(distancias_loop_diag)
  mean <- mean(distancias_loop_diag_df$distancias_loop_diag)
  p[i,1] <- mean
}

# Representation
hist(p, xlim = c(0, 600))
abline(v = mean(real_distances$distancias_diag), col = "red")

# Put it in a dataframe
tipo_p <- rep("boots", 10000)
data_p <- cbind(p, tipo_p)


# Bootstrapping of the dataset
shuffle_statquest <- real_distances[,-1]
s <- array(0, dim = c(10000, 1))
for(i in 1:nrow(s)){
  shuffle_loop <- slice_sample(shuffle_statquest, n = 139, replace = T)
  mean <- mean(shuffle_loop$distancias_diag)
  s[i,1] <- mean
}

# Representation
hist(s, xlim = c(0, 600))
abline(v = mean(real_distances$distancias_diag), col = "red")

# Put it in a dataset
tipo_s <- rep("real", 10000)
data_s <- cbind(s, tipo_s)

# Unifying both datasets
all <- rbind(data_p, data_s)
all <- as.data.frame(all)
all$tipo_p <- as.factor(all$tipo_p)
all$V1 <- as.numeric(all$V1)
boxplot(V1 ~ tipo_p, data = all)


data_p <- as.data.frame(data_p)
data_p$V1 <- as.numeric(data_p$V1)
median(data_p$V1)

data_s <- as.data.frame(data_s)
data_s$V1 <- as.numeric(data_s$V1)
median(data_s$V1)

under_random <- sum(data_p$V1 < 168.57)
over_random <- sum(data_p$V1 > 168.57)
p_nac_flot_random <- (under_random/over_random)

layout(mat = matrix(c(1,2,1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(data_p$V1 , horizontal=TRUE , ylim=c(100, 500), xaxt="n" , col="#7fcdbb" , frame=F)
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(data_s$V1 , horizontal=TRUE , ylim=c(100, 500), xaxt="n" , col="#756bb1" , frame=F, add =T)
par(mar=c(4, 3.1, 1.1, 2.1)) 
hist(data_p$V1 , breaks=40 , col="#7fcdbb" , border=F , main="" , xlab="Distance (m)", xlim=c(100,500))
par(mar=c(4, 3.1, 1.1, 2.1)) 
hist(data_s$V1 , breaks=40 , col= "#756bb1", border=F , main="" , xlab="Distance (m)", xlim=c(100,500), add = T)
abline(v = mean(real_distances$distancias_diag), col = "black", lwd = 3, lty = "dashed", add = T)
legend(170, 1200, legend=c("Bootstrapped dataset", "Bootstrapped dataset with random distances"), fill=c("#756bb1", "#7fcdbb"), cex=0.8)

png(file="./saving_plot2.png",
    width=600, height=600)
layout(mat = matrix(c(1,2,1,2),2,1, byrow=TRUE),  height = c(1,8))
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(data_p$V1 , horizontal=TRUE , ylim=c(100, 500), xaxt="n" , col="#7fcdbb" , frame=F)
par(mar=c(0, 3.1, 1.1, 2.1))
boxplot(data_s$V1 , horizontal=TRUE , ylim=c(100, 500), xaxt="n" , col="#756bb1" , frame=F, add =T)
par(mar=c(4, 3.1, 1.1, 2.1)) 
hist(data_p$V1 , breaks=40 , col="#7fcdbb" , border=F , main="" , xlab="Distance (m)", xlim=c(100,500))
par(mar=c(4, 3.1, 1.1, 2.1)) 
hist(data_s$V1 , breaks=40 , col= "#756bb1", border=F , main="" , xlab="Distance (m)", xlim=c(100,500), add = T)
abline(v = mean(real_distances$distancias_diag), col = "black", lwd = 3, lty = "dashed", add = T)
dev.off()


