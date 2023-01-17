#' Script for calculating the temporal autocorrelation within each nestbox.
#' I used these resources
#' 
#' https://aosmith.rbind.io/2018/06/27/uneven-grouped-autocorrelation/
#' https://sakai.unc.edu/access/content/group/2842013b-58f5-4453-aa8d-3e01bacbfc3d/public/Ecol562_Spring2012/docs/lectures/lecture17.htm#testing
#' 

#--Library
library(tidyverse)
library(lattice)
library(stringr)
library(lme4)
library(ggpubr)
library(DHARMa)
library(lmerTest)


#--Data
nests <- read.csv("./data/GPS.csv", sep = ";", head = T, dec = ".")
rs <- read.csv("./data/exito_rep.csv", sep = ",", head = T, dec = ".")

# Customizing
rs <- rs[,-1]
rs <- rs %>% filter(! YEAR_NEST == "#N/D") 
rs <- rs %>% rename(first_broods = X1, second_broods = X2)
rs <- rs %>% mutate(total_vols = first_broods + second_broods) 
rs_join <- rs %>% dplyr::select(1, 4,6, 7) 
rs_join <- rs_join %>% distinct(YEAR_NEST, .keep_all = T)

# Nest data
nests_data <- array(0, dim = c(nrow(nests), ncol = 10)) 
nests_data <- as.data.frame(nests_data) 
nests_data$CAJA <- nests$CAJA 
nests_data <- nests_data %>% relocate(CAJA, .before = everything()) 
colnames(nests_data)[2:11] <- c(2012:2021) 

# Changing the dataframe to longformat to work with repeated measures.
nest_data_longer <- nests_data %>% pivot_longer(!CAJA, names_to = "year", values_to = "reproductive_success")
nest_data_longer <- nest_data_longer %>% mutate(year_nest = paste(year, CAJA, sep = "_"))
nest_data_longer <- nest_data_longer %>% mutate(bred = ifelse(nest_data_longer$year_nest %in% rs $YEAR_NEST, 1, 0))

# Together
nestdata_joined_rs <- nest_data_longer %>% left_join(rs_join, by = c("year_nest" = "YEAR_NEST"))
nestdata_joined_rs$total_vols[nestdata_joined_rs$bred==0] <- 0 
# We input an NA in those nest-boxes that have been subjected to experimentation
nestdata_joined_rs$total_vols[nestdata_joined_rs$EXPERIMENTACION_IRAIDA==1] <- NA 


# Calculation of the temporal autocorrelation
# following: #' https://aosmith.rbind.io/2018/06/27/uneven-grouped-autocorrelation/
# We calculate the temporal autocorrelation up to 5 time lags

temporal_autocor <- nestdata_joined_rs %>% dplyr::select(1, 2, 8)
temporal_autocor <- temporal_autocor %>% rename(unit = CAJA, time = year, rs = total_vols) # renaming
temporal_autocor <- temporal_autocor %>% filter(time %in% c(2012:2017))

temporal_autocor <- temporal_autocor %>% mutate(time2 = case_when(time == 2012 ~ 1,   
                                                                  time == 2013 ~ 2,
                                                                  time == 2014 ~ 3,
                                                                  time == 2015 ~ 4,
                                                                  time == 2016 ~ 5,
                                                                  time == 2017 ~ 6))

temporal_autocor <- temporal_autocor %>% relocate(time2, .after = time) 
temporal_autocor <- temporal_autocor %>% dplyr::select(1, 3, 4) 

data <- temporal_autocor 
data$unit <- as.factor(data$unit) 

# Checking
# Number of nest boxes
data %>%
  count(unit)           

# Number of observation per nest 
data %>%
  count(unit) %>%       
  filter(n == max(n))

# We arrange the dataset to do a correct calculation of the temporal autocorrelation.
data = data %>%           
  arrange(unit, time2)

# We add 'fake observations' to be able to separate each nest box.
dat_expand = data %>%     
  group_by(unit) %>%
  complete(time2 = 1:12)


filter(dat_expand, unit == 1) # Example

(nall = map_df(1:5, 
               ~data %>% 
                 drop_na(rs) %>%
                 group_by(unit) %>%
                 arrange(unit, time2) %>%
                 summarise(lag = list( diff(time2, lag = .x ) ) )
) %>%
    unnest(lag) %>%
    group_by(lag) %>% 
    summarise(n = n()))


# Calculation of autocorrelation
acf_g <- acf(dat_expand$rs, lag.max = 5, na.action = na.pass, ci = 0, ylim=c(-0.4, 0.4), xlab = "Year Lags", ylab = "Autocorrelation", main = "Reproductive Success' Temporal Autocorrelation ")
ci_upper <- c(qnorm(1-.025)/sqrt(nall$n))
ci_lower <- c(qnorm(1-.025)/sqrt(nall$n))

# Putting it in a dataframe
acf <- acf_g$acf
acf <- acf[,,1] # 3rd dimension of the array
acf # this


acd <- data.frame(acf = acf[2:6] , ci_upper = ci_upper, ci_lower = ci_lower, lag = c(1:5))


# Plot
acf_plot <- acd %>% 
  #as_tibble() %>% mutate(lags = 1:n()) %>% 
  ggplot(aes(x=lag, y = acf)) + scale_x_continuous(breaks=seq(0,6,1)) +
  #geom_point(aes(y=ci_upper,  col='blue')) +
  #geom_point(aes(y = ci_lower , col = "red"))+
  geom_line(aes(y=ci_upper), linetype = "dashed", col = "purple", lwd = 1)+
  geom_line(aes(y=ci_lower), linetype = "dashed", col = "purple", lwd = 1)+
  geom_hline(yintercept = 0)+
  labs(y="Autocorrelations", x="Lag") +
  scale_y_continuous("Autocorrelation", limits = c(-0.4,0.4), breaks = seq(-0.4,0.4, by = 0.10
  )) +
  geom_segment(aes(xend=lag, yend=0), lwd = 0.8) +geom_point(aes(y = acf), size = 3) + 
  theme(axis.title.x= element_text(face="plain", size=16,  margin = margin(t=20, r=0, l=0, b=0)),
        axis.title.y= element_text(face="plain", size=16, margin = margin(t=0, r=20, l=0, b=0)),
        axis.text.x = element_text(face="plain", size =14, vjust = 1, hjust = 1),
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
        panel.spacing.y = unit(2, "lines"))

