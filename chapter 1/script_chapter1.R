# Chapter 1 - Sexual differences in phenotypical predictors of floating status:
# body condition influences male but not female reproductive status in a wild passerine      


# Adult female analysis. 
# Dataset can be found on the following link: https://figshare.com/articles/dataset/Data_rar/19705033

##--Set working directory
setwd("C:/Users/iraid/Desktop/females/females")
##--Load data
data <- read.csv("adult_females.csv", header= T, sep = ";", dec = ".")

##--Checking variable types
str(data)
data$status <- as.factor(as.character(data$status))
data$cert <- as.numeric(data$cert)
data$capture_year <- as.factor(as.character(data$capture_year))
data$ring <- as.factor(as.character(data$ring))
data$observer <- as.factor(data$observer)


##--Considering the variability among observers at measuring
tarsus <- lm(tarsus ~ observer, data = data) # for tarsus length
beak<- lm(beak ~ observer, data = data)      # for beak length
wing <- lm(wing ~ observer, data = data)     # for wing length



##--Joining  the residuals of the last models to the dataframe

# Tarsus length 
ID_birds <- data %>% 
  select(ring, tarsus) %>% # select columns
  na.omit() # remove NAs
  nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring, "tarso_cor" =  residuals(tarsus)) 
data <- left_join(data, residuals_data_a) # left join


# Wing length 
ID_birds <- data %>% 
  select(ring, wing) %>% # select columns
  na.omit() # remove NAs
  nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring,"ala_cor" =  residuals(wing)) 
data <- left_join(data, residuals_data_a) # left join


# Beak length
ID_birds <- data %>% 
  select(ring, beak) %>% # select columns
  na.omit() # remove NAs
  nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring, "beak_cor" =  residuals(beak)) 
data <- left_join(data, residuals_data_a) # left join



##--Calculation of a size index by means of PCA
data_pc <- data # duplicate

# PCA 
pca1 <- prcomp(data_pc[,16:18], center = TRUE, scale = TRUE) # selection of the measurements corrected by observer.

# To see the eigenvalues, and the % variance explained
eig <- (pca1$sdev)^2
variance <-eig*100/sum(eig)
cumvar <- cumsum(variance)
tabla.eigenv <- data.frame(eig = eig, variance = variance, cumvar = cumvar)
print(tabla.eigenv)

# Plot for eigenvalues
barplot(tabla.eigenv[, 2], names.arg=1:nrow(tabla.eigenv),
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variance explained",
        col = "red")

#To see the variable loadings
unclass(pca1$rotation)[,1:3]
data_pc$size <- pca1$x[,1]*(-1) # adding to the dataframe as a measure of 'size'.



## Model to calculate a body condition index including potential influencing variables such as the timing of the capture (hour of day),
# and capture date. Residuals of this model were used as a body condition index. 

# First model 
condition <- lmer(weight ~ scale(size) + scale(time_capture) + scale(capture_date)  + (1|capture_year), data = data_pc) 
summary(condition)

# We remove time of capture from the model. 
condition <- lmer(weight ~ scale(size)  + scale(capture_date)  + (1|capture_year), data = data_pc) 
summary(condition)

# Checking the normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(condition), ylab="Residuals")
qqline(residuals(condition))
shapiro.test(residuals(condition))
hist(residuals(condition))  

# Homocedasticity 
fitted.values=fitted(condition)
residual.values=residuals(condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  #

# Overall checking.
check_model(condition)

# Join body condition (residuals) to the dataframe
ID_birds <- data_pc %>% 
  select(ring, weight, tarsus) %>% 
  na.omit() # remove NAs
  nrow(ID_birds) # # checking

residuals_data_c <- data.frame("ring" = ID_birds$ring,"adult_condition" =  residuals(condition))
dim(residuals_data_c) # checking
data_pc <- left_join(data_pc, residuals_data_c) # left join


# Floater is coded as the reference level.
data_pc$status <- relevel(data_pc$status, ref = "F")


##--Multivariate model for  adult females

adult_females <- glmer(status ~ scale(adult_condition) +scale(size) + scale(feather_length) + scale(spottiness) + (1|capture_year), family = "binomial", data = data_pc)
summary(adult_females)


## Use of DHARMa package to check the residuals of the GLMM model. 
simulationOutput <- simulateResiduals(fittedModel = adult_females) 
plot(simulationOutput)  # Fine. 


# Using drop1 to get the p-values. 
drop1(adult_females, test = "Chisq")


# Calculation of the odd ratios (OR)
se <- sqrt(diag(vcov(adult_females)))

# Table of estimates in the logit scale with 95% CI
(tab <- cbind(Est = fixef(adult_females), LL = fixef(adult_females) - 1.96 * se, UL = fixef(adult_females) + 1.96 * se))

# 
exp(tab)   # if OR = 1, no association is detected; OR > 1, positive association; OR < 1, negative 
# association; if the 95% CI includes the 1 or values near one = no significant. 


rm(list=ls())


# Nestling female analysis. 
# Dataset can be found on the following link: https://figshare.com/articles/dataset/Data_rar/19705033

##--Set working directory
##--Load data-
data <- read.csv("nestling_females.csv", header= T, sep = ";", dec = ".")

##--Checking variable types
str(data)
data$status <- as.factor(as.character(data$status))
data$cert <- as.numeric(data$cert)
data$cohort <- as.factor(as.character(data$cohort))
data$ring <- as.factor(as.character(data$ring))
data$observer2 <- as.factor(data$observer2)

##--Considering the variability among observers at measuring
tarsus <- lm(tarsus ~ observer2, data = data) # for tarsus length
check_model(tarsus)
summary(tarsus)
##--Joining  the residuals of the last models to the dataframe

# Tarsus length 
ID_birds <- data %>% 
  select(ring, tarsus) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring, "tarso_cor" =  residuals(tarsus)) 
data <- left_join(data, residuals_data_a) # left join


## Model to calculate a body condition index.

# First model 
condition <- lmer(weight ~ scale(tarso_cor) + (1|cohort), data = data) 
summary(condition)


# Checking the normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(condition), ylab="Residuals")
qqline(residuals(condition))
shapiro.test(residuals(condition))
hist(residuals(condition))  

# Homocedasticity 
fitted.values=fitted(condition)
residual.values=residuals(condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  #

# Overall checking.
check_model(condition)


# Join body condition (residuals) to the dataframe
ID_birds <- data %>% 
  select(ring, weight, tarsus) %>% 
  na.omit() # remove NAs
  nrow(ID_birds) # # checking

residuals_data_c <- data.frame("ring" = ID_birds$ring,"nestling_condition" =  residuals(condition))
dim(residuals_data_c) # checking
data <- left_join(data, residuals_data_c) # left join


# Floater is coded as the reference level.
data$status <- relevel(data$status, ref = "F")


##--Multivariate model for  adult females

nestling_females <- glmer(status ~ scale(run_day_eclosion) + scale(nestling_condition) + scale(brood_size)+ (1|cohort), family = "binomial", data = data)
summary(nestling_females)


## Use of DHARMa package to check the residuals of the GLMM model. 
simulationOutput <- simulateResiduals(fittedModel = nestling_females) 
plot(simulationOutput)  # Fine. 


# Using drop1 to get the p-values. 
drop1(nestling_females, test = "Chisq")


# Calculation of the odd ratios (OR)
se <- sqrt(diag(vcov(nestling_females)))

# Table of estimates in the logit scale with 95% CI
(tab <- cbind(Est = fixef(nestling_females), LL = fixef(nestling_females) - 1.96 * se, UL = fixef(nestling_females) + 1.96 * se))

# 
exp(tab)   # if OR = 1, no association is detected; OR > 1, positive association; OR < 1, negative 
# association; if the 95% CI includes the 1 or values near one = no significant. 


rm(list = ls())

# Adult Male analysis. 
# Dataset can be found on the following link: https://figshare.com/articles/dataset/Data_rar/19705033
##--Set working directory
setwd("C:/Users/iraid/Desktop/males/males")
##--Load data
data <- read.csv("adult_males.csv", header= T, sep = ";", dec = ".")
data <- data[-(281:285),] # delete 5 empty rows that are added in the csv

##--Checking variable types
str(data)
data$status <- as.factor(as.character(data$status))
data$cert <- as.numeric(data$cert)
data$capture_year <- as.factor(as.character(data$capture_year))
data$ring <- as.factor(as.character(data$ring))
data$observer <- as.factor(data$observer)


##--Considering the variability among observers at measuring
tarsus <- lm(tarsus ~ observer, data = data) # for tarsus length
beak<- lm(beak ~ observer, data = data)      # for beak length
wing <- lm(wing ~ observer, data = data)     # for wing length



##--Joining  the residuals of the last models to the dataframe

# Tarsus length 
ID_birds <- data %>% 
  select(ring, tarsus) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring, "tarso_cor" =  residuals(tarsus)) 
data <- left_join(data, residuals_data_a) # left join


# Wing length 
ID_birds <- data %>% 
  select(ring, wing) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring,"ala_cor" =  residuals(wing)) 
data <- left_join(data, residuals_data_a) # left join


# Beak length
ID_birds <- data %>% 
  select(ring, beak) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring, "beak_cor" =  residuals(beak)) 
data <- left_join(data, residuals_data_a) # left join



##--Calculation of a size index by means of PCA
data_pc <- data # duplicate

# PCA 
pca1 <- prcomp(data_pc[,15:17], center = TRUE, scale = TRUE) # selection of the measurements corrected by observer.

# To see the eigenvalues, and the % variance explained
eig <- (pca1$sdev)^2
variance <-eig*100/sum(eig)
cumvar <- cumsum(variance)
tabla.eigenv <- data.frame(eig = eig, variance = variance, cumvar = cumvar)
print(tabla.eigenv)

# Plot for eigenvalues
barplot(tabla.eigenv[, 2], names.arg=1:nrow(tabla.eigenv),
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variance explained",
        col = "red")

#To see the variable loadings
unclass(pca1$rotation)[,1:3]
data_pc$size <- pca1$x[,1] # adding to the dataframe as a measure of 'size'.



## Model to calculate a body condition index including potential influencing variables such as the timing of the capture (hour of day),
# and capture date. Residuals of this model were used as a body condition index. 

# First model 
condition <- lmer(weight ~ scale(size) + scale(time_capture) + scale(date_january) + (1|capture_year), data = data_pc) 
summary(condition)

# We remove capture date from the model. 
condition <- lmer(weight ~ scale(size)  + scale(time_capture)  + (1|capture_year), data = data_pc) 
summary(condition)

# Checking the normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(condition), ylab="Residuals")
qqline(residuals(condition))
shapiro.test(residuals(condition))
hist(residuals(condition))  

# Homocedasticity 
fitted.values=fitted(condition)
residual.values=residuals(condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  #

# Overall checking.
check_model(condition)

# Join body condition (residuals) to the dataframe
ID_birds <- data_pc %>% 
  select(ring, weight, tarsus) %>% 
  na.omit() # remove NAs
nrow(ID_birds) # # checking

residuals_data_c <- data.frame("ring" = ID_birds$ring,"adult_condition" =  residuals(condition))
dim(residuals_data_c) # checking
data_pc <- left_join(data_pc, residuals_data_c) # left join


# Floater is coded as the reference level.
data_pc$status <- relevel(data_pc$status, ref = "F")


##--Multivariate model for  adult females

adult_males <- glmer(status ~ scale(adult_condition) +scale(size) + scale(feather_length) + (1|capture_year), family = "binomial", data = data_pc)
summary(adult_males)


## Use of DHARMa package to check the residuals of the GLMM model. 
simulationOutput <- simulateResiduals(fittedModel = adult_males) 
plot(simulationOutput)  # Fine. 


# Using drop1 to get the p-values. 
drop1(adult_males, test = "Chisq")


# Calculation of the odd ratios (OR)
se <- sqrt(diag(vcov(adult_males)))

# Table of estimates in the logit scale with 95% CI
(tab <- cbind(Est = fixef(adult_males), LL = fixef(adult_males) - 1.96 * se, UL = fixef(adult_females) + 1.96 * se))

# 
exp(tab)   # if OR = 1, no association is detected; OR > 1, positive association; OR < 1, negative 
# association; if the 95% CI includes the 1 or values near one = no significant. 


rm(list=ls())




# Nestling male analysis. 
# Dataset can be found on the following link: https://figshare.com/articles/dataset/Data_rar/19705033

##--Set working directory
##--Load data
data <- read.csv("nestling_males.csv", header= T, sep = ";", dec = ".")

##--Checking variable types
str(data)
data$status <- as.factor(as.character(data$status))
data$cert <- as.numeric(data$cert)
data$cohort <- as.factor(as.character(data$cohort))
data$ring <- as.factor(as.character(data$ring))
data$observer <- as.factor(data$observer)

##--Considering the variability among observers at measuring
tarsus <- lm(tarsus ~ observer, data = data) # for tarsus length
check_model(tarsus)
summary(tarsus)
##--Joining  the residuals of the last models to the dataframe

# Tarsus length 
ID_birds <- data %>% 
  select(ring, tarsus) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring, "tarso_cor" =  residuals(tarsus)) 
data <- left_join(data, residuals_data_a) # left join


## Model to calculate a body condition index.

# First model 
condition <- lmer(weight ~ scale(tarso_cor) + (1|cohort), data = data) 
summary(condition)


# Checking the normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(condition), ylab="Residuals")
qqline(residuals(condition))
shapiro.test(residuals(condition))
hist(residuals(condition))  

# Homocedasticity 
fitted.values=fitted(condition)
residual.values=residuals(condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  #

# Overall checking.
check_model(condition)


# Join body condition (residuals) to the dataframe
ID_birds <- data %>% 
  select(ring, weight, tarsus) %>% 
  na.omit() # remove NAs
nrow(ID_birds) # # checking

residuals_data_c <- data.frame("ring" = ID_birds$ring,"nestling_condition" =  residuals(condition))
dim(residuals_data_c) # checking
data <- left_join(data, residuals_data_c) # left join


# Floater is coded as the reference level.
data$status <- relevel(data$status, ref = "F")


##--Multivariate model for nestling males.

nestling_males <- glmer(status ~ scale(run_day_eclosion) + scale(nestling_condition) + scale(brood_size)+ (1|cohort), family = "binomial", data = data)
summary(nestling_males)


## Use of DHARMa package to check the residuals of the GLMM model. 
simulationOutput <- simulateResiduals(fittedModel = nestling_males) 
plot(simulationOutput)  # Fine. 

# Using drop1 to get the p-values. 
drop1(nestling_males, test = "Chisq")

# Calculation of the odd ratios (OR)
se <- sqrt(diag(vcov(nestling_males)))

# Table of estimates in the logit scale with 95% CI
(tab <- cbind(Est = fixef(nestling_males), LL = fixef(nestling_males) - 1.96 * se, UL = fixef(nestling_males) + 1.96 * se))

# 
exp(tab)   # if OR = 1, no association is detected; OR > 1, positive association; OR < 1, negative 
# association; if the 95% CI includes the 1 or values near one = no significant.  


rm(list = ls())


## Path analysis - Females

##--Set working directory
setwd("C:/Users/iraid/Desktop/females/females")
##--Load data-
data <- read.csv("path_analysis_females.csv", header= T, sep = ";", dec = ".")


#--Checking point
str(data)
data$status <- as.factor(as.character(data$status))
data$status_cat <- as.factor(as.character(data$status_cat))
data$cert <- as.numeric(data$cert)
data$cohort <- as.factor(as.character(data$cohort))
data$capture_year <- as.factor(as.character(data$capture_year))
data$ring <- as.factor(as.character(data$ring))
data$observer_n <- as.factor(as.character(data$observer_n))
data$observer_a <- as.factor(as.character(data$observer_a))



##--Considering the variability among observers at measuring
tarsus_a <- lm(tarsus_a ~ observer_a, data = data)
beak_a <- lm(beak ~ observer_a, data = data)
wing_a <- lm(wing_a ~ observer_a, data = data)


##--Joining  the residuals of the last models to the dataframe
ID_birds <- data %>% 
  select(ring, tarsus_a) %>% # select columns
  na.omit() # remove NAs
  nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring,
                               "tarsus_cor" =  residuals(tarsus_a)) #adding our size var. corrected.
data <- left_join(data, residuals_data_a) # join


ID_birds <- data %>% 
  select(ring, wing_a) %>% # select columns
  na.omit() # remove NAs
  nrow(ID_birds) 

residuals_data_a <- data.frame("ring" = ID_birds$ring,
                               "wing_cor" =  residuals(wing_a)) #adding our size var. corrected.
data <- left_join(data, residuals_data_a) # join


ID_birds <- data %>% 
  select(ring, beak) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?

residuals_data_a <- data.frame("ring" = ID_birds$ring,
                               "beak_cor" =  residuals(beak_a)) #adding our size var. corrected.
data <- left_join(data, residuals_data_a) # join



##--Calculation of a size index by means of PCA
pca1 <- prcomp(data[,22:24], center = TRUE, scale = TRUE)

# To see the Eigenvalues, and the % variance explained
eig <- (pca1$sdev)^2
variance <-eig*100/sum(eig)
cumvar <- cumsum(variance)
tabla.eigenv <- data.frame(eig = eig, variance = variance, cumvar = cumvar)
print(tabla.eigenv)

#Scree plot for eigenvalues
barplot(tabla.eigenv[, 2], names.arg=1:nrow(tabla.eigenv),
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variance explained",
        col ="red")

#to see the variable loadings
unclass(pca1$rotation)[,1:3]
data$size <- pca1$x[,1]*(-1)


## Model to calculate a body condition index including potential influencing variables such as the timing of the capture (hour of day),
# and capture date. Residuals of this model were used as a body condition index. 

condition <- lmer(weight_a ~ scale(size) + scale(capture_date) + (1|capture_year), data = data) 
summary(condition) 


# Normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(condition), ylab="Residuals")
qqline(residuals(condition))
shapiro.test(residuals(condition))
hist(residuals(condition))  

# Homocedasticity 
fitted.values=fitted(condition)
residual.values=residuals(condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  #


# Join body condition (residuals) to the dataframe
ID_birds <- data %>% 
  select(ring, weight_a, tarsus_a) %>% 
  na.omit() # remove NAs
  nrow(ID_birds) #

residuals_data_b <- data.frame("ring" = ID_birds$ring,
                               "adult_condition" =  residuals(condition))

dim(residuals_data_b) #check
data <- left_join(data, residuals_data_b) # left join

##--Considering the variability among observers at measuring in nestlings
tarsus_n <- lm(tarsus_n ~ observer_n, data = data)

# to check normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(tarsus_n), ylab="Residuals")
qqline(residuals(tarsus_n))
shapiro.test(residuals(tarsus_n))
hist(residuals(tarsus_n))  # nop

# to check homocedasticity
fitted.values=fitted(tarsus_n)
residual.values=residuals(tarsus_n)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  # ok

# Joining to the dataframe

ID_birds <- data %>% 
  select(ring, tarsus_n) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?

residuals_data_n <- data.frame("ring" = ID_birds$ring,
                               "tarsus_cor_n" =  residuals(tarsus_n)) #adding our tarsus var. corrected.
head(residuals_data_n, 20)

data <- left_join(data, residuals_data_n) # join

## Model to calculate a body condition index 
nestling_condition <- lm(weight_n ~ scale(tarsus_cor_n) + (1|cohort), data = data) 

# to check normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(nestling_condition), ylab="Residuals")
qqline(residuals(nestling_condition))
shapiro.test(residuals(nestling_condition))
hist(residuals(nestling_condition))  # ew

# to check homocedasticity
fitted.values=fitted(nestling_condition)
residual.values=residuals(nestling_condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  


# Join body condition (residuals) to the dataframe
ID_birds <- data %>% 
  select(ring, weight_n, tarsus_n) %>% # select columns
  na.omit() # remove NAs
  nrow(ID_birds) 

residuals_data_n <- data.frame("ring" = ID_birds$ring,
                               "nestling_condition" =  residuals(nestling_condition))
head(residuals_data_n, 20)

data <- left_join(data, residuals_data_n) # left join

## Path analysis

str(data)

write.csv(data, "path_revision_females.csv")
data <- read.csv("path_revision_females.csv") # this is the only way I made psem() worked with categorical response variables.

## Path analysis for females: 

path_females <- psem(
  M1 <- glm(status_cat ~ nestling_condition + adult_condition + feather_length  + spottiness + run_day_eclosion + brood_size , family = binomial, data = data),
  M2 <- lm(adult_condition ~ nestling_condition + run_day_eclosion + brood_size , data = data) ,
  M3 <- lm(feather_length ~ adult_condition , data = data),
  M5 <- lm(spottiness ~run_day_eclosion , data = data),
  M4 <- lm(nestling_condition ~ run_day_eclosion + brood_size , data = data),data = data) 

summary(path_females, standardize = "scale")


## Path analysis - Males

##--Set working directory
setwd("C:/Users/iraid/Desktop/males/males")
setwd("C:/Users/iraid/Desktop/comparacion")

data <- read.csv("path_analysis_males.csv", header= T, sep = ";", dec = ".")
data2 <- read.csv("data_machos_2907.csv", header= T, sep = ";", dec = ".")

##--Load data
data <- read.csv("path_analysis_males.csv", header= T, sep = ";", dec = ".")

#--Checking point
str(data)
data$status <- as.factor(as.character(data$status))
data$status_cat <- as.factor(as.character(data$status_cat))
data$cert <- as.numeric(data$cert)
data$cohort <- as.factor(as.character(data$cohort))
data$capture_year <- as.factor(as.character(data$capture_year))
data$ring <- as.factor(as.character(data$ring))
data$observer_n <- as.factor(as.character(data$observer_n))
data$observer_a2 <- as.factor(as.character(data$observer_a2))
str(data)

##--Considering the variability among observers at measuring
tarsus_ad <- lm(tarsus_a ~ observer_a2, data = data)
beak_ad <- lm(beak ~ observer_a2, data = data)
wing_ad <- lm(wing_a ~ observer_a2, data = data)

summary(wing_ad)


##--Joining  the residuals of the last models to the dataframe
ID_birds <- data %>% 
  select(ring, tarsus_a) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # checking

residuals_data_a <- data.frame("ring" = ID_birds$ring,
                               "tarsus_cor" =  residuals(tarsus_ad)) #adding our size var. corrected.
data <- left_join(data, residuals_data_a) # join


ID_birds <- data %>% 
  select(ring, wing_a) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) 

residuals_data_a <- data.frame("ring" = ID_birds$ring,
                               "wing_cor" =  residuals(wing_ad)) #adding our size var. corrected.
data <- left_join(data, residuals_data_a) # join


ID_birds <- data %>% 
  select(ring, beak) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?

residuals_data_a <- data.frame("ring" = ID_birds$ring,
                               "beak_cor" =  residuals(beak_ad)) #adding our size var. corrected.
data <- left_join(data, residuals_data_a) # join



##--Calculation of a size index by means of PCA
pca1 <- prcomp(data[,23:25], center = TRUE, scale = TRUE)

# To see the Eigenvalues, and the % variance explained
eig <- (pca1$sdev)^2
variance <-eig*100/sum(eig)
cumvar <- cumsum(variance)
tabla.eigenv <- data.frame(eig = eig, variance = variance, cumvar = cumvar)
print(tabla.eigenv)

#Scree plot for eigenvalues
barplot(tabla.eigenv[, 2], names.arg=1:nrow(tabla.eigenv),
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variance explained",
        col ="red")

#to see the variable loadings
unclass(pca1$rotation)[,1:3]
data$size <- pca1$x[,1]


## Model to calculate a body condition index including potential influencing variables such as the timing of the capture (hour of day),
# and capture date. Residuals of this model were used as a body condition index. 

condition <- lmer(weight_a ~ scale(size) + (1|capture_year), data = data) 
summary(condition) 


# Normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(condition), ylab="Residuals")
qqline(residuals(condition))
shapiro.test(residuals(condition))
hist(residuals(condition))  

# Homocedasticity 
fitted.values=fitted(condition)
residual.values=residuals(condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  #


# Join body condition (residuals) to the dataframe
ID_birds <- data %>% 
  select(ring, weight_a) %>% 
  na.omit() # remove NAs
nrow(ID_birds) #

residuals_data_b <- data.frame("ring" = ID_birds$ring,
                               "adult_condition" =  residuals(condition))

dim(residuals_data_b) #check
data <- left_join(data, residuals_data_b) # left join

##--Considering the variability among observers at measuring in nestlings
tarsus_nt <- lm(tarsus_n ~ observer_n, data = data)

# to check normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(tarsus_nt), ylab="Residuals")
qqline(residuals(tarsus_nt))
shapiro.test(residuals(tarsus_nt))
hist(residuals(tarsus_nt))  # nop

# to check homocedasticity
fitted.values=fitted(tarsus_nt)
residual.values=residuals(tarsus_nt)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  # ok

# Joining to the dataframe

ID_birds <- data %>% 
  select(ring, tarsus_n) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?

residuals_data_n <- data.frame("ring" = ID_birds$ring,
                               "tarsus_cor_n" =  residuals(tarsus_nt)) #adding our tarsus var. corrected.
head(residuals_data_n, 20)

data <- left_join(data, residuals_data_n) # join

## Model to calculate a body condition index 
nestling_condition <- lmer(weight_n ~ scale(tarsus_cor_n) + (1|cohort), data = data) 

# To check normality of residuals
par(mfrow=c(1,2))
qqnorm(residuals(nestling_condition), ylab="Residuals")
qqline(residuals(nestling_condition))
shapiro.test(residuals(nestling_condition))
hist(residuals(nestling_condition))  # ew

# To check homocedasticity
fitted.values=fitted(nestling_condition)
residual.values=residuals(nestling_condition)
plot(residual.values~fitted.values)
abline(mean(residual.values),0)  


# Join body condition (residuals) to the dataframe
ID_birds <- data %>% 
  select(ring, weight_n) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) 

residuals_data_n <- data.frame("ring" = ID_birds$ring,
                               "nestling_condition" =  residuals(nestling_condition))
head(residuals_data_n, 20)

data <- left_join(data, residuals_data_n) # left join

## Path analysis

str(data)

write.csv(data, "path_revision_males.csv")

data <- read.csv("path_revision_males.csv") # this is the only way I made psem() worked with categorical response variables.

## Path analysis for males 

path_machos_new <- psem(
  H1 <- glm(status_cat ~ nestling_condition + adult_condition + feather_length   +run_day_eclosion + brood_size, family = "binomial", data = data),
  H2 <- lm(adult_condition ~ nestling_condition + run_day_eclosion + brood_size, data = data),
  H3 <- lm(feather_length ~ adult_condition, data = data),
  H5 <- lm(nestling_condition ~ run_day_eclosion + brood_size, data = data), data = data)

summary(path_machos_new, standardize = "scale")
