

# R script containing the code used to analyze the relationship between variables of prospected nestboxes and the nest that 
# was chosen on the following year. 


##--Libraries 
library(tidyverse)
library(lattice)
library(stringr)
library(lme4)
library(ggpubr)
library(DHARMa)
library(lmerTest)
library(glmmTMB)
library(performance)
library(GLMMadaptive)

options(contrasts=c(factor="contr.sum", ordered="contr.poly"))


#--Data 
data <- read.csv("chapter 2/data.csv", sep = ";", head = T )
length(unique(data$ring)) # 139 unique males. 
str(data)


#--Checking variable type.
data$chosen_dummy <- as.factor(data$chosen_dummy)
data$ring <- as.factor(data$ring)
data$nest_dummy <- as.factor(data$nest_dummy)
data$visited <- as.factor(data$visited)
data$chosen <- as.factor(data$chosen)
data$fuera_del_circulo <- as.factor(data$fuera_del_circulo)
data$año_recaptura <- as.factor(data$año_recaptura)
data$YEAR_NEST <- as.factor(data$YEAR_NEST)
data$excluir <- as.factor(data$excluir)
data$EXPERIMENTACION_IRAIDA <- as.factor(data$EXPERIMENTACION_IRAIDA)
data$YEAR_NEST_rep <- as.factor(data$YEAR_NEST_rep)
data$nest_dummy_rep <- as.factor(data$nest_dummy_rep)
data$rep_prospected <- as.factor(data$rep_prospected)
data$rep_inside_circle <- as.factor(data$rep_inside_circle)
data$YEAR_NEST_rep_previous <- as.factor(data$YEAR_NEST_rep_previous)
data$EXPERIMENTACION_IRAIDA_rep <- as.factor(data$EXPERIMENTACION_IRAIDA_rep)
data$tipo_macho <- as.factor(data$tipo_macho)
data$macho_a_buscarT <- as.factor(data$macho_a_buscarT)
data$macho_a_buscarT1 <- as.factor(data$macho_a_buscarT1)
data$observT <- as.factor(data$observT)
data$observT1 <- as.factor(data$observT1)
data$destino_CERT <- as.factor(data$destino_CERT)
data$next_YEAR_NEST_CERT <- as.factor(data$next_YEAR_NEST_CERT)
data$next_NEST_CERT <- as.factor(data$next_NEST_CERT)
data$same_nest_CERT <- as.factor(data$same_nest_CERT)
data$destino_sin_CERT <- as.factor(data$destino_sin_CERT)
data$next_YEAR_NEST_sin_CERT <- as.factor(data$next_YEAR_NEST_sin_CERT)
data$next_NEST_sin_CERT <- as.factor(data$next_NEST_sin_CERT)
data$same_nest_sin_CERT <- as.factor(data$same_nest_sin_CERT)
data$observT_focal <- as.factor(data$observT_focal)
data$observT1_focal <- as.factor(data$observT1_focal)
data$codigo_destin_CERT <- as.factor(data$codigo_destin_CERT)
data$dummy_destin <- as.factor(data$dummy_destin)
data$mean_feeding_rate <- as.numeric(data$mean_feeding_rate)
str(data)

# Barplot showing the percentage of nest-boxes in which floaters bred that were
# visited (green bar) or were included inside the projected circle (purple bar) the previous year.

elegidas <- data %>% filter(chosen == 1)
descriptive_visitadas <- elegidas %>% summarise(count = sum(elegidas$visited ==1)/139) # 13%
descriptive_dentro_circulo <- elegidas %>% summarise(count = sum(elegidas$fuera_del_circulo ==0)/139) # 39%

descriptive_dataframe <- data.frame(porcentaje = c("13%", "39%"), tipo = c("Visited", "Inside prospected circle"))
descriptive_dataframe$tipo <- factor(descriptive_dataframe$tipo, levels = c("Visited", "Inside prospected circle"))

percentage <- ggplot(data=descriptive_dataframe, aes(x=tipo, y=porcentaje, fill=tipo)) +
  geom_bar(stat="identity", position=position_dodge(),  width = 0.5)+
  geom_text(aes(label=porcentaje),  vjust=-0.3, size=3, color="black",
            position = position_dodge(0.9))+
  scale_fill_brewer(palette="Accent") + 
  theme(axis.title.x = element_blank(),
        axis.title.y= element_blank(),
        axis.text.x = element_text(face="plain", size =8),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 1, face = "bold", size = 16),
        legend.position = "none", 
        legend.title = element_text(colour="black", size=14, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.3, linetype = "solid"),
        panel.spacing.y = unit(2, "lines"))

plot(percentage)


##-- Logistic regressions
# Response variable: chosen vs. not chosen (1 vs. 0). 
# Explanatory variable: total number of fledged nestlings in that nestbox during the prospecting year . 
# Not chosen used as the reference level.
data$chosen_dummy <- relevel(data$chosen_dummy, ref="N")
data <- data %>% filter(excluir == 0)

# Simplifying dataset
volantones <- data %>% dplyr::select(2, 3, 5, 6, 9, 13, 14) 

# Model
mod1 <- glmer(chosen ~ vols_año_prospeccion + (1|ring/año_recaptura) , na.action = na.exclude, data = volantones, family = binomial)
summary(mod1) # Singular fit!

## Residuals
simulationOutput <- simulateResiduals(fittedModel = mod1) 
plot(simulationOutput)  # Fine. 
testDispersion(mod1)
plot(simulationOutput, data$vols_año_prospeccion)

## P-value
drop1(mod1, test = "Chisq")

# From logit scale to odds ratio
se <- sqrt(diag(vcov(mod1)))
# Table of estimates with 95% CI
(tab <- cbind(Est = fixef(mod1), LL = fixef(mod1) - 1.96 * se, UL = fixef(mod1) + 1.96 * se))
exp(tab)   

# The same as the latter but only considering the fledged nestlings of the second clutches.

# Model
mod1_bis <- glmer(chosen ~ vols_segundas + (1|ring/año_recaptura) , na.action = na.exclude, data = volantones, family = binomial)
summary(mod1_bis) # Singular fit!

## Residuals
simulationOutput <- simulateResiduals(fittedModel = mod1_bis) 
plot(simulationOutput)  # Fine. 
testDispersion(mod1_bis)
plot(simulationOutput, data$vols_segundas)

# From logit scale to odds ratio
se <- sqrt(diag(vcov(mod1_bis)))

# Table of estimates with 95% CI
(tab <- cbind(Est = fixef(mod1_bis), LL = fixef(mod1_bis) - 1.96 * se, UL = fixef(mod1_bis) + 1.96 *se))
exp(tab)   



# Age
# Known age. We filtered those individuals for which the exact age was known.
elegida <- data %>% filter(chosen == 1) %>% drop_na(ageT) %>% filter(!str_detect(ageT, "min"))  
no_elegida <- data  %>% filter(chosen == 0) %>% drop_na(ageT) %>% filter(!str_detect(ageT, "min")) 

no_elegida_limpia <- no_elegida %>% filter(ring %in% elegida$ring)
no_elegida_limpia2 <- no_elegida_limpia %>% group_by(ring) %>% summarise(count = n_distinct(nest_dummy))
data_edad <- rbind(elegida, no_elegida_limpia)
data_edad$ageT <- as.numeric(data_edad$ageT)
edad <- data_edad %>% dplyr::select(2, 3, 5, 6, 9, 37)

# Model  
mod2 <- glmer(chosen ~ ageT + (1|ring/año_recaptura), data = edad, family = binomial) # Singular fit. 
summary(mod2)

# Residuals
simulationOutput <- simulateResiduals(fittedModel = mod2) 
plot(simulationOutput) 
plot(simulationOutput, edad$ageT)
testDispersion(mod2)

# From logit scale to odds ratio
se <- sqrt(diag(vcov(mod2)))
(tab <- cbind(Est = fixef(mod2), LL = fixef(mod2) - 1.96 * se, UL = fixef(mod2) + 1.96 * se))
exp(tab)   


# Ornamentation of the last owner
elegida <- data %>% filter(chosen == 1) %>% drop_na(long_plumasT)
no_elegida <- data %>% filter(chosen == 0) %>% drop_na(long_plumasT) 
no_elegida <- no_elegida %>% filter(ring %in% elegida$ring) 

data_plumas <- rbind(elegida, no_elegida)
length(unique(data_plumas$ring)) # Solo 49 machos.

plumas <- data_plumas %>% dplyr::select(2, 3, 5, 6, 9, 43)

# Model
mod5 <- glmer(chosen ~ long_plumasT + (1|ring/año_recaptura), data = plumas, family = binomial)
summary(mod5)

# Residuals
simulationOutput <- simulateResiduals(fittedModel = mod5) 
plot(simulationOutput) 
testDispersion(mod5)


# From logit scale to odds ratio
se <- sqrt(diag(vcov(mod5)))
(tab <- cbind(Est = fixef(mod5), LL = fixef(mod5) - 1.96 * se, UL = fixef(mod5) + 1.96 * se))
exp(tab)   

drop1(mod5, test = "Chisq")


# Owners' return 
elegidas <- data %>% filter(chosen == 1) %>% drop_na(dummy_destin)
no_elegidas <- data %>% filter(chosen ==0) %>% drop_na(dummy_destin)
no_elegida_limpia <- no_elegidas %>% filter(ring %in% elegidas$ring)

data_survival <- rbind(elegidas, no_elegida_limpia)
str(data_survival$dummy_destin)

survival <- data_survival %>% dplyr::select(2, 3, 5, 6, 9, 62)

# Modelo. 

mod6 <- glmer(chosen ~ dummy_destin + (1|ring/año_recaptura), data = survival, family = binomial)
summary(mod6)

# Residuals
simulationOutput <- simulateResiduals(fittedModel = mod6) 
plot(simulationOutput)  # Fine. 
testDispersion(mod6)
drop1(mod6, test = "Chisq")

se <- sqrt(diag(vcov(mod6)))
(tab <- cbind(Est = fixef(mod6), LL = fixef(mod6) - 1.96 * se, UL = fixef(mod6) + 1.96 * se))
exp(tab)   

plot(elegidas$dummy_destin)
cuenta <- elegidas %>% group_by(dummy_destin) %>% summarise(count = n())
cuenta


survival_barplot <- data.frame(porcentaje = c("32%", "68%"), tipo = c("Owner returned", "Owner did not return"))
percentage_survival <- ggplot(data=survival_barplot, aes(x=tipo, y=porcentaje, fill=tipo)) +
  geom_bar(stat="identity", position=position_dodge(),  width = 0.5)+
  geom_text(aes(label=porcentaje),  vjust=-0.3, size=3, color="black",
            position = position_dodge(0.9))+
  scale_fill_brewer(palette="Set2") + 
  theme(axis.title.x = element_blank(),
        axis.title.y= element_blank(),
        axis.text.x = element_text(face="plain", size =8),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 1, face = "bold", size = 16),
        legend.position = "none", 
        legend.title = element_text(colour="black", size=14, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.3, linetype = "solid"),
        panel.spacing.y = unit(2, "lines"))

plot(percentage_survival)

rm(list = ls())

# We clean the environment to continue with condition and feeding rate variables. This way (at least for me) is clearer.

#--Data 
data <- read.csv("chapter 2/data.csv", sep = ";", head = T )
length(unique(data$ring)) # 139 unique males.
str(data)

#--Checking variable type.
data$chosen_dummy <- as.factor(data$chosen_dummy)
data$ring <- as.factor(data$ring)
data$nest_dummy <- as.factor(data$nest_dummy)
data$visited <- as.factor(data$visited)
data$chosen <- as.factor(data$chosen)
data$fuera_del_circulo <- as.factor(data$fuera_del_circulo)
data$año_recaptura <- as.factor(data$año_recaptura)
data$YEAR_NEST <- as.factor(data$YEAR_NEST)
data$excluir <- as.factor(data$excluir)
data$EXPERIMENTACION_IRAIDA <- as.factor(data$EXPERIMENTACION_IRAIDA)
data$YEAR_NEST_rep <- as.factor(data$YEAR_NEST_rep) # Hay dos pájaros que acaban criando en el mismo
# nido (¡son hermanos!)
data$nest_dummy_rep <- as.factor(data$nest_dummy_rep)
data$rep_prospected <- as.factor(data$rep_prospected)
data$rep_inside_circle <- as.factor(data$rep_inside_circle)
data$YEAR_NEST_rep_previous <- as.factor(data$YEAR_NEST_rep_previous)
data$EXPERIMENTACION_IRAIDA_rep <- as.factor(data$EXPERIMENTACION_IRAIDA_rep)
data$tipo_macho <- as.factor(data$tipo_macho)
data$macho_a_buscarT <- as.factor(data$macho_a_buscarT)
data$macho_a_buscarT1 <- as.factor(data$macho_a_buscarT1)
data$observT <- as.factor(data$observT)
data$observT1 <- as.factor(data$observT1)
data$destino_CERT <- as.factor(data$destino_CERT)
data$next_YEAR_NEST_CERT <- as.factor(data$next_YEAR_NEST_CERT)
data$next_NEST_CERT <- as.factor(data$next_NEST_CERT)
data$same_nest_CERT <- as.factor(data$same_nest_CERT)
data$destino_sin_CERT <- as.factor(data$destino_sin_CERT)
data$next_YEAR_NEST_sin_CERT <- as.factor(data$next_YEAR_NEST_sin_CERT)
data$next_NEST_sin_CERT <- as.factor(data$next_NEST_sin_CERT)
data$same_nest_sin_CERT <- as.factor(data$same_nest_sin_CERT)
data$observT_focal <- as.factor(data$observT_focal)
data$observT1_focal <- as.factor(data$observT1_focal)
data$codigo_destin_CERT <- as.factor(data$codigo_destin_CERT)
data$dummy_destin <- as.factor(data$dummy_destin)
str(data)
data <- data %>% filter(excluir == 0)


# To calculate the body condition, we need to use the global dataset of captures. 
captures_global <- read.csv("./chapter 2/captures_global.csv", sep = ";", head = T ) # long-term dataset of captures.
str(captures_global)

# From the datum, we filter those with complete measurements of tarsus, beak and wing length to calculate size. 
data <- data[complete.cases(data$macho_a_buscarT, data$pesoT, data$alaT, data$pesoT, data$observT), ]

# We work with a subset.
elegidas <- data %>% filter(chosen == 1)
no_elegidas <- data %>% filter(chosen == 0)
no_elegidas <- no_elegidas %>% filter(ring %in% elegidas$ring)

# Make the dataframe
data_bis <- rbind(elegidas, no_elegidas)

# Dataframe to calculate size.
males <- data_bis %>% dplyr::select(9, 34, 39:42, 47)

# Only unique identities to not mess up with the joins
males_complete <- males[!duplicated(males$macho_a_buscarT), ] 

# Searching the males.
anillas <- captures_global  %>% dplyr::select(anilla, RING_YEAR)
anillas <- anillas[!duplicated(anillas$RING_YEAR),]
males_complete <- males_complete %>% left_join(anillas, by = c("macho_a_buscarT" = "RING_YEAR"))

# Remove NAs
males_complete <- males_complete %>% drop_na()

# We operate in the same way as we did in the chapter 1.
tarso <- lm(tarsoT ~ observT , data = males_complete)
pico <- lm(picoT ~ observT , data = males_complete)
ala <- lm(alaT ~ observT , data = males_complete)

summary(tarso)
summary(pico)
summary(ala)

# Join the residuals.

# Tarsus
ID_birds <- males_complete %>% 
  dplyr::select(macho_a_buscarT, tarsoT) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?
residuals_data_a <- data.frame("macho_a_buscarT" = ID_birds$macho_a_buscarT,
                               "tarso_corr" =  residuals(tarso)) #adding our size var. corrected.
males_complete <- left_join(males_complete, residuals_data_a) # join

# Wing
ID_birds <- males_complete%>% 
  dplyr::select(macho_a_buscarT, alaT) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?
residuals_data_a <- data.frame("macho_a_buscarT" = ID_birds$macho_a_buscarT,
                               "ala_cor" =  residuals(ala)) #adding our size var. corrected.
males_complete <- left_join(males_complete, residuals_data_a) # join

# Beak
ID_birds <- males_complete %>% 
  dplyr::select(macho_a_buscarT, picoT) %>% # select columns
  na.omit() # remove NAs
nrow(ID_birds) # rows left?
residuals_data_a <- data.frame("macho_a_buscarT" = ID_birds$macho_a_buscarT,
                               "beak_cor" =  residuals(pico)) #adding our size var. corrected.
males_complete <- left_join(males_complete, residuals_data_a) # join



# Principal Component Analysis
pca1 <- prcomp(males_complete[,9:11], center = TRUE, scale = TRUE)

# To see the eigenvalues, and the % variance explained
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
males_complete$tamano <- pca1$x[,1]


# Joining 
captures <- captures_global %>% arrange(RING_YEAR, (first_cap_year))
captures_focal <- captures %>% group_by(RING_YEAR) %>% slice_head(n=1) %>% filter(RING_YEAR %in% males_complete$macho_a_buscarT)
join <- captures_focal %>% dplyr::select(RING_YEAR, anilla,hora_corregida, fecha_enero)
males_complete <- males_complete %>% left_join(join, by = c("macho_a_buscarT" = "RING_YEAR"))
males_complete$hora_corregida <- as.numeric(males_complete$hora_corregida)

# Removing Nas-
analysis <- males_complete %>% drop_na()
analysis2 <- analysis %>% filter(pesoT < 102)
analysis$anilla.x <- as.factor(analysis$anilla.x)

# Condition model

condition <- lmer(pesoT ~ scale(tamano)  + scale(hora_corregida) +  scale(fecha_enero) + (1|año_recaptura) + (1|anilla.x)  , data = analysis) 
summary(condition) #  
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
check_model(condition) 

# Joining body condition index.
str(analysis)

ID_birds <- analysis %>% select(macho_a_buscarT, pesoT, tarsoT) %>% na.omit() # remove NAs
nrow(ID_birds) # rows left

residuals_data_b <- data.frame("macho_a_buscarT" = ID_birds$macho_a_buscarT,
                               "condition" =  residuals(condition))

dim(residuals_data_b)
males_complete <- left_join(males_complete, residuals_data_b) # left join


# Joining
condicion <- males_complete%>% dplyr::select(macho_a_buscarT, condition)
data_bis <- data_bis %>% left_join(condicion, by = "macho_a_buscarT")

# Simplication of the dataset 

data_condicion <- data_bis %>% dplyr::select(2, 3, 5, 6, 9, 87)

# Model
mod7 <- glmer(chosen ~ condition + (1|ring/año_recaptura) , data = data_condicion, family = binomial)
summary(mod7)

# Residuals
simulationOutput <- simulateResiduals(fittedModel = mod7) 
plot(simulationOutput)  # Fine. 
testDispersion(mod7)

# From logit scale to odd ratio
se <- sqrt(diag(vcov(mod7)))
(tab <- cbind(Est = fixef(mod7), LL = fixef(mod7) - 1.96 * se, UL = fixef(mod7) + 1.96 * se))
exp(tab)   

drop1(mod7, test = "Chisq")



# Feeding rate:
elegidas <- data %>% filter(chosen == 1) %>% drop_na(mean_feeding_rate)
no_elegidas <- data %>% filter(chosen ==0) %>% drop_na(mean_feeding_rate)
no_elegida_limpia <- no_elegidas %>% filter(ring %in% elegidas$ring)

data_cebas <- rbind(elegidas, no_elegida_limpia)

cebas <- data_cebas %>% dplyr::select(2, 3, 5, 6, 9, 36)

# Modelo. 
mod8 <- glmer(chosen ~ mean_feeding_rate + (1|ring/año_recaptura), data = cebas, family = binomial)
summary(mod8)

# Residuals
simulationOutput <- simulateResiduals(fittedModel = mod8) 
plot(simulationOutput)  # Fine. 
testDispersion(mod8)

# From logit scale to odds ratio
se <- sqrt(diag(vcov(mod8)))
(tab <- cbind(Est = fixef(mod8), LL = fixef(mod8) - 1.96 * se, UL = fixef(mod8) + 1.96 * se))
exp(tab)   

drop1(mod7, test = "Chisq")


# Success of the last clutch as categorical (0 = zero fledglings; 1 = at least one nestling fledged)
volantones <- data %>% dplyr::select(2, 3, 5, 6, 9, 13, 14) 
volantones <- volantones %>% mutate(exito = ifelse(vols_segundas == 0, 0, 1))
volantones$exito <- as.factor(volantones$exito)

# Modelo
mod1_bis2 <- glmer(chosen ~ exito + (1|ring/año_recaptura) , na.action = na.exclude, data = volantones, family = binomial)
summary(mod1_bis2) # Singular fit!

# Residuals
simulationOutput <- simulateResiduals(fittedModel = mod1_bis2) 
plot(simulationOutput)  # Fine. 
testDispersion(mod1_bis2)
plot(simulationOutput, data$vols_año_prospeccion)


# From logit to odds ratio.
se <- sqrt(diag(vcov(mod1_bis2)))
(tab <- cbind(Est = fixef(mod1_bis2), LL = fixef(mod1_bis2) - 1.96 * se, UL = fixef(mod1_bis2) + 1.96 *
                se))
exp(tab)   

