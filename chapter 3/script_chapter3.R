
# This R script contains the code for running a multievent capture-recapture model.
#In this multievent capture-recapture model we examine:
# breeding status related differences in survival probability, 
# probability of transition between breeding states 
# and recapture probability of 
# males and females from a wild population of spotless starlings. 


# Our work is based on the CMR workshop and materials delivered by 
# Olivier Gimenez, Chloe Nater, Sarah Cubaynes, Maraud Quéroué and Perry de Valpine 
# https://oliviergimenez.github.io/bayesian-cr-workshop/

# Recommendations:
# 1 - opening the document outline (Ctrl + Shift + O)
# 2 - having the matrices near! 


##--Libraries--------------------------------------------------------------------
library(tidyverse)
library(nimble) 
library(MCMCvis)
library(abind)
library(coda)


setwd("/mnt/lustre/scratch/nlsas/home/csic/ecf/dgi/Aberdeen_New/Aberdeen/IMEDEA/FINAL_ANALYSIS")


##--Load data--------------------------------------------------------------------
starlings <- read.csv("./data/event_hist_new_cert1.csv", sep = ",", head = T) # dataset with the event histories
#sex <- starlings$sex # vector for sex; 1 = males; 2 = females
sex_dummy <- readRDS("./data/sex_new_vector.rds")
starlings <- starlings[,-c(1:4)] # I get rid of useless columns
names(starlings)<-sapply(str_remove_all(colnames(starlings),"X"),"[")
names(starlings)
y <- as.matrix(starlings) #dataset as a matrix named 'y'.
head(y)



##--1. Model code --------------------------------------------------------------------
multievent_random_simplified <- nimbleCode({  
  # --Define states (z) :-------------------------------------------
  # 1 Alive as male nestling
  # 2 Alive as female nestling
  # 3 Alive as male breeder
  # 4 Alive as female breeder
  # 5 Alive as male floater
  # 6 Alive as female floater
  # 7 Dead
  # --Define events/observataions (y -> y+1):-------------------------------------------
  # 0 = Not seen         -> 1      # y+1 is necessary as nimble doesn't like zeros
  # 1 = Seen as nestling -> 2
  # 2 = Seen as breeder  -> 3
  # 3 = Seen as floater  -> 4
  # 4 = Seen as unknown  -> 5
  
  
  # IS state probabilities 
  
  
  delta[1,1] <- 1    # prob. of being in initial state: male nestling
  delta[1,2] <- 0    # prob. of being in initial state: female nestling
  delta[1,3] <- 0    # prob. of being in initial state: male breeder  
  delta[1,4] <- 0    # prob. of being in initial state: female breeder
  delta[1,5] <- 0    # prob. of being in initial state: male floater
  delta[1,6] <- 0    # prob. of being in initial state: female floater
  delta[1,7] <- 0    # prob. of being in initial state: dead
  delta[2,1] <- 0    # prob. of being in initial state: male nestling
  delta[2,2] <- 1    # prob. of being in initial state: female nestling
  delta[2,3] <- 0    # prob. of being in initial state: male breeder
  delta[2,4] <- 0    # prob. of being in initial state: female breeder
  delta[2,5] <- 0    # prob. of being in initial state: male floater
  delta[2,6] <- 0    # prob. of being in initial state: female floater
  delta[2,7] <- 0    # prob. of being in initial state: dead
  
  # ├ Transition matrix (Pr = probability) ----
  
  for(t in 1:(T-1)){
    
    # Note: I have specified one beta for each of the phi and psi parameters. Betas are the mean value of 
    # phi or psi. Then we add a variation as a random year effect, represented by eps (~ epsilon). Eps is drawn for 
    # same distribution for all reproductive status and sexes.
    # We have not specified sex-specific random effects. 
    # Since we specify survival and transition as function of two parameters, we must use the logit link function.
    
    # IMPORTANT to keep track of the number inside the brackets to avoid making mistakes in the prior and transformation stages. 
    
    # beta[1]  is mean phi for male juveniles.
    # beta[2]  is mean psi for the transition from juvenile to breeder in males.
    # beta[3]  is mean phi for female juveniles.
    # beta[4]  is mean psi for the transition from juvenile to breeder in females.
    # beta[5]  is mean phi for male breeders.
    # beta[6]  is mean psi for the transition from breeder to floater in males.
    # beta[7]  is mean phi for female breeders.
    # beta[8]  is mean psi for the transition from breeder to floater in females.
    # beta[9]  is mean phi for male floaters.
    # beta[10] is mean psi for the transition from floater to breeder in males.
    # beta[11] is mean phi for female floaters.
    # beta[12] is mean psi for the transition from floater to breeder in females.
    
    
    logit(phiJ_male[t])    <- beta[1]  + eps_phi[t]
    logit(psiJB_male[t])   <- beta[2]  + eps_psi[t]
    
    logit(phiJ_female[t])  <- beta[3]  + eps_phi[t]
    logit(psiJB_female[t]) <- beta[4]  + eps_psi[t]
    
    logit(phiB_male[t])    <- beta[5]  + eps_phi[t]
    logit(psiBF_male[t])   <- beta[6]  + eps_psi[t]
    
    logit(phiB_female[t])  <- beta[7]  + eps_phi[t]
    logit(psiBF_female[t]) <- beta[8]  + eps_psi[t]
    
    logit(phiF_male[t])    <- beta[9]  + eps_phi[t]
    logit(psiFB_male[t])   <- beta[10] + eps_psi[t]
    
    logit(phiF_female[t])  <- beta[11] + eps_phi[t]
    logit(psiFB_female[t]) <- beta[12] + eps_psi[t]
    
    # Probabilities of state z(t+1) given z(t) (lines 47-55)
    
    gamma[1,1,t] <-  0                                        # Pr. of being alive as a nestling male at (t)   -> Pr. of being alive as a nestling male at (t+1)
    gamma[1,2,t] <-  0                                        # Pr. of being alive as a nestling male at (t)   -> Pr. of being alive as a nestling female at (t+1)
    gamma[1,3,t] <-  phiJ_male[t] * psiJB_male[t]             # Pr. of being alive as a nestling male at (t)   -> Pr. of being alive as a male breeder at (t+1)  
    gamma[1,4,t] <-  0                                        # Pr. of being alive as a nestling male at (t)   -> Pr. of being alive as a female breeder at (t+1) 
    gamma[1,5,t] <-  phiJ_male[t] * (1 - psiJB_male[t])       # Pr. of being alive as a nestling male at (t)   -> Pr. of being alive as a male floater at (t+1)
    gamma[1,6,t] <-  0                                        # Pr. of being alive as a nestling male at (t)   -> Pr. of being alive as a female floater at (t+1)
    gamma[1,7,t] <-  1 - phiJ_male[t]                         # Pr. of being alive as a nestling male at (t)   -> Pr. of being dead at (t+1) 
    gamma[2,1,t] <-  0                                        # Pr. of being alive as a nestling female at (t) -> Pr. of being alive as a nestling male at (t+1)
    gamma[2,2,t] <-  0                                        # Pr. of being alive as a nestling female at (t) -> Pr. of being alive as a nestling female at (t+1)
    gamma[2,3,t] <-  0                                        # Pr. of being alive as a nestling female at (t) -> Pr. of being alive as a male breeder at (t+1)
    gamma[2,4,t] <-  phiJ_female[t] * psiJB_female[t]         # Pr. of being alive as a nestling female at (t) -> Pr. of being alive as a female breeder at (t+1) 
    gamma[2,5,t] <-  0                                        # Pr. of being alive as a nestling female at (t) -> Pr. of being alive as a male floater at (t+1)  
    gamma[2,6,t] <-  phiJ_female[t] * (1 - psiJB_female[t])   # Pr. of being alive as a nestling female at (t) -> Pr. of being alive as a female floater at (t+1)
    gamma[2,7,t] <-  1 - phiJ_female[t]                       # Pr. of being alive as a nestling female at (t) -> Pr. of being dead at (t+1) 
    gamma[3,1,t] <-  0                                        # Pr. of being alive as a breeder male at (t)    ->Pr. of being alive as a nestling male at (t+1)
    gamma[3,2,t] <-  0                                        # Pr. of being alive as a breeder male at (t)    -> Pr. of being alive as a nestling female at (t+1)
    gamma[3,3,t] <-  phiB_male[t] * (1 - psiBF_male[t])       # Pr. of being alive as a breeder male at (t)    -> Pr. of being alive as a male breeder at (t+1)
    gamma[3,4,t] <-  0                                        # Pr. of being alive as a breeder male at (t)    -> Pr. of being alive as a female breeder at (t+1) 
    gamma[3,5,t] <-  phiB_male[t] * psiBF_male[t]             # Pr. of being alive as a breeder male at (t)    -> Pr. of being alive as a male floater at (t+1) 
    gamma[3,6,t] <-  0                                        # Pr. of being alive as a breeder male at (t)    -> Pr. of being alive as a female floater at (t+1)
    gamma[3,7,t] <-  1 - phiB_male[t]                         # Pr. of being alive as a breeder male at (t)    ->Pr. of being dead at (t+1) 
    gamma[4,1,t] <-  0                                        # Pr. of being alive as a breeder female at (t)  -> Pr. of being alive as a nestling male at (t+1)
    gamma[4,2,t] <-  0                                        # Pr. of being alive as a breeder female at (t)  -> Pr. of being alive as a nestling female at (t+1)
    gamma[4,3,t] <-  0                                        # Pr. of being alive as a breeder female at (t)  -> Pr. of being alive as a male breeder at (t+1)
    gamma[4,4,t] <-  phiB_female[t] * (1 - psiBF_female[t])   # Pr. of being alive as a breeder female at (t)  -> Pr. of being alive as a female breeder at (t+1)    
    gamma[4,5,t] <-  0                                        # Pr. of being alive as a breeder female at (t)  -> Pr. of being alive as a male floater at (t+1)
    gamma[4,6,t] <-  phiB_female[t] * psiBF_female[t]         # Pr. of being alive as a breeder female at (t)  -> Pr. of being alive as a female floater at (t+1)
    gamma[4,7,t] <-  1 - phiB_female[t]                       # Pr. of being alive as a breeder female at (t)  -> Pr. of being dead at (t+1) 
    gamma[5,1,t] <-  0                                        # Pr. of being alive as a floater male at (t)    -> Pr. of being alive as a nestling male at (t+1)
    gamma[5,2,t] <-  0                                        # Pr. of being alive as a floater male at (t)    -> Pr. of being alive as a nestling female at (t+1)   
    gamma[5,3,t] <-  phiF_male[t] * psiFB_male[t]             # Pr. of being alive as a floater male at (t)    -> Pr. of being alive as a male breeder at (t+1) 
    gamma[5,4,t] <-  0                                        # Pr. of being alive as a floater male at (t)    -> Pr. of being alive as a female breeder at (t+1)
    gamma[5,5,t] <-  phiF_male[t] * (1 - psiFB_male[t])       # Pr. of being alive as a floater male at (t)    -> Pr. of being alive as a male floater at (t+1)
    gamma[5,6,t] <-  0                                        # Pr. of being alive as a floater male at (t)    -> Pr. of being alive as a female floater at (t+1)  
    gamma[5,7,t] <-  1 - phiF_male[t]                         # Pr. of being alive as a floater male at (t)    -> Pr. of being dead at (t+1) 
    gamma[6,1,t] <-  0                                        # Pr. of being alive as a floater female at (t)  -> Pr. of being alive as a nestling male at (t+1)
    gamma[6,2,t] <-  0                                        # Pr. of being alive as a floater female at (t)  -> Pr. of being alive as a nestling female at (t+1) 
    gamma[6,3,t] <-  0                                        # Pr. of being alive as a floater female at (t)  -> Pr. of being alive as a male breeder at (t+1)
    gamma[6,4,t] <-  phiF_female[t] * psiFB_female[t]         # Pr. of being alive as a floater female at (t)  -> Pr. of being alive as a female breeder at (t+1)
    gamma[6,5,t] <-  0                                        # Pr. of being alive as a floater female at (t)  -> Pr. of being alive as a male floater at (t+1)
    gamma[6,6,t] <-  phiF_female[t] * (1 - psiFB_female[t])   # Pr. of being alive as a floater female at (t)  -> Pr. of being alive as a female floater at (t+1)     
    gamma[6,7,t] <-  1 - phiF_female[t]                       # Pr. of being alive as a floater female at (t)  -> Pr. of being dead at (t+1) 
    gamma[7,1,t] <-  0                                        # Pr. of being dead at (t)                       -> Pr. of being alive as a nestling male at (t+1)
    gamma[7,2,t] <-  0                                        # Pr. of being dead at (t)                       -> Pr. of being alive as a nestling female at (t+1) 
    gamma[7,3,t] <-  0                                        # Pr. of being dead at (t)                       -> Pr. of being alive as a male breeder at (t+1)
    gamma[7,4,t] <-  0                                        # Pr. of being dead at (t)                       -> Pr. of being alive as a female breeder at (t+1)  
    gamma[7,5,t] <-  0                                        # Pr. of being dead at (t)                       -> Pr. of being alive as a male floater at (t+1)
    gamma[7,6,t] <-  0                                        # Pr. of being dead at (t)                       -> Pr. of being alive as a female floater at (t+1)
    gamma[7,7,t] <-  1                                        # Pr. of being dead at (t)                       -> Pr. of being dead at (t+1) 
    
  }  
  
  
  # ├ Priors for the transition matrix ----
  
  # Priors for the beta terms: we draw 12 different betas from a normal distribution.
  for(i in 1:12){
    beta[i] ~ dnorm(0, sd = 1.5) # uninformative prior
  }
  
  
  # ├ Survival random effects ----
  
  # Uninformative priors for the random year effect for phi (survival)
  
  lambda_phi ~ dunif(0, 10) 
  
  for(t in 1:(T-1)){
    eps_phi[t] ~ dnorm(0, sd = lambda_phi)
  }
  
  
  # ├ Transition random effects ----
  
  # Uninformative priors for the random year effect for psi (transition)
  
  lambda_psi ~ dunif(0, 10)
  
  for(t in 1:(T-1)){
    eps_psi[t] ~ dnorm(0, sd = lambda_psi)
    
  }
  
  
  # Througout this part of the script  I recommend having lines 168-182 on sight to being able of following 
  # which betas correspond to each state.
  
  # ├ MALES ----
  
  # T before every parameter makes reference to 'transformed', as we apply the ilogit to get our estimates  
  # in the 0-1 scale.
  
  # SURVIVAL: phi
  
  # Juveniles
  
  TphiJ_male[1] <- ilogit(beta[1] + eps_phi[1])          # t = 1
  for(t in 2:5){                                          
    TphiJ_male[t] <- 0                                    # from t = 2 to t = 5 there are
  }                                                       # no input of new juveniles so
  # I fix this value to 0   
  TphiJ_male[6]  <- ilogit(beta[1] + eps_phi[6])         # t = 6
  TphiJ_male[7]  <- ilogit(beta[1] + eps_phi[7])         # t = 7
  TphiJ_male[8]  <- ilogit(beta[1] + eps_phi[8])         # t = 8
  TphiJ_male[9]  <- ilogit(beta[1] + eps_phi[9])         # t = 9
  TphiJ_male[10] <- ilogit(beta[1] + eps_phi[10])        # t = 10
  
  # Adult breeder 
  # All individuals start as chicks
  TphiB_male[1] <- 0                                    # t = 1 there are no male breeders. I fix it to 0. 
  
  for(t in 2:(T-1)){                                    
    TphiB_male[t] <- ilogit(beta[5] + eps_phi[t])      # male breeders survival from t = 2 to t = 10 
  }
  
  # Adult floater
  
  # All individuals start as chicks
  TphiF_male[1] <- 0                                    # t = 1 there are no male floaters I fix it to 0. 
  
  for(t in 2:(T-1)){
    TphiF_male[t] <- ilogit(beta[9] + eps_phi[t])       # male floaters survival from t = 2 to t = 10  from t = 2 to t = 10
  }
  
  # TRANSITION: psi
  
  #Transition FROM: Juveniles TO: Breeder
  TpsiJB_male[1] <- ilogit(beta[2] + eps_psi[1])     # This transition is only possible for the years where there 
  # are juveniles (t = 1, 6, 7, 8, 9 and 10)
  TpsiJB_male[2] <- 0                                 # Not possible. Fixed to 0.
  TpsiJB_male[3] <- 0                                 # Not possible. Fixed to 0. 
  TpsiJB_male[4] <- 0                                 # Not possible. Fixed to 0. 
  TpsiJB_male[5] <- 0                                 # Not possible. Fixed to 0. 
  
  for(t in 6:10){
    TpsiJB_male[t] <- ilogit(beta[2] + eps_psi[t])
  }
  
  # Transition FROM: Juveniles TO: Floater -> not transitioning to B but to F 
  #same assumptions as above about cohorts without new individuals. 
  TpsiJF_male[1] <- 1 - TpsiJB_male[1]
  TpsiJF_male[2] <- 0
  TpsiJF_male[3] <- 0
  TpsiJF_male[4] <- 0
  TpsiJF_male[5] <- 0
  
  for(t in 6:10){
    TpsiJF_male[t] <- 1 - TpsiJB_male[t]
  }
  
  # TRANSITION FROM: Breeder TO: Floater
  
  TpsiBF_male[1] <- 0                                # This transition is not possible at t = 1, because there are no breeders yet.
  
  # This transition is only possible from t = 2 onwards.
  for(t in 2:(T-1)){
    TpsiBF_male[t] <- ilogit(beta[6] + eps_psi[t])
  }  
  
  # TRANSITION FROM: Breeder TO: Breeder (i.e. remaining as a breeder) (substraction)
  
  TpsiBB_male[1] <- 0
  
  for(t in 2:(T-1)){
    TpsiBB_male[t] <- 1 -  TpsiBF_male[t]
  }
  
  # TRANSITION FROM: Floater  TO: Breeder
  
  TpsiFB_male[1] <- 0                                # This transition is not possible at t = 1, because there are no floaters yet.
  
  for(t in 2:(T-1)){                                 # This transition is  only possible from t=2 onwards.
    TpsiFB_male[t] <- ilogit(beta[10] + eps_psi[t])
  }  
  
  
  # TRANSITION FROM: Floater  TO: Floater -> not transitioning to F but to B (i.e. remaining as a floater) (substraction)
  TpsiFF_male[1] <- 0
  for(t in 2:(T-1)){
    TpsiFF_male[t] <- 1 -  TpsiFB_male[t]
  }
  
  
  
  # ├ FEMALES ----
  
  # SURVIVAL: phi
  
  # Juveniles
  
  TphiJ_female[1] <- ilogit(beta[3] + eps_phi[1])    # t = 1
  
  for(t in 2:5){                                      # from t = 2 to t = 5 there are
    TphiJ_female[t] <- 0                              # no input of new juveniles so
  }                                                   # I fix this value to 0   
  
  TphiJ_female[6]  <- ilogit(beta[3] + eps_phi[6])    # t = 6
  TphiJ_female[7]  <- ilogit(beta[3] + eps_phi[7])    # t = 7
  TphiJ_female[8]  <- ilogit(beta[3] + eps_phi[8])    # t = 8
  TphiJ_female[9]  <- ilogit(beta[3] + eps_phi[9])    # t = 9
  TphiJ_female[10] <- ilogit(beta[3] + eps_phi[10])   # t = 10
  
  # Adult breeder 
  
  TphiB_female[1] <- 0                                # At t = 1 there no female breeders yet. I fix it to 0
  for(t in 2:(T-1)){
    TphiB_female[t] <- ilogit(beta[7] + eps_phi[t])
  }
  
  # Adult floater
  
  TphiF_female[1] <- 0                                # At t = 1 there no female floaters yet. I fix it to 0   
  for(t in 2:(T-1)){
    TphiF_female[t] <- ilogit(beta[11] + eps_phi[t])
  }
  
  # TRANSITION: psi
  
  # FROM: Juveniles TO: Breeder
  
  TpsiJB_female[1] <- ilogit(beta[4] + eps_psi[1]) # This transition is only possible for the years where there 
  # are juveniles (t = 1, 6, 7, 8, 9 and 10)
  TpsiJB_female[2] <- 0                             # Not possible. Fixed to 0.
  TpsiJB_female[3] <- 0                             # Not possible. Fixed to 0.
  TpsiJB_female[4] <- 0                             # Not possible. Fixed to 0.
  TpsiJB_female[5] <- 0                             # Not possible. Fixed to 0.
  
  for(t in 6:10){
    TpsiJB_female[t] <- ilogit(beta[4] + eps_psi[t])
  }
  
  # FROM: Juveniles TO: Floater ->  not transitioning to B but to F (substraction)
  
  TpsiJF_female[1] <- 1 - TpsiJB_female[1]
  TpsiJF_female[2] <- 0 
  TpsiJF_female[3] <- 0
  TpsiJF_female[4] <- 0
  TpsiJF_female[5] <- 0
  
  for(t in 6:10){
    TpsiJF_female[t] <- 1 - TpsiJB_female[t]
  }
  
  # TRANSITION FROM: Breeder TO: Floater 
  
  TpsiBF_female[1] <- 0                              # This transition is not possible at t = 1, because there are no breeders yet.
  
  for(t in 2:(T-1)){                                 # This transition is possible onwards.
    TpsiBF_female[t] <- ilogit(beta[8] + eps_psi[t])
  }  
  
  # TRANSITION FROM: Breeder TO: Breeder -> not transitioning to F but to B (i.e. remaining as a breeder) (substraction)
  
  TpsiBB_female[1] <- 0
  for(t in 2:(T-1)){
    TpsiBB_female[t] <- 1 -  TpsiBF_female[t]
  }
  
  # TRANSITION FROM: Floater  TO: Breeder 
  
  TpsiFB_female[1] <- 0                             # This transition is not possible at t = 1, because there are no floaters yet.
  
  for(t in 2:(T-1)){                                # This transition is possible onwards.
    TpsiFB_female[t] <- ilogit(beta[12] + eps_psi[t])
  }  
  
  # TRANSITION FROM: Floater  TO: Floater -> not transitioning to F but to B (i.e. remaining as a floater) (substraction)
  
  TpsiFF_female[1] <- 0
  for(t in 2:(T-1)){
    TpsiFF_female[t] <- 1 -  TpsiFB_female[t]
  }
  
  # ├ Observation matrix (Pr = probability) ----
  
  #In this model, recapture probability is sex dependent but it does not change with time. 
  
  # Probabilities of y(t) given z(t) (lines 48-61) #CHECK LINE NUMBERS
  
  omega[1,1] <-  0                            # Pr of being alive as a male nestling at t -> non-detected t)
  omega[1,2] <-  1                            # Pr of being alive as a male nestling at t -> being detected as a nestling t)
  omega[1,3] <-  0                            # Pr of being alive as a male nestling at t -> being detected and ascertained as a breeder t)
  omega[1,4] <-  0                            # Pr of being alive as a male nestling at t -> being detected and ascertained as a floater t)
  omega[1,5] <-  0                            # Pr of being alive as a male nestling at t -> being detected and not ascertained   t)
  omega[2,1] <-  0                            # Pr of being alive as a female nestling at t -> non-detected t)
  omega[2,2] <-  1                            # Pr of being alive as a female nestling at t -> being detected as a nestling t)   
  omega[2,3] <-  0                            # Pr of being alive as a female nestling at t -> being detected and ascertained as a breeder t)  
  omega[2,4] <-  0                            # Pr of being alive as a female nestling at t -> being detected and ascertained as a floater t)
  omega[2,5] <-  0                            # Pr of being alive as a female nestling at t -> being detected and not ascertained   t)
  omega[3,1] <-  1 - pB_male                  # Pr of being alive as a male breeder at t -> non-detected t)
  omega[3,2] <-  0                            # Pr of being alive as a male breeder at t -> being detected as a nestling t)
  omega[3,3] <-  aB_male * pB_male            # Pr of being alive as a male breeder at t -> being detected and ascertained as a breeder t)
  omega[3,4] <-  0                            # Pr of being alive as a male breeder at t -> being detected and ascertained as a floater t)
  omega[3,5] <-  (1 - aB_male) * pB_male      # Pr of being alive as a male breeder at t -> being detected and not ascertained t)
  omega[4,1] <-  1 - pB_female                # Pr of being alive as a female breeder at t -> non-detected t)
  omega[4,2] <-  0                            # Pr of being alive as a female breeder at t -> being detected as a nestling t)
  omega[4,3] <-  aB_female * pB_female        # Pr of being alive as a female breeder at t -> being detected and ascertained as a breeder t)
  omega[4,4] <-  0                            # Pr of being alive as a female breeder at t -> being detected and ascertained as a floater t) 
  omega[4,5] <-  (1 - aB_female) * pB_female  # Pr of being alive as a female breeder at t -> being detected and not ascertained t)
  omega[5,1] <-  1 - pF_male                  # Pr of being alive as a male floater at t -> non-detected t)
  omega[5,2] <-  0                            # Pr of being alive as a male floater at t -> being detected as a nestling t)
  omega[5,3] <-  0                            # Pr of being alive as a male floater at t -> being detected and ascertained as a breeder t)
  omega[5,4] <-  aF_male * pF_male            # Pr of being alive as a male floater at t -> being detected and ascertained as a floater t) 
  omega[5,5] <-  (1 - aF_male) * pF_male      # Pr of being alive as a male floater at t -> non-detected and not ascertained t)
  omega[6,1] <-  1 - pF_female                # Pr of being alive as a female floater at t -> non-detected t)
  omega[6,2] <-  0                            # Pr of being alive as a female floater at t -> being detected as a nestling t)
  omega[6,3] <-  0                            # Pr of being alive as a female floater at t -> being detected and ascertained as a breeder t)
  omega[6,4] <-  aF_female * pF_female        # Pr of being alive as a female floater at t -> being detected and ascertained as a floater t) 
  omega[6,5] <-  (1 - aF_female) * pF_female  # Pr of being alive as a female floater at t -> non-detected and not ascertained t)
  omega[7,1] <-  1                            # Pr of being dead at t -> non-detected t)
  omega[7,2] <-  0                            # Pr of being dead at t -> being detected as a nestling t)
  omega[7,3] <-  0                            # Pr of being dead at t -> being detected and ascertained as a breeder t)
  omega[7,4] <-  0                            # Pr of being dead at t -> being detected and ascertained as a floater t) 
  omega[7,5] <-  0                            # Pr of being dead at t -> non-detected and not ascertained t) 
  
  # Inital observation matrix (Pr = probability): t = 1. 
  # pB and pF are fixed to 1 because at the first encounter, recapture prob = 1. 
  
  omega.init[1,1] <-  0                            # Pr of being alive as a male nestling at t = 1 -> non-detected t = 1)
  omega.init[1,2] <-  1                            # Pr of being alive as a male nestling at t = 1 -> being detected as a nestling t = 1)
  omega.init[1,3] <-  0                            # Pr of being alive as a male nestling at t = 1 -> being detected and ascertained as a breeder t = 1)
  omega.init[1,4] <-  0                            # Pr of being alive as a male nestling at t = 1 -> being detected and ascertained as a floater t = 1)
  omega.init[1,5] <-  0                            # Pr of being alive as a male nestling at t = 1 -> being detected and not ascertained   t = 1 )
  omega.init[2,1] <-  0                            # Pr of being alive as a female nestling at t = 1 -> non-detected t = 1)
  omega.init[2,2] <-  1                            # Pr of being alive as a female nestling at t = 1 -> being detected as a nestling t = 1)   
  omega.init[2,3] <-  0                            # Pr of being alive as a female nestling at t = 1 -> being detected and ascertained as a breeder t = 1)  
  omega.init[2,4] <-  0                            # Pr of being alive as a female nestling at t = 1 -> being detected and ascertained as a floater t = 1)
  omega.init[2,5] <-  0                            # Pr of being alive as a female nestling at t = 1 -> being detected and not ascertained  t = 1)
  omega.init[3,1] <-  0                            # Pr of being alive as a male breeder at t = 1 -> non-detected t = 1)
  omega.init[3,2] <-  0                            # Pr of being alive as a male breeder at t = 1 -> being detected as a nestling t = 1)
  omega.init[3,3] <-  aB_male                      # Pr of being alive as a male breeder at t = 1 -> being detected and ascertained as a breeder t = 1)
  omega.init[3,4] <-  0                            # Pr of being alive as a male breeder at t = 1 -> being detected and ascertained as a floater t = 1)
  omega.init[3,5] <-  1 - aB_male                  # Pr of being alive as a male breeder at t = 1 -> being detected and not ascertained t = 1)
  omega.init[4,1] <-  0                            # Pr of being alive as a female breeder at t = 1 -> non-detected t = 1)
  omega.init[4,2] <-  0                            # Pr of being alive as a female breeder at t = 1 -> being detected as a nestling t = 1)
  omega.init[4,3] <-  aB_female                    # Pr of being alive as a female breeder at t = 1 -> being detected and ascertained as a breeder t = 1)
  omega.init[4,4] <-  0                            # Pr of being alive as a female breeder at t = 1 -> being detected and ascertained as a floater t = 1) 
  omega.init[4,5] <-  1 - aB_female                # Pr of being alive as a female breeder at t = 1 -> being detected and not ascertained t = 1)
  omega.init[5,1] <-  0                            # Pr of being alive as a male floater at t = 1 -> non-detected t = 1)
  omega.init[5,2] <-  0                            # Pr of being alive as a male floater at t = 1 -> being detected as a nestling t = 1)
  omega.init[5,3] <-  0                            # Pr of being alive as a male floater at t = 1 -> being detected and ascertained as a breeder t = 1)
  omega.init[5,4] <-  aF_male                      # Pr of being alive as a male floater at t = 1-> being detected and ascertained as a floater t = 1) 
  omega.init[5,5] <-  1 - aF_male                  # Pr of being alive as a male floater at t = 1 -> non-detected and not ascertained t = 1)
  omega.init[6,1] <-  0                            # Pr of being alive as a female floater at t = 1 -> non-detected t = 1)
  omega.init[6,2] <-  0                            # Pr of being alive as a female floater at t = 1 -> being detected as a nestling t = 1)
  omega.init[6,3] <-  0                            # Pr of being alive as a female floater at t = 1 -> being detected and ascertained as a breeder t = 1)
  omega.init[6,4] <-  aF_female                    # Pr of being alive as a female floater at t = 1-> being detected and ascertained as a floater t = 1) 
  omega.init[6,5] <-  1 - aF_female                # Pr of being alive as a female floater at t = 1 -> non-detected and not ascertained t = 1)
  omega.init[7,1] <-  1                            # Pr of being dead at t = 1 -> non-detected t = 1)
  omega.init[7,2] <-  0                            # Pr of being dead at t = 1 -> being detected as a nestling t = 1)
  omega.init[7,3] <-  0                            # Pr of being dead at t = 1 -> being detected and ascertained as a breeder t = 1)
  omega.init[7,4] <-  0                            # Pr of being dead at t = 1 -> being detected and ascertained as a floater t = 1) 
  omega.init[7,5] <-  0                            # Pr of being dead at t = 1 -> non-detected and not ascertained t = 1) 
  
  # Priors for assignment and recapture. 
  # In this case, we do not use the logit so we can draw priors directly from a uniform dist. 
  
  pB_male   ~ dunif(0,1)
  pB_female ~ dunif(0,1)
  pF_male   ~ dunif(0,1)
  pF_female ~ dunif(0,1)
  aB_male   ~ dunif(0,1)
  aB_female ~ dunif(0,1)
  aF_male   ~ dunif(0,1)
  aF_female ~ dunif(0,1)
  
  # ├ Calculating differences between parameters for each time step. ----
  
  # Survival
  
  # Juveniles
  juvenile_survival_diff_male_female[1]  <- TphiJ_male[1]  - TphiJ_female[1]
  juvenile_survival_diff_male_female[2]  <- 0
  juvenile_survival_diff_male_female[3]  <- 0
  juvenile_survival_diff_male_female[4]  <- 0
  juvenile_survival_diff_male_female[5]  <- 0
  juvenile_survival_diff_male_female[6]  <- TphiJ_male[6]  - TphiJ_female[6]
  juvenile_survival_diff_male_female[7]  <- TphiJ_male[7]  - TphiJ_female[7]
  juvenile_survival_diff_male_female[8]  <- TphiJ_male[8]  - TphiJ_female[8]
  juvenile_survival_diff_male_female[9]  <- TphiJ_male[9]  - TphiJ_female[9]
  juvenile_survival_diff_male_female[10] <- TphiJ_male[10] - TphiJ_female[10]
  
  # Adults
  
  # Breeders: male survival vs. female survival. 
  
  diff_survival_male_female_breeders[1]  <- 0
  
  for(t in 2:(T-1)){
    diff_survival_male_female_breeders[t]  <- TphiB_male[t]  - TphiB_female[t]
  }
  
  
  # Floaters: male survival vs. female survival.
  
  diff_survival_male_female_floaters[1] <- 0
  
  for(t in 2:(T-1)){
    diff_survival_male_female_floaters[t]  <- TphiF_male[t]  - TphiF_female[t]
  }
  
  
  # Males: breeder vs. floater survival. 
  
  diff_survival_malesB_malesF[1]  <- 0
  
  for(t in 2:(T-1)){
    diff_survival_malesB_malesF[t]  <- TphiB_male[t]  - TphiF_male[t]
  }
  
  
  # Females: breeder vs. floater survival. 
  diff_survival_femalesB_femalesF[1]  <- 0
  
  for(t in 2:(T-1)){
    diff_survival_femalesB_femalesF[t]  <- TphiB_female[t]  - TphiF_female[t]
  }
  
  
  # Transition
  # FROM: Juvenile -> Breeder; male vs. female.
  
  juvenile_transition_to_breed_male_female[1]  <- TpsiJB_male[1]  - TpsiJB_female[1]
  juvenile_transition_to_breed_male_female[2]  <- 0
  juvenile_transition_to_breed_male_female[3]  <- 0
  juvenile_transition_to_breed_male_female[4]  <- 0
  juvenile_transition_to_breed_male_female[5]  <- 0
  juvenile_transition_to_breed_male_female[6]  <- TpsiJB_male[6]  - TpsiJB_female[6]
  juvenile_transition_to_breed_male_female[7]  <- TpsiJB_male[7]  - TpsiJB_female[7]
  juvenile_transition_to_breed_male_female[8]  <- TpsiJB_male[8]  - TpsiJB_female[8]
  juvenile_transition_to_breed_male_female[9]  <- TpsiJB_male[9]  - TpsiJB_female[9]
  juvenile_transition_to_breed_male_female[10] <- TpsiJB_male[10] - TpsiJB_female[10]
  
  
  # FROM: Juvenile -> Floater; male vs. female.
  
  juvenile_transition_to_float_male_female[1]  <- TpsiJF_male[1]  - TpsiJF_female[1]
  juvenile_transition_to_float_male_female[2]  <- 0
  juvenile_transition_to_float_male_female[3]  <- 0
  juvenile_transition_to_float_male_female[4]  <- 0
  juvenile_transition_to_float_male_female[5]  <- 0
  juvenile_transition_to_float_male_female[6]  <- TpsiJF_male[6]  - TpsiJF_female[6]
  juvenile_transition_to_float_male_female[7]  <- TpsiJF_male[7]  - TpsiJF_female[7]
  juvenile_transition_to_float_male_female[8]  <- TpsiJF_male[8]  - TpsiJF_female[8]
  juvenile_transition_to_float_male_female[9]  <- TpsiJF_male[9]  - TpsiJF_female[9]
  juvenile_transition_to_float_male_female[10] <- TpsiJF_male[10] - TpsiJF_female[10]
  
  # FROM: Breeder ->  Floater; male vs. female.
  adult_transition_from_breeder_to_floater_male_female[1]  <- 0
  
  for(t in 2:(T-1)){
    adult_transition_from_breeder_to_floater_male_female[t]  <- TpsiBF_male[t]  - TpsiBF_female[t]    
  }
  
  
  # FROM: Floater -> Breeder; male vs. female.
  adult_transition_from_floater_to_breeder_male_female[1]  <- 0
  
  for(t in 2:(T-1)){
    adult_transition_from_floater_to_breeder_male_female[t]  <- TpsiFB_male[t]  - TpsiFB_female[t]    
  }
  
  
  # FROM: Breeder -> Breeder: male vs. female
  adult_remaining_as_breeder_male_female[1]  <- 0
  
  for(t in 2:(T-1)){
    adult_remaining_as_breeder_male_female[t]  <- TpsiBB_male[t]  - TpsiBB_female[t]    
  }
  
  
  # FROM: Floater -> Floater: male vs. female
  adult_remaining_as_floater_male_female[1]  <- 0
  
  for(t in 2:(T-1)){
    adult_remaining_as_floater_male_female[t]  <- TpsiFF_male[t]  - TpsiFF_female[t]    
  }
  
  
  # Recapture:
  
  # Males: breeder vs. floater
  diff_recapture_malesB_malesF     <- pB_male   - pF_male
  
  # Females: breeder vs. floater
  diff_recapture_femalesB_femalesF <- pB_female - pF_female
  
  # Breeders: male vs. females
  diff_recapture_male_female_breeders <- pB_male - pB_female
  
  # Floaters: male vs. females
  diff_recapture_male_female_floaters <- pF_male - pF_female
  
  
  # likelihood 
  for (i in 1:N){
    # latent state at first capture
    z[i,first[i]] ~ dcat(delta[sex[i], 1:7]) # depending of the sex, each individual will begin at the first (males) or second (females) row of this vector.
    y[i,first[i]] ~ dcat(omega.init[z[i,first[i]], 1:5])
    for (t in (first[i]+1):T){
      # z(t) given z(t-1)
      z[i,t] ~ dcat(gamma[z[i,t-1], 1:7, t-1])  
      # y(t) given z(t)
      y[i,t] ~ dcat(omega[z[i,t],1:5])         
    }
  }
})

## --2. Get the data of first capture. ----------------------------------------------------------------------------------------
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
first # this vector tells us the first occasion of capture. Its length = number of individuals.


## --3. Data in a list. ---------------------------------------------------------------------------------------
my.data <- list(y = y + 1) # We add 1 to the whole matrix to get rid of the 0s,
# If you take a look at the likelihood section of the model, we use the categorical
# distribution: dcat(). 

# This distribution cannot work with zeros, so we must add 1 to our capture-recapture histories


## --4. Initial values.----------------------------------------------------------------------------------------

# In this section we have to initialize the model. 
#set zinit equal to the observations
zinit <- y  
sex <- sex_dummy

# In our dataset, we have a lot of individuals with unknown sex. Using the data from a whole cohort that was sexed and which 
# showed a sex ratio of 1:1. We use that value to randomly assign sex to individuals to initialize the model.  
# sex[sex==3] <- sample(c(1,2), sum(sex_dummy==3), replace = TRUE, prob = c(1/2,1/2))

# This will replace individuals of unknown sex for 1 (males) or 2 (females). To avoid having a different assignation of sex every time
# we run this model, we can set a seed or store the sex vector and use it later. This is very useful for the part of the 
# Posterior Predictive Checks, where we must use the same data set that is used in nimbleMCMC.

for(i in 1:nrow(y)){                 # for each individual from 1 to 7158
  if(sex[i]== 1){                    # if it is a MALE
    for(j in 1:ncol(y)){             # for each column (= occasion) from 1 to 11
      if (j  < first[i])              {zinit[i,j] <- 0} # if the occasion is previous to the first encounter of the ind <- put in a 0 (an individual cannot be seen before their first encounter) 
      if (j == first[i])              {zinit[i,j] <- 1} # if the occasion is the same as the first encounter occasion < put in a 1 (seen as a MALE nestling)
      if (j  > first[i] & y[i,j] == 0){zinit[i,j] <- which(rmultinom(1, 1, c(0, 0, 1/2, 0, 1/2, 0))==1)}# if the occasion is posterior to the first encounter and there is a 0, 
      # We should draw a possible latent state for adult males (maybe they were not seen but  
      # their latent state was breeder or floater). The position on the vector in the multinom 
      # represents this: position 3 = male breeder and position 5 = male floater. Lines 47-55. 
      if (j  > first[i] & y[i,j] ==2) {zinit[i,j] <- 3}                                                  # Since our capture recapture histories are coded as the observations of breeding status and not sex,
      # We have to make a correspondence between those males that are breeders and their observation in our data only as 'breeders'.
      # For example: we only have 5 events (lines 56-61). The information about the sex is stored in a vector.
      # If we know that individual 85 was seen as a breeder (event = 2 in our) and that is was male (sex = 1), his latent state (z) is equal 
      # to = 3 (lines 47-55). The same goes for the rest of events and states.  
      if (j  > first[i] & y[i,j] ==3) {zinit[i,j] <- 5}
      if (j  > first[i] & y[i,j] ==4) {zinit[i,j] <- which(rmultinom(1, 1, c(0, 0, 1/2, 0, 1/2, 0))==1)}   # And if the individual was seen with unknown breeding status, the latent state could be either breeder or floater! We must draw one!
    }
  }else{    # if the sex is not equal to 1, ergo, is equal to 2 (females), we do the same but changing the numbers that correspond to the latent states of females. 
    for(j in 1:ncol(y)){
      if (j  < first[i])              {zinit[i,j] <- 0}
      if (j == first[i])              {zinit[i,j] <- 2}
      if (j  > first[i] & y[i,j] ==0) {zinit[i,j] <- which(rmultinom(1, 1, c(0, 0, 0, 1/2, 0, 1/2))==1)}
      if (j  > first[i] & y[i,j] ==2) {zinit[i,j] <- 4}
      if (j  > first[i] & y[i,j] ==3) {zinit[i,j] <- 6}
      if (j  > first[i] & y[i,j] ==4) {zinit[i,j] <- which(rmultinom(1, 1, c(0, 0, 0, 1/2, 0, 1/2))==1)}
    } 
  }
  
}

zinit <- as.matrix(zinit)

## --5. Constants in a list.----------------------------------------------------------------------------------------

# Constants are those elements that can never be changed and must be provided when a model is defined. 
my.constants <- list(first = first, # first encounter occasion
                     T = ncol(y),   # number of time occasions
                     N = nrow(y),   # number of individuals
                     sex = sex)     # vector of sex. At this point there are no individuals of unknown sex, remember!


# Initial values for the model.
# Special attention to the distribution of the parameters.
set.seed(1375) # for reproducibility
initial.values <- list(beta = rnorm(12, 0, 1),
                       eps_phi = rep(0, my.constants$T-1),
                       eps_psi = rep(0, my.constants$T-1),
                       lambda_phi = runif(1, 0, 3),
                       lambda_psi = runif(1, 0, 3),
                       pB_male   = runif(1, 0, 1),
                       pB_female = runif(1, 0, 1),
                       pF_male   = runif(1, 0, 1),
                       pF_female = runif(1, 0, 1),
                       aB_male   = runif(1, 0, 1),
                       aB_female = runif(1, 0, 1),
                       aF_male   = runif(1, 0, 1),
                       aF_female = runif(1, 0, 1),
                       z = zinit)



# We build the model with nimbleModel(). This creates an uncompiled object.
survival <- nimbleModel(multievent_random_simplified, data = list(y = y+1), 
                        constants = list(N = nrow(y), T = ncol(y), first = first, sex = sex),
                        inits = initial.values)

Cmodel <- compileNimble(survival)
survivalConf <- configureMCMC(survival, monitors = c("beta", "eps_phi", "eps_psi",
                                                     "TphiJ_male",
                                                     "TphiJ_female",
                                                     "TphiB_male",
                                                     "TphiB_female",
                                                     "TphiF_male",
                                                     "TphiF_female",
                                                     "TpsiJB_male",
                                                     "TpsiJB_female",
                                                     "TpsiBF_male",
                                                     "TpsiBF_female",
                                                     "TpsiFB_male",
                                                     "TpsiFB_female",
                                                     "TpsiJF_male",
                                                     "TpsiBB_male",
                                                     "TpsiFF_male",
                                                     "TpsiJF_female",
                                                     "TpsiBB_female",
                                                     "TpsiFF_female",
                                                     "lambda_phi",
                                                     "lambda_psi",
                                                     "pB_male",
                                                     "pB_female",
                                                     "pF_male",
                                                     "pF_female",
                                                     "aB_male",
                                                     "aB_female",
                                                     "aF_male",
                                                     "aF_female"))  

# Choosing the sampler and compiling the model.

survivalConf$removeSamplers(c("beta", "eps_phi", "eps_psi"))
survivalConf$addSampler(target = c("beta", "eps_phi", "eps_psi"), type = "AF_slice")
survivalMCMC <- buildMCMC(survivalConf)
CsurvivalMCMC <- compileNimble(survivalMCMC, 
                               project = survival)


# Running the model. 
samples <- runMCMC(mcmc = CsurvivalMCMC, 
                   niter = 30000, 
                   nburnin = 5000, 
                   nchains = 2)


saveRDS(samples, "./output.rds") #output in rds format
