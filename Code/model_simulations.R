rm(list=ls())

# load packages
library(dplyr)
library(ggplot2)
library(marginaleffects)
library(ggpubr)

# define palettes
palette_landscape <- c("#4ea33a80","#357029ff","#307ecd80","#1d4d7cff")
palette_community <- c("black","#4da43aff","#2f7eceff")
col_green <- "#4da43aff"
col_blue <- "#2f7eceff"

# MODEL PARAMETERS ----

# import dataframe with model parameters
# determined from parametrization experiments
model_parameters <- read.csv("Data/model_parameters.csv")

# BB
r1_CI = as.numeric(model_parameters[model_parameters$parameter=="r1",3:4]) # growth rate confidence interval
alpha11_CI = as.numeric(model_parameters[model_parameters$parameter=="alpha11",3:4]) # intraspecific competition confidence interval
alpha12_CI = as.numeric(model_parameters[model_parameters$parameter=="alpha12",3:4]) # interspecific competition confidence interval
e1_CI = as.numeric(model_parameters[model_parameters$parameter=="e1",3:4]) # emigration rate confidence interval
Ae1 = as.numeric(model_parameters[model_parameters$parameter=="Ae1",2]) # minimum density for emigration
beta1_CI = as.numeric(model_parameters[model_parameters$parameter=="beta1",3:4]) # parasitoid rate confidence interval (type I response)

# LE
r2_CI = as.numeric(model_parameters[model_parameters$parameter=="r2",3:4]) # growth rate confidence interval
alpha22_CI = as.numeric(model_parameters[model_parameters$parameter=="alpha22",3:4]) # intraspecific competition confidence interval
alpha21_CI = as.numeric(model_parameters[model_parameters$parameter=="alpha21",3:4]) # interspecific competition confidence interval
e2_CI = as.numeric(model_parameters[model_parameters$parameter=="e2",3:4]) # emigration rate confidence interval
Ae2 = as.numeric(model_parameters[model_parameters$parameter=="Ae2",2]) # minimum density for emigration
beta2_CI = as.numeric(model_parameters[model_parameters$parameter=="beta2",3:4]) # parasitoid rate confidence interval (type I response)

# DR
pmax = as.numeric(model_parameters[model_parameters$parameter=="pmax",2]) # maximum parasitization
tau = round(as.numeric(model_parameters[model_parameters$parameter=="tau",2]),0) # time to emergence
f = as.numeric(model_parameters[model_parameters$parameter=="f",2]) # fraction of females
lambda = 7 # lifespan 
eP = as.numeric(model_parameters[model_parameters$parameter=="eP",2]) # emigration rate
Pe = as.numeric(model_parameters[model_parameters$parameter=="Pe",2]) # minimum density for emigration

# simulation inputs
n_rep = 100 # number of replicas
dt = 1 # timestep size
tmax = 26 # maxmimum number of timesteps
M_adj = matrix(c(0,1,1,1,1,
                 1,0,0,0,0,
                 1,0,0,0,0,
                 1,0,0,0,0,
                 1,0,0,0,0), 5,5) # patch connectivity matrix
A10 = 10 # initial number of BB
A20 = 10 # initial number of LE
P0 = 1 # initial number of DR


# FUNCTIONS ----

# dynamics functions

# single aphid species
dt_one_aphid = function(r, alpha, Ae, e, A0, 
                        dt, tmax, patch_state, M_adj){
  
  # time vector
  time = seq(0,tmax,dt)
  
  # number of time steps
  n_dt = length(time)
  
  # nummber of patches
  n_patch = nrow(M_adj)
  
  # number of connections of patches
  c = rowSums(M_adj)
  
  # initialise aphid density matrix
  A = matrix(0,n_dt,n_patch)
  A[1,patch_state==1] = A0
  
  # initialise matrix of community and dispersal contributions
  dA_comm = matrix(0,n_dt,n_patch)
  dA_disp = matrix(0,n_dt,n_patch)
  
  for(t in 1:(n_dt-1)){
    
    for(k in 1:n_patch){
      
      # GROWTH & COMPETITION
      
      # change in A due to growth
      dA_growth_comp = A[t,k]*exp((r-alpha*A[t,k])*dt) - A[t,k]
      
      # EMIGRATION
      
      # delta for emigration
      delta = ifelse(A[t,k]<Ae, 0, 1)
      
      # change in A due to emigration
      dA_emig = delta*e*(A[t,k]-Ae) * dt
      
      # IMIGRATION
      
      # neighbours of patch k
      k_neigh = which(M_adj[k,]==1)
      
      # initialise dA due to immigration
      dA_immig = 0
      
      # loop through neighbouring patches
      for(n in k_neigh){
        
        # delta of neighbouring patch
        delta_neigh = ifelse(A[t,n]<Ae, 0, 1)
        
        # update dA due to immigration
        dA_immig = dA_immig + delta_neigh*e*(A[t,n]-Ae)/c[n] * dt
      }
      
      # new population size
      A[t+1,k] = A[t,k] + dA_growth_comp - dA_emig + dA_immig
      
      if(is.nan(A[t+1,k])){
        print("Error: A=NaN")
        break}
      if(A[t+1,k]<0){A[t+1,k]=0}
      
      # store community and dispersal contributions
      dA_comm[t+1,k] = dA_growth_comp
      dA_disp[t+1,k] = - dA_emig + dA_immig
      
    } # k loop
    
  } # t loop
  
  # output dataframe
  df_out = data.frame(expand.grid(t=time, patch=1:n_patch),
                      population_size=as.vector(A),
                      dN_comm=as.vector(dA_comm),
                      dN_disp=as.vector(dA_disp))
  
  return(df_out)
}

# two aphid species
dt_two_aphid = function(r1, r2, alpha11, alpha22, alpha12, alpha21, 
                        Ae1, Ae2, e1, e2, A10, A20,
                        dt, tmax, patch_state, M_adj){
  
  # time vector
  time = seq(0,tmax,dt)
  
  # number of time steps
  n_dt = length(time)
  
  # nummber of patches
  n_patch = nrow(M_adj)
  
  # number of connections of patches
  c = rowSums(M_adj)
  
  # initialise aphid density matrix
  A1 = matrix(0,n_dt,n_patch)
  A1[1,patch_state==1] = A10
  A2 = matrix(0,n_dt,n_patch)
  A2[1,patch_state==1] = A20
  
  # initialise matrix of community and dispersal contributions
  dA1_comm = matrix(0,n_dt,n_patch)
  dA1_disp = matrix(0,n_dt,n_patch)
  dA2_comm = matrix(0,n_dt,n_patch)
  dA2_disp = matrix(0,n_dt,n_patch)
  
  for(t in 1:(n_dt-1)){
    
    for(k in 1:n_patch){
      
      # GROWTH & COMPETITION
      
      # change in A due to growth
      dA1_growth_comp = A1[t,k]*exp((r1-alpha11*A1[t,k]-alpha12*A2[t,k])*dt) - A1[t,k]
      dA2_growth_comp = A2[t,k]*exp((r2-alpha22*A2[t,k]-alpha21*A1[t,k])*dt) - A2[t,k]
      
      # EMIGRATION
      
      # delta for emigration
      delta1 = ifelse(A1[t,k]<Ae1, 0, 1)
      delta2 = ifelse(A2[t,k]<Ae2, 0, 1)
      
      # change in A due to emigration
      dA1_emig = delta1*e1*(A1[t,k]-Ae1) * dt
      dA2_emig = delta2*e2*(A2[t,k]-Ae2) * dt
      
      # IMIGRATION
      
      # neighbours of patch k
      k_neigh = which(M_adj[k,]==1)
      
      # initialise dA due to immigration
      dA1_immig = 0
      dA2_immig = 0
      
      # loop through neighbouring patches
      for(n in k_neigh){
        
        # delta of neighbouring patch
        delta1_neigh = ifelse(A1[t,n]<Ae1, 0, 1)
        delta2_neigh = ifelse(A2[t,n]<Ae2, 0, 1)
        
        # update dA due to immigration
        dA1_immig = dA1_immig + delta1_neigh*e1*(A1[t,n]-Ae1)/c[n] * dt
        dA2_immig = dA2_immig + delta2_neigh*e2*(A2[t,n]-Ae2)/c[n] * dt
      }
      
      # new population size
      A1[t+1,k] = A1[t,k] + dA1_growth_comp - dA1_emig + dA1_immig
      A2[t+1,k] = A2[t,k] + dA2_growth_comp - dA2_emig + dA2_immig
      
      if(is.nan(A1[t+1,k]) || is.nan(A2[t+1,k])){
        print("Error: A=NaN")
        break}
      if(A1[t+1,k]<0){A1[t+1,k]=0}
      if(A2[t+1,k]<0){A2[t+1,k]=0}
      
      # store community and dispersal contributions
      dA1_comm[t+1,k] = dA1_growth_comp
      dA1_disp[t+1,k] = - dA1_emig + dA1_immig
      dA2_comm[t+1,k] = dA2_growth_comp
      dA2_disp[t+1,k] = - dA2_emig + dA2_immig
      
    } # k loop
    
  } # t loop
  
  # output dataframe
  df_out = rbind(data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="BB",
                            population_size=as.vector(A1),
                            dN_comm=as.vector(dA1_comm),
                            dN_disp=as.vector(dA1_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="LE",
                            population_size=as.vector(A2),
                            dN_comm=as.vector(dA2_comm),
                            dN_disp=as.vector(dA2_disp)))
  
  return(df_out)
}

# two aphid species & parasitoid
dt_two_aphid_ptoid = function(r1, r2, alpha11, alpha22, alpha12, alpha21, 
                              Ae1, Ae2, e1, e2, A10, A20,
                              beta1, beta2, pmax, f, tau, lambda, Pe, eP,
                              dt, tmax, patch_state, M_adj){
  
  # time vector
  time = seq(0,tmax,dt)
  
  # number of time steps
  n_dt = length(time)
  
  # nummber of patches
  n_patch = nrow(M_adj)
  
  # number of connections of patches
  c = rowSums(M_adj)
  
  # initialise aphid density matrix
  A1 = matrix(0,n_dt,n_patch)
  A1[1,patch_state==1] = A10
  A2 = matrix(0,n_dt,n_patch)
  A2[1,patch_state==1] = A20
  
  # initialise matrix of community and dispersal contributions
  dA1_comm = matrix(0,n_dt,n_patch)
  dA1_disp = matrix(0,n_dt,n_patch)
  dA2_comm = matrix(0,n_dt,n_patch)
  dA2_disp = matrix(0,n_dt,n_patch)
  dP_comm = matrix(0,n_dt,n_patch)
  dP_disp = matrix(0,n_dt,n_patch)
  
  # initialise parasitoid density matrix
  P = matrix(0,n_dt,n_patch)
  P[1:(lambda-1),patch_state==1] = P0
  
  # initialise number of parasitized aphids at each timestep
  A_parasit = matrix(0,n_dt,n_patch)
  
  # max number of parasitized aphids per timestep
  pmax_dt = pmax / ((lambda-1)/dt)
  
  for(t in 1:(n_dt-1)){
    
    for(k in 1:n_patch){
      
      # GROWTH & COMPETITION - APHID
      
      # change in A due to growth
      dA1_growth_comp = A1[t,k]*exp((r1-alpha11*A1[t,k]-alpha12*A2[t,k])*dt) - A1[t,k]
      dA2_growth_comp = A2[t,k]*exp((r2-alpha22*A2[t,k]-alpha21*A1[t,k])*dt) - A2[t,k]
      
      # PREDATION
      
      # potential number of parasitized aphids
      if(t*dt<tau){
        # type I response
        A1_parasit = beta1*A1[t,k]*P[t,k] * dt
        A2_parasit = beta2*A2[t,k]*P[t,k] * dt
        
      } else {
        # type I response
        A1_parasit = beta1*A1[t,k]*f*P[t,k] * dt
        A2_parasit = beta2*A2[t,k]*f*P[t,k] * dt
      }
      # check if number of parasitized aphids > total number of aphids
      if(A1_parasit > A1[t,k]){A1_parasit = A1[t,k]}
      if(A2_parasit > A2[t,k]){A2_parasit = A2[t,k]}
      
      # total number of aphids parasitized
      n_parasit = A1_parasit + A2_parasit
      
      # potential number of parasitized aphids
      parasit_max = pmax_dt*P[t,k]
      
      # check if maximum number of parasitized aphids exceeded
      if(n_parasit > parasit_max){
        # total number of parasitized aphids
        A_parasit[t,k] = parasit_max
        # change in A due to parasitism
        dA1_parasit = parasit_max * A1_parasit / (A1_parasit+A2_parasit)
        dA2_parasit = parasit_max * A2_parasit / (A1_parasit+A2_parasit)
      } else {
        # total number of parasitized aphids
        A_parasit[t,k] = n_parasit
        # change in A due to parasitism
        dA1_parasit = A1_parasit
        dA2_parasit = A2_parasit
      }
      
      # BIRTHS - PARASITOID
      
      # dP due to births
      if(t*dt>=tau){
        dP_birth = A_parasit[t-(tau-1)/dt,k]
      } else {
        dP_birth = 0
      }
      
      # DEATHS - PARASITOID
      
      # total number of ptoids across all patches
      P_total = sum(P[t,])
      
      # dP due to deaths
      if(t*dt==(lambda-1)){
        dP_death = P[t,k]
      } else if(t*dt>(tau+lambda-1) & P_total>0) {
        # total ptoids emerged at t-(lambda-1) = total deaths across all patches
        P_death_total = sum(A_parasit[t-(lambda-1)/dt,])
        # number of deaths in current patch (proportional to number of ptoids)
        dP_death = P[t,k]/P_total * P_death_total
      } else {
        dP_death = 0
      }
      
      # EMIGRATION
      
      # delta for emigration
      delta1 = ifelse(A1[t,k]<Ae1, 0, 1)
      delta2 = ifelse(A2[t,k]<Ae2, 0, 1)
      deltaP = ifelse(P[t,k]<Pe, 0, 1)
      
      # change in A due to emigration
      dA1_emig = delta1*e1*(A1[t,k]-Ae1) * dt
      dA2_emig = delta2*e2*(A2[t,k]-Ae2) * dt
      dP_emig = deltaP*eP*(P[t,k]-Pe) * dt
      
      # IMIGRATION
      
      # neighbours of patch k
      k_neigh = which(M_adj[k,]==1)
      
      # initialise dA due to immigration
      dA1_immig = 0
      dA2_immig = 0
      dP_immig = 0
      
      # loop through neighbouring patches
      for(n in k_neigh){
        
        # delta of neighbouring patch
        delta1_neigh = ifelse(A1[t,n]<Ae1, 0, 1)
        delta2_neigh = ifelse(A2[t,n]<Ae2, 0, 1)
        deltaP_neigh = ifelse(P[t,n]<Pe, 0, 1)
        
        # update dA due to immigration
        dA1_immig = dA1_immig + delta1_neigh*e1*(A1[t,n]-Ae1)/c[n] * dt
        dA2_immig = dA2_immig + delta2_neigh*e2*(A2[t,n]-Ae2)/c[n] * dt
        dP_immig = dP_immig + deltaP_neigh*eP*(P[t,n]-Pe)/c[n] * dt
      }
      
      # new population size
      A1[t+1,k] = A1[t,k] + dA1_growth_comp - dA1_parasit - dA1_emig + dA1_immig
      A2[t+1,k] = A2[t,k] + dA2_growth_comp - dA2_parasit - dA2_emig + dA2_immig
      P[t+1,k] = P[t,k] + dP_birth - dP_death - dP_emig + dP_immig
      
      if(is.nan(A1[t+1,k]) || is.nan(A2[t+1,k]) || is.nan(P[t+1,k])){
        print("Error: A=NaN")
        break}
      if(A1[t+1,k]<0){A1[t+1,k]=0}
      if(A2[t+1,k]<0){A2[t+1,k]=0}
      if(P[t+1,k]<0){P[t+1,k]=0}
      
      # store community and dispersal contributions
      dA1_comm[t+1,k] = dA1_growth_comp - dA1_parasit
      dA1_disp[t+1,k] = - dA1_emig + dA1_immig
      dA2_comm[t+1,k] = dA2_growth_comp - dA2_parasit
      dA2_disp[t+1,k] = - dA2_emig + dA2_immig
      dP_comm[t+1,k] = dP_birth - dP_death
      dP_disp[t+1,k] = - dP_emig + dP_immig
      
    } # k loop
    
  } # t loop
  
  # output dataframe
  df_out = rbind(data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="BB",
                            population_size=as.vector(A1),
                            dN_comm=as.vector(dA1_comm),
                            dN_disp=as.vector(dA1_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="LE",
                            population_size=as.vector(A2),
                            dN_comm=as.vector(dA2_comm),
                            dN_disp=as.vector(dA2_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="DR",
                            population_size=as.vector(P),
                            dN_comm=as.vector(dP_comm),
                            dN_disp=as.vector(dP_disp)))
  
  return(df_out)
}

# three aphid species & parasitoid
dt_three_aphid_ptoid = function(r1, r2, alpha11, alpha22, alpha12, alpha21, 
                                Ae1, Ae2, e1, e2, A10, A20,
                                beta1, beta2, pmax, f, tau, lambda, Pe, eP,
                                dt, tmax, patch_state, M_adj){
  
  # time vector
  time = seq(0,tmax,dt)
  
  # number of time steps
  n_dt = length(time)
  
  # nummber of patches
  n_patch = nrow(M_adj)
  
  # number of connections of patches
  c = rowSums(M_adj)
  
  # aphid 3 parameters
  A30 = (A10+A20)/2
  r3 = (r1+r2)/2
  alpha33 = (alpha11+alpha22)/2
  alpha3j = (alpha12+alpha21)/2
  alphai3 = (alpha12+alpha21)/2
  Ae3 = (Ae1+Ae2)/2
  e3 = (e1+e2)/2
  beta3 = (beta1+beta2)/2
  
  # initialise aphid density matrix
  A1 = matrix(0,n_dt,n_patch)
  A1[1,patch_state==1] = A10
  A2 = matrix(0,n_dt,n_patch)
  A2[1,patch_state==1] = A20
  A3 = matrix(0,n_dt,n_patch)
  A3[1,patch_state==1] = A30
  
  # initialise matrix of community and dispersal contributions
  dA1_comm = matrix(0,n_dt,n_patch)
  dA1_disp = matrix(0,n_dt,n_patch)
  dA2_comm = matrix(0,n_dt,n_patch)
  dA2_disp = matrix(0,n_dt,n_patch)
  dA3_comm = matrix(0,n_dt,n_patch)
  dA3_disp = matrix(0,n_dt,n_patch)
  dP_comm = matrix(0,n_dt,n_patch)
  dP_disp = matrix(0,n_dt,n_patch)
  
  # initialise parasitoid density matrix
  P = matrix(0,n_dt,n_patch)
  P[1:(lambda-1),patch_state==1] = P0
  
  # initialise number of parasitized aphids at each timestep
  A_parasit = matrix(0,n_dt,n_patch)
  
  # max number of parasitized aphids per timestep
  pmax_dt = pmax / ((lambda-1)/dt)
  
  for(t in 1:(n_dt-1)){
    
    for(k in 1:n_patch){
      
      # GROWTH & COMPETITION - APHID
      
      # change in A due to growth
      dA1_growth_comp = A1[t,k]*exp((r1-alpha11*A1[t,k]-alpha12*A2[t,k]-alphai3*A3[t,k])*dt) - A1[t,k]
      dA2_growth_comp = A2[t,k]*exp((r2-alpha22*A2[t,k]-alpha21*A1[t,k]-alphai3*A3[t,k])*dt) - A2[t,k]
      dA3_growth_comp = A3[t,k]*exp((r3-alpha33*A3[t,k]-alpha3j*A1[t,k]-alpha3j*A2[t,k])*dt) - A3[t,k]
      
      # PREDATION
      
      # potential number of parasitized aphids
      if(t*dt<tau){
        # type I response
        A1_parasit = beta1*A1[t,k]*P[t,k] * dt
        A2_parasit = beta2*A2[t,k]*P[t,k] * dt
        A3_parasit = beta3*A3[t,k]*P[t,k] * dt
        
      } else {
        # type I response
        A1_parasit = beta1*A1[t,k]*f*P[t,k] * dt
        A2_parasit = beta2*A2[t,k]*f*P[t,k] * dt
        A3_parasit = beta3*A3[t,k]*f*P[t,k] * dt
      }
      # check if number of parasitized aphids > total number of aphids
      if(A1_parasit > A1[t,k]){A1_parasit = A1[t,k]}
      if(A2_parasit > A2[t,k]){A2_parasit = A2[t,k]}
      if(A3_parasit > A3[t,k]){A3_parasit = A3[t,k]}
      
      # total number of aphids parasitized
      n_parasit = A1_parasit + A2_parasit + A3_parasit
      
      # potential number of parasitized aphids
      parasit_max = pmax_dt*P[t,k]
      
      # check if maximum number of parasitized aphids exceeded
      if(n_parasit > parasit_max){
        # total number of parasitized aphids
        A_parasit[t,k] = parasit_max
        # change in A due to parasitism
        dA1_parasit = parasit_max * A1_parasit / (A1_parasit+A2_parasit+A3_parasit)
        dA2_parasit = parasit_max * A2_parasit / (A1_parasit+A2_parasit+A3_parasit)
        dA3_parasit = parasit_max * A3_parasit / (A1_parasit+A2_parasit+A3_parasit)
      } else {
        # total number of parasitized aphids
        A_parasit[t,k] = n_parasit
        # change in A due to parasitism
        dA1_parasit = A1_parasit
        dA2_parasit = A2_parasit
        dA3_parasit = A3_parasit
      }
      
      # BIRTHS - PARASITOID
      
      # dP due to births
      if(t*dt>=tau){
        dP_birth = A_parasit[t-(tau-1)/dt,k]
      } else {
        dP_birth = 0
      }
      
      # DEATHS - PARASITOID
      
      # total number of ptoids across all patches
      P_total = sum(P[t,])
      
      # dP due to deaths
      if(t*dt==(lambda-1)){
        dP_death = P[t,k]
      } else if(t*dt>(tau+lambda-1) & P_total>0) {
        # total ptoids emerged at t-(lambda-1) = total deaths across all patches
        P_death_total = sum(A_parasit[t-(lambda-1)/dt,])
        # number of deaths in current patch (proportional to number of ptoids)
        dP_death = P[t,k]/P_total * P_death_total
      } else {
        dP_death = 0
      }
      
      # EMIGRATION
      
      # delta for emigration
      delta1 = ifelse(A1[t,k]<Ae1, 0, 1)
      delta2 = ifelse(A2[t,k]<Ae2, 0, 1)
      delta3 = ifelse(A3[t,k]<Ae3, 0, 1)
      deltaP = ifelse(P[t,k]<Pe, 0, 1)
      
      # change in A due to emigration
      dA1_emig = delta1*e1*(A1[t,k]-Ae1) * dt
      dA2_emig = delta2*e2*(A2[t,k]-Ae2) * dt
      dA3_emig = delta3*e3*(A3[t,k]-Ae3) * dt
      dP_emig = deltaP*eP*(P[t,k]-Pe) * dt
      
      # IMIGRATION
      
      # neighbours of patch k
      k_neigh = which(M_adj[k,]==1)
      
      # initialise dA due to immigration
      dA1_immig = 0
      dA2_immig = 0
      dA3_immig = 0
      dP_immig = 0
      
      # loop through neighbouring patches
      for(n in k_neigh){
        
        # delta of neighbouring patch
        delta1_neigh = ifelse(A1[t,n]<Ae1, 0, 1)
        delta2_neigh = ifelse(A2[t,n]<Ae2, 0, 1)
        delta3_neigh = ifelse(A3[t,n]<Ae3, 0, 1)
        deltaP_neigh = ifelse(P[t,n]<Pe, 0, 1)
        
        # update dA due to immigration
        dA1_immig = dA1_immig + delta1_neigh*e1*(A1[t,n]-Ae1)/c[n] * dt
        dA2_immig = dA2_immig + delta2_neigh*e2*(A2[t,n]-Ae2)/c[n] * dt
        dA3_immig = dA3_immig + delta3_neigh*e3*(A3[t,n]-Ae3)/c[n] * dt
        dP_immig = dP_immig + deltaP_neigh*eP*(P[t,n]-Pe)/c[n] * dt
      }
      
      # new population size
      A1[t+1,k] = A1[t,k] + dA1_growth_comp - dA1_parasit - dA1_emig + dA1_immig
      A2[t+1,k] = A2[t,k] + dA2_growth_comp - dA2_parasit - dA2_emig + dA2_immig
      A3[t+1,k] = A3[t,k] + dA3_growth_comp - dA3_parasit - dA3_emig + dA3_immig
      P[t+1,k] = P[t,k] + dP_birth - dP_death - dP_emig + dP_immig
      
      if(is.nan(A1[t+1,k]) || is.nan(A2[t+1,k]) || is.nan(A3[t+1,k]) || is.nan(P[t+1,k])){
        print("Error: A=NaN")
        break}
      if(A1[t+1,k]<0){A1[t+1,k]=0}
      if(A2[t+1,k]<0){A2[t+1,k]=0}
      if(A3[t+1,k]<0){A3[t+1,k]=0}
      if(P[t+1,k]<0){P[t+1,k]=0}
      
      # store community and dispersal contributions
      dA1_comm[t+1,k] = dA1_growth_comp - dA1_parasit
      dA1_disp[t+1,k] = - dA1_emig + dA1_immig
      dA2_comm[t+1,k] = dA2_growth_comp - dA2_parasit
      dA2_disp[t+1,k] = - dA2_emig + dA2_immig
      dA3_comm[t+1,k] = dA3_growth_comp - dA3_parasit
      dA3_disp[t+1,k] = - dA3_emig + dA3_immig
      dP_comm[t+1,k] = dP_birth - dP_death
      dP_disp[t+1,k] = - dP_emig + dP_immig
      
    } # k loop
    
  } # t loop
  
  # output dataframe
  df_out = rbind(data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="BB",
                            population_size=as.vector(A1),
                            dN_comm=as.vector(dA1_comm),
                            dN_disp=as.vector(dA1_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="LE",
                            population_size=as.vector(A2),
                            dN_comm=as.vector(dA2_comm),
                            dN_disp=as.vector(dA2_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="A3",
                            population_size=as.vector(A3),
                            dN_comm=as.vector(dA3_comm),
                            dN_disp=as.vector(dA3_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="DR",
                            population_size=as.vector(P),
                            dN_comm=as.vector(dP_comm),
                            dN_disp=as.vector(dP_disp)))
  
  return(df_out)
}

# two aphid species, parasitoid & hyperparasitoid
dt_two_aphid_ptoid_hyper = function(r1, r2, alpha11, alpha22, alpha12, alpha21,
                                    Ae1, Ae2, e1, e2, A10, A20,
                                    beta1, beta2, pmax, f, tau, lambda, Pe, eP,
                                    dt, tmax, patch_state, M_adj){
  
  # time vector
  time = seq(0,tmax,dt)
  
  # number of time steps
  n_dt = length(time)
  
  # nummber of patches
  n_patch = nrow(M_adj)
  
  # number of connections of patches
  c = rowSums(M_adj)
  
  # initialise aphid density matrix
  A1 = matrix(0,n_dt,n_patch)
  A1[1,patch_state==1] = A10
  A2 = matrix(0,n_dt,n_patch)
  A2[1,patch_state==1] = A20
  
  # initialise parasitoid density matrix
  P = matrix(0,n_dt,n_patch)
  P[1:(lambda-1),patch_state==1] = P0
  
  # initialise hyper-parasitoid density matrix
  H = matrix(0,n_dt,n_patch)
  H[1:(lambda-1),patch_state==1] = P0
  
  # initialise matrix of community and dispersal contributions
  dA1_comm = matrix(0,n_dt,n_patch)
  dA1_disp = matrix(0,n_dt,n_patch)
  dA2_comm = matrix(0,n_dt,n_patch)
  dA2_disp = matrix(0,n_dt,n_patch)
  dP_comm = matrix(0,n_dt,n_patch)
  dP_disp = matrix(0,n_dt,n_patch)
  dH_comm = matrix(0,n_dt,n_patch)
  dH_disp = matrix(0,n_dt,n_patch)
  
  # initialise number of parasitized aphids at each timestep
  A_parasit = matrix(0,n_dt,n_patch)
  
  # initialise number of parasitized parasitoids at each timestep
  P_parasit = matrix(0,n_dt,n_patch)
  
  # max number of parasitized aphids per timestep
  pmax_dt = pmax / ((lambda-1)/dt)
  
  for(t in 1:(n_dt-1)){
    
    for(k in 1:n_patch){
      
      # GROWTH & COMPETITION - APHID
      
      # change in A due to growth
      dA1_growth_comp = A1[t,k]*exp((r1-alpha11*A1[t,k]-alpha12*A2[t,k])*dt) - A1[t,k]
      dA2_growth_comp = A2[t,k]*exp((r2-alpha22*A2[t,k]-alpha21*A1[t,k])*dt) - A2[t,k]
      
      # PREDATION - ON APHID
      
      # potential number of parasitized aphids
      if(t*dt<tau){
        # type I response
        A1_parasit = beta1*A1[t,k]*P[t,k] * dt
        A2_parasit = beta2*A2[t,k]*P[t,k] * dt
        
      } else {
        # type I response
        A1_parasit = beta1*A1[t,k]*f*P[t,k] * dt
        A2_parasit = beta2*A2[t,k]*f*P[t,k] * dt
      }
      # check if number of parasitized aphids > total number of aphids
      if(A1_parasit > A1[t,k]){A1_parasit = A1[t,k]}
      if(A2_parasit > A2[t,k]){A2_parasit = A2[t,k]}
      
      # total number of aphids parasitized
      n_parasit = A1_parasit + A2_parasit
      
      # potential number of parasitized aphids
      parasit_max = pmax_dt*P[t,k]
      
      # check if maximum number of parasitized aphids exceeded
      if(n_parasit > parasit_max){
        # total number of parasitized aphids
        A_parasit[t,k] = parasit_max
        # change in A due to parasitism
        dA1_parasit = parasit_max * A1_parasit / (A1_parasit+A2_parasit)
        dA2_parasit = parasit_max * A2_parasit / (A1_parasit+A2_parasit)
      } else {
        # total number of parasitized aphids
        A_parasit[t,k] = n_parasit
        # change in A due to parasitism
        dA1_parasit = A1_parasit
        dA2_parasit = A2_parasit
      }
      
      # PREDATION - ON PARASITOID
      
      # potential number of parasitized parasitoids
      if(t*dt<tau){
        # type I response
        dP_parasit = (beta1+beta2)/2*A_parasit[t,k]*H[t,k] * dt
        
      } else {
        # type I response
        dP_parasit = (beta1+beta2)/2*A_parasit[t,k]*f*H[t,k] * dt
      }
      # check if number of parasitized parasitoids > total number of mummies
      if(dP_parasit > A_parasit[t,k]){dP_parasit = A_parasit[t,k]}
      
      # potential number of parasitized parasitoids
      parasit_max = pmax_dt*H[t,k]
      
      # check if maximum number of parasitized parasitoids exceeded
      if(dP_parasit > parasit_max){dP_parasit = parasit_max}
      
      # update number of parasitized parasitoids
      P_parasit[t,k] = dP_parasit
      
      # BIRTHS - PARASITOID
      
      # dP due to births
      if(t*dt>=tau){
        dP_birth = A_parasit[t-(tau-1)/dt,k]
      } else {
        dP_birth = 0
      }
      
      # BIRTHS - HYPER PARASITOID
      
      # dH due to births
      if(t*dt>=tau){
        dH_birth = P_parasit[t-(tau-1)/dt,k]
      } else {
        dH_birth = 0
      }
      
      # DEATHS - PARASITOID
      
      # total number of ptoids across all patches
      P_total = sum(P[t,])
      
      # dP due to deaths
      if(t*dt==(lambda-1)){
        dP_death = P[t,k]
      } else if(t*dt>(tau+lambda-1) & P_total>0) {
        # total ptoids emerged at t-(lambda-1) = total deaths across all patches
        P_death_total = sum(A_parasit[t-(lambda-1)/dt,])
        # number of deaths in current patch (proportional to number of ptoids)
        dP_death = P[t,k]/P_total * P_death_total
      } else {
        dP_death = 0
      }
      
      # DEATHS - HYPER PARASITOID
      
      # total number of ptoids across all patches
      H_total = sum(H[t,])
      
      # dP due to deaths
      if(t*dt==(lambda-1)){
        dH_death = H[t,k]
      } else if(t*dt>(tau+lambda-1) & H_total>0) {
        # total ptoids emerged at t-(lambda-1) = total deaths across all patches
        H_death_total = sum(P_parasit[t-(lambda-1)/dt,])
        # number of deaths in current patch (proportional to number of ptoids)
        dH_death = H[t,k]/H_total * H_death_total
      } else {
        dH_death = 0
      }
      
      # EMIGRATION
      
      # delta for emigration
      delta1 = ifelse(A1[t,k]<Ae1, 0, 1)
      delta2 = ifelse(A2[t,k]<Ae2, 0, 1)
      deltaP = ifelse(P[t,k]<Pe, 0, 1)
      deltaH = ifelse(H[t,k]<Pe, 0, 1)
      
      # change in A due to emigration
      dA1_emig = delta1*e1*(A1[t,k]-Ae1) * dt
      dA2_emig = delta2*e2*(A2[t,k]-Ae2) * dt
      dP_emig = deltaP*eP*(P[t,k]-Pe) * dt
      dH_emig = deltaP*eP*(H[t,k]-Pe) * dt
      
      # IMIGRATION
      
      # neighbours of patch k
      k_neigh = which(M_adj[k,]==1)
      
      # initialise dA due to immigration
      dA1_immig = 0
      dA2_immig = 0
      dP_immig = 0
      dH_immig = 0
      
      # loop through neighbouring patches
      for(n in k_neigh){
        
        # delta of neighbouring patch
        delta1_neigh = ifelse(A1[t,n]<Ae1, 0, 1)
        delta2_neigh = ifelse(A2[t,n]<Ae2, 0, 1)
        deltaP_neigh = ifelse(P[t,n]<Pe, 0, 1)
        deltaH_neigh = ifelse(H[t,n]<Pe, 0, 1)
        
        # update dA due to immigration
        dA1_immig = dA1_immig + delta1_neigh*e1*(A1[t,n]-Ae1)/c[n] * dt
        dA2_immig = dA2_immig + delta2_neigh*e2*(A2[t,n]-Ae2)/c[n] * dt
        dP_immig = dP_immig + deltaP_neigh*eP*(P[t,n]-Pe)/c[n] * dt
        dH_immig = dH_immig + deltaH_neigh*eP*(H[t,n]-Pe)/c[n] * dt
      }
      
      # new population size
      A1[t+1,k] = A1[t,k] + dA1_growth_comp - dA1_parasit - dA1_emig + dA1_immig
      A2[t+1,k] = A2[t,k] + dA2_growth_comp - dA2_parasit - dA2_emig + dA2_immig
      P[t+1,k] = P[t,k] + dP_birth - dP_death - dP_parasit - dP_emig + dP_immig
      H[t+1,k] = H[t,k] + dH_birth - dH_death - dH_emig + dH_immig
      
      if(is.nan(A1[t+1,k]) || is.nan(A2[t+1,k]) || is.nan(P[t+1,k]) || is.nan(H[t+1,k])){
        print("Error: A=NaN")
        break}
      if(A1[t+1,k]<0){A1[t+1,k]=0}
      if(A2[t+1,k]<0){A2[t+1,k]=0}
      if(P[t+1,k]<0){P[t+1,k]=0}
      if(H[t+1,k]<0){H[t+1,k]=0}
      
      # store community and dispersal contributions
      dA1_comm[t+1,k] = dA1_growth_comp - dA1_parasit
      dA1_disp[t+1,k] = - dA1_emig + dA1_immig
      dA2_comm[t+1,k] = dA2_growth_comp - dA2_parasit
      dA2_disp[t+1,k] = - dA2_emig + dA2_immig
      dP_comm[t+1,k] = dP_birth - dP_death - dP_parasit
      dP_disp[t+1,k] = - dP_emig + dP_immig
      dH_comm[t+1,k] = dH_birth - dH_death
      dH_disp[t+1,k] = - dH_emig + dH_immig
      
    } # k loop
    
  } # t loop
  
  # output dataframe
  df_out = rbind(data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="BB",
                            population_size=as.vector(A1),
                            dN_comm=as.vector(dA1_comm),
                            dN_disp=as.vector(dA1_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="LE",
                            population_size=as.vector(A2),
                            dN_comm=as.vector(dA2_comm),
                            dN_disp=as.vector(dA2_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="DR",
                            population_size=as.vector(P),
                            dN_comm=as.vector(dP_comm),
                            dN_disp=as.vector(dP_disp)),
                 data.frame(expand.grid(t=time, patch=1:n_patch),
                            aphid="HP",
                            population_size=as.vector(H),
                            dN_comm=as.vector(dH_comm),
                            dN_disp=as.vector(dH_disp)))
  
  return(df_out)
}


# replica functions 
# (run multiple replicas sampling parameters from confidence intervals)

dt_one_aphid_rep = function(r_CI, alpha_CI, Ae, e_CI,
                            A0, dt, tmax, patch_state, M_adj, n_rep){
  
  # sample parameter values between confidence intervals
  set.seed(1)
  r_vals = runif(n_rep, min=r_CI[1], max=r_CI[2])
  alpha_vals = runif(n_rep, min=alpha_CI[1], max=alpha_CI[2])
  e_vals = runif(n_rep, min=e_CI[1], max=e_CI[2])

  # initialise dataframe for storing results
  df_out = data.frame(replica=integer(),
                      t=double(),
                      patch=integer(),
                      population_size=double())
  
  # loop through replicas
  for(i in 1:n_rep){
    
    # simulate dynamics
    df_dt = dt_one_aphid(r=r_vals[i], alpha=alpha_vals[i], Ae, e=e_vals[i],
                         A0, dt, tmax, patch_state, M_adj) %>%
      mutate(replica=i)
    
    # combine results
    df_out = rbind(df_out, df_dt)
  }
  
  return(df_out)
}

dt_two_aphid_rep = function(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                            Ae1, Ae2, e1_CI, e2_CI, A10, A20, 
                            dt, tmax, patch_state, M_adj, n_rep){
  
  # sample parameter values between confidence intervals
  set.seed(1)
  r1_vals = runif(n_rep, min=r1_CI[1], max=r1_CI[2])
  alpha11_vals = runif(n_rep, min=alpha11_CI[1], max=alpha11_CI[2])
  e1_vals = runif(n_rep, min=e1_CI[1], max=e1_CI[2])
  alpha12_vals = runif(n_rep, min=alpha12_CI[1], max=alpha12_CI[2])

  set.seed(1)
  r2_vals = runif(n_rep, min=r2_CI[1], max=r2_CI[2])
  alpha22_vals = runif(n_rep, min=alpha22_CI[1], max=alpha22_CI[2])
  e2_vals = runif(n_rep, min=e2_CI[1], max=e2_CI[2])
  alpha21_vals = runif(n_rep, min=alpha21_CI[1], max=alpha21_CI[2])

  # initialise dataframe for storing results
  df_out = data.frame(replica=integer(),
                      t=double(),
                      patch=integer(),
                      aphid=character(),
                      population_size=double())
  
  # loop through replicas
  for(i in 1:n_rep){
    
    # simulate dynamics
    df_dt = dt_two_aphid(r1=r1_vals[i], r2=r2_vals[i], 
                         alpha11=alpha11_vals[i], alpha22=alpha22_vals[i], 
                         alpha12=alpha12_vals[i], alpha21=alpha21_vals[i], 
                         Ae1, Ae2, e1=e1_vals[i], e2=e2_vals[i],
                         A10, A20, dt, tmax, patch_state, M_adj) %>%
      mutate(replica=i)
    
    # combine results
    df_out = rbind(df_out, df_dt)
  }
  
  return(df_out)
}

dt_two_aphid_ptoid_rep = function(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                  Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                  beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                  dt, tmax, patch_state, M_adj, n_rep){
  
  # sample parameter values between confidence intervals
  set.seed(1)
  r1_vals = runif(n_rep, min=r1_CI[1], max=r1_CI[2])
  alpha11_vals = runif(n_rep, min=alpha11_CI[1], max=alpha11_CI[2])
  e1_vals = runif(n_rep, min=e1_CI[1], max=e1_CI[2])
  alpha12_vals = runif(n_rep, min=alpha12_CI[1], max=alpha12_CI[2])
  beta1_vals = runif(n_rep, min=beta1_CI[1], max=beta1_CI[2])
  
  set.seed(1)
  r2_vals = runif(n_rep, min=r2_CI[1], max=r2_CI[2])
  alpha22_vals = runif(n_rep, min=alpha22_CI[1], max=alpha22_CI[2])
  e2_vals = runif(n_rep, min=e2_CI[1], max=e2_CI[2])
  alpha21_vals = runif(n_rep, min=alpha21_CI[1], max=alpha21_CI[2])
  beta2_vals = runif(n_rep, min=beta2_CI[1], max=beta2_CI[2])
  
  # initialise dataframe for storing results
  df_out = data.frame(replica=integer(),
                      t=double(),
                      patch=integer(),
                      aphid=character(),
                      population_size=double())
  
  # loop through replicas
  for(i in 1:n_rep){
    
    # simulate dynamics
    df_dt = dt_two_aphid_ptoid(r1=r1_vals[i], r2=r2_vals[i], 
                               alpha11=alpha11_vals[i], alpha22=alpha22_vals[i], 
                               alpha12=alpha12_vals[i], alpha21=alpha21_vals[i], 
                               Ae1, Ae2, e1=e1_vals[i], e2=e2_vals[i], A10, A20, 
                               beta1=beta1_vals[i], beta2=beta2_vals[i], 
                               pmax, f, tau, lambda, Pe, eP,
                               dt, tmax, patch_state, M_adj) %>%
      mutate(replica=i)
    
    # combine results
    df_out = rbind(df_out, df_dt)
  }
  
  return(df_out)
}

dt_three_aphid_ptoid_rep = function(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                    Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                    beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                    dt, tmax, patch_state, M_adj, n_rep){
  
  # sample parameter values between confidence intervals
  set.seed(1)
  r1_vals = runif(n_rep, min=r1_CI[1], max=r1_CI[2])
  alpha11_vals = runif(n_rep, min=alpha11_CI[1], max=alpha11_CI[2])
  e1_vals = runif(n_rep, min=e1_CI[1], max=e1_CI[2])
  alpha12_vals = runif(n_rep, min=alpha12_CI[1], max=alpha12_CI[2])
  beta1_vals = runif(n_rep, min=beta1_CI[1], max=beta1_CI[2])
  
  set.seed(1)
  r2_vals = runif(n_rep, min=r2_CI[1], max=r2_CI[2])
  alpha22_vals = runif(n_rep, min=alpha22_CI[1], max=alpha22_CI[2])
  e2_vals = runif(n_rep, min=e2_CI[1], max=e2_CI[2])
  alpha21_vals = runif(n_rep, min=alpha21_CI[1], max=alpha21_CI[2])
  beta2_vals = runif(n_rep, min=beta2_CI[1], max=beta2_CI[2])
  
  # initialise dataframe for storing results
  df_out = data.frame(replica=integer(),
                      t=double(),
                      patch=integer(),
                      aphid=character(),
                      population_size=double())
  
  # loop through replicas
  for(i in 1:n_rep){
    
    # simulate dynamics
    df_dt = dt_three_aphid_ptoid(r1=r1_vals[i], r2=r2_vals[i], 
                                 alpha11=alpha11_vals[i], alpha22=alpha22_vals[i], 
                                 alpha12=alpha12_vals[i], alpha21=alpha21_vals[i], 
                                 Ae1, Ae2, e1=e1_vals[i], e2=e2_vals[i], A10, A20, 
                                 beta1=beta1_vals[i], beta2=beta2_vals[i], 
                                 pmax, f, tau, lambda, Pe, eP,
                                 dt, tmax, patch_state, M_adj) %>%
      mutate(replica=i)
    
    # combine results
    df_out = rbind(df_out, df_dt)
  }
  
  return(df_out)
}

dt_two_aphid_ptoid_hyper_rep = function(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                        Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                        beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                        dt, tmax, patch_state, M_adj, n_rep){
  
  # sample parameter values between confidence intervals
  set.seed(1)
  r1_vals = runif(n_rep, min=r1_CI[1], max=r1_CI[2])
  alpha11_vals = runif(n_rep, min=alpha11_CI[1], max=alpha11_CI[2])
  e1_vals = runif(n_rep, min=e1_CI[1], max=e1_CI[2])
  alpha12_vals = runif(n_rep, min=alpha12_CI[1], max=alpha12_CI[2])
  beta1_vals = runif(n_rep, min=beta1_CI[1], max=beta1_CI[2])
  
  set.seed(1)
  r2_vals = runif(n_rep, min=r2_CI[1], max=r2_CI[2])
  alpha22_vals = runif(n_rep, min=alpha22_CI[1], max=alpha22_CI[2])
  e2_vals = runif(n_rep, min=e2_CI[1], max=e2_CI[2])
  alpha21_vals = runif(n_rep, min=alpha21_CI[1], max=alpha21_CI[2])
  beta2_vals = runif(n_rep, min=beta2_CI[1], max=beta2_CI[2])
  
  # initialise dataframe for storing results
  df_out = data.frame(replica=integer(),
                      t=double(),
                      patch=integer(),
                      aphid=character(),
                      population_size=double())
  
  # loop through replicas
  for(i in 1:n_rep){
    
    # simulate dynamics
    df_dt = dt_two_aphid_ptoid_hyper(r1=r1_vals[i], r2=r2_vals[i],
                                     alpha11=alpha11_vals[i], alpha22=alpha22_vals[i], 
                                     alpha12=alpha12_vals[i], alpha21=alpha21_vals[i], 
                                     Ae1, Ae2, e1=e1_vals[i], e2=e2_vals[i], A10, A20, 
                                     beta1=beta1_vals[i], beta2=beta2_vals[i], 
                                     pmax, f, tau, lambda, Pe, eP,
                                     dt, tmax, patch_state, M_adj) %>%
      mutate(replica=i)
    
    # combine results
    df_out = rbind(df_out, df_dt)
  }
  
  return(df_out)
}

# plotting functions
# (plot recovery credit and anova predictions)

plot_model_prediction_partfig <- function(data_test, fig_col){
  
  # linear anova, recovery log(x+1) transformed (to deal with 0-values)
  anova_m2 <- lm(log1p(recovery) ~ landscape_patches*landscape_type*community, data=data_test)
  
  # effect of number of patches
  pred_landscape_patches <- avg_predictions(anova_m2, variables="landscape_patches", by="landscape_patches") %>%
    mutate(sig=ifelse(anova(anova_m2)[1,5]<0.05,"sig","NS"))
  pred_landscape_patches$sig <- factor(pred_landscape_patches$sig, levels=c("NS", "sig"))
  p1 <- ggplot(data=NULL, aes(x=landscape_patches)) +
    geom_point(data=data_test, aes(y=log1p(recovery), col=recovery), 
               position=position_jitter(width=0.1, height=0), alpha=1) +
    geom_errorbar(data=pred_landscape_patches, aes(ymin=conf.low, ymax=conf.high), width=0) +
    geom_point(data=pred_landscape_patches, aes(y=estimate, shape=sig), size=4, fill="white") +
    coord_cartesian(ylim=c(log1p(min(data_test$recovery)),log1p(max(data_test$recovery)))) +
    scale_color_gradient(low="lightgrey", high=fig_col) +
    scale_shape_manual(values=c(21,19), drop=FALSE, guide="none") +
    labs(subtitle="number of communities", y="ln(recovery credit +1)", col="recovery\ncredit") +
    theme(panel.background=element_rect(fill="white", colour="grey"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), 
          axis.title=element_text(size=8), axis.title.x=element_blank(),
          legend.text=element_text(size=8), legend.title=element_text(size=8),
          plot.subtitle=element_text(size=10))
  
  # effect of landscape type
  pred_landscape_type <- avg_predictions(anova_m2, variables="landscape_type", by="landscape_type") %>%
    mutate(sig=ifelse(anova(anova_m2)[2,5]<0.05,"sig","NS"))
  pred_landscape_type$sig <- factor(pred_landscape_type$sig, levels=c("NS", "sig"))
  p2 <- ggplot(data=NULL, aes(x=landscape_type)) +
    geom_point(data=data_test, aes(y=log1p(recovery), col=recovery), 
               position=position_jitter(width=0.1, height=0), alpha=1) +
    geom_errorbar(data=pred_landscape_type, aes(ymin=conf.low, ymax=conf.high), width=0) +
    geom_point(data=pred_landscape_type, aes(y=estimate, shape=sig), size=4, fill="white") +
    coord_cartesian(ylim=c(log1p(min(data_test$recovery)),log1p(max(data_test$recovery)))) +
    scale_color_gradient(low="lightgrey", high=fig_col) +
    scale_shape_manual(values=c(21,19), drop=FALSE, guide="none") +
    labs(subtitle="location of communities", y="ln(recovery credit +1)", col="recovery\ncredit") +
    theme(panel.background=element_rect(fill="white", colour="grey"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), 
          axis.title=element_blank(),
          legend.text=element_text(size=8), legend.title=element_text(size=8),
          plot.subtitle=element_text(size=10))
  
  # effect of community type
  pred_community <- avg_predictions(anova_m2, variables="community", by="community") %>%
    mutate(sig=ifelse(anova(anova_m2)[3,5]<0.05,"sig","NS"))
  pred_community$sig <- factor(pred_community$sig, levels=c("NS", "sig"))
  p3 <- ggplot(data=NULL, aes(x=community)) +
    geom_point(data=data_test, aes(y=log1p(recovery), col=recovery), 
               position=position_jitter(width=0.1, height=0), alpha=1) +
    geom_errorbar(data=pred_community, aes(ymin=conf.low, ymax=conf.high), width=0) +
    geom_point(data=pred_community, aes(y=estimate), size=4) +
    coord_cartesian(ylim=c(log1p(min(data_test$recovery)),log1p(max(data_test$recovery)))) +
    scale_color_gradient(low="lightgrey", high=fig_col) +
    labs(subtitle="community food web", y="ln(recovery credit +1)", col="recovery\ncredit") +
    theme(panel.background=element_rect(fill="white", colour="grey"),
          panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          axis.text.x=element_text(size=10), axis.text.y=element_text(size=8), 
          axis.title=element_blank(),
          legend.text=element_text(size=8), legend.title=element_text(size=8),
          plot.subtitle=element_text(size=10))
  
  # combined plot
  plot_out = ggarrange(p1, p2, p3, nrow=1, common.legend=TRUE, legend="right")
  
  return(plot_out)
}

plot_model_prediction_fullfig <- function(data_pop, data_metapop, aphid_plot, land_plot, comms){
  
  # EMPTY PATCHES
  
  # average across equivalent patches
  data_test <- data_pop %>%
    filter(patch_type=="empty", aphid==aphid_plot, land==land_plot, community %in% comms) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    droplevels()
  
  # plot
  plot_empty = plot_model_prediction_partfig(data_test, fig_col="black")
  
  
  # POPULATED PATCHES
  
  # average across equivalent patches
  data_test <- data_pop %>%
    filter(patch_type=="populated", aphid==aphid_plot, land==land_plot, community %in% comms) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    droplevels()
  
  # plot
  plot_populated = plot_model_prediction_partfig(data_test, fig_col=col_green)
  
  
  # METAPOPULATION
  
  # subset data for analysis
  data_test <- data_metapop %>%
    filter(aphid==aphid_plot, land==land_plot, community %in% comms) %>%
    droplevels()
  
  # plot
  plot_meta = plot_model_prediction_partfig(data_test, fig_col=col_blue)
  
  # combine plots
  plot_out = ggarrange(plot_empty, plot_populated, plot_meta,
                       nrow=3, labels=c("A","B","C"))
  
  return(plot_out)
}


# SIMULATIONS - experiment ----

# ONE APHID

df_1A_1C = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                            A0=A10, dt, tmax, patch_state=c(1,0,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", aphid="BB", community="BB")
df_1A_1P = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI,
                            A0=A10, dt, tmax, patch_state=c(0,1,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", aphid="BB", community="BB")
df_1A_4C = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI,
                            A0=A10, dt, tmax, patch_state=c(1,1,1,1,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", aphid="BB", community="BB")
df_1A_4P = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                            A0=A10, dt, tmax, patch_state=c(0,1,1,1,1), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", aphid="BB", community="BB")


# TWO APHIDS

df_2A_1C = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                           Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                           dt, tmax, patch_state=c(1,0,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE")
df_2A_1P = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                           Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                           dt, tmax, patch_state=c(0,1,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE")
df_2A_4C = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                           Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                           dt, tmax, patch_state=c(1,1,1,1,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE")
df_2A_4P = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                           Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                           dt, tmax, patch_state=c(0,1,1,1,1), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE")


# TWO APHIDS & PARASITOID

df_2AP_1C = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                 Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                 beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                 dt, tmax, patch_state=c(1,0,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE-DR")
df_2AP_1P = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                 Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                 beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                 dt, tmax, patch_state=c(0,1,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE-DR")
df_2AP_4C = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                 Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                 beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                 dt, tmax, patch_state=c(1,1,1,1,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE-DR")
df_2AP_4P = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                 Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                 beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                 dt, tmax, patch_state=c(0,1,1,1,1), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE-DR")


# COMBINE MODEL RESULTS

# dataframe with population sizes
df_pop = rbind(df_1A_1C, df_1A_1P, df_1A_4C, df_1A_4P, 
               df_2A_1C, df_2A_1P, df_2A_4C, df_2A_4P, 
               df_2AP_1C, df_2AP_1P, df_2AP_4C, df_2AP_4P) %>%
  mutate(patch=paste0("patch",patch),
         landscape_patches=as.factor(landscape_patches),
         landscape_type=as.factor(landscape_type),
         aphid=as.factor(aphid),
         community=as.factor(community),
         landscape=ifelse(landscape_type=="central", paste0(landscape_patches,"C"), paste0(landscape_patches,"P")),
         land="L0")

# dataframe with metapopulation sizes
df_metapop = df_pop %>%
  group_by(replica, t, landscape_patches, landscape_type, aphid, community, landscape, land) %>%
  summarise(metapopulation_size=sum(population_size)) %>%
  ungroup()

# dataframe with population recovery credit
df_RC_pop = df_pop %>%
  group_by(replica, patch, landscape_patches, landscape_type, aphid, community, land) %>%
  mutate(A_diff=(population_size+lead(population_size))/2*dt,
         A_dN_comm=(dN_comm+lead(dN_comm))/2*dt,
         A_dN_disp=(dN_disp+lead(dN_disp))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE),
            A_dN_comm=sum(A_dN_comm,na.rm=TRUE),
            A_dN_disp=sum(A_dN_disp,na.rm=TRUE)) %>%
  left_join(., df_pop %>% filter(t==0) %>% select(-t, -dN_comm, -dN_disp) %>%
              rename(pop_t0=population_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-pop_t0*tmax) %>%
  ungroup() %>%
  mutate(patch_type=ifelse(pop_t0!=0, "populated", "empty"))

# dataframe with metapopulation recovery credit
df_RC_metapop = df_metapop %>%
  group_by(replica, landscape_patches, landscape_type, aphid, community, landscape, land) %>%
  mutate(A_diff=(metapopulation_size+lead(metapopulation_size))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE)) %>%
  left_join(., df_metapop %>% filter(t==0) %>% select(-t) %>%
              rename(metapop_t0=metapopulation_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-metapop_t0*tmax) %>%
  ungroup()


# SIMULATIONS - larger landscapes ----

# L1             1 2 3 4 5 6 7 8 9 10 11 12 13
M_adj_L1 = matrix(c(0,1,1,1,1,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,1,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,1,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,1,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,1, 0, 0, 0, 0,
                    0,1,0,0,0,0,0,0,0, 1, 0, 0, 0,
                    0,0,1,0,0,0,0,0,0, 0, 1, 0, 0,
                    0,0,0,1,0,0,0,0,0, 0, 0, 1, 0,
                    0,0,0,0,1,0,0,0,0, 0, 0, 0, 1,
                    0,0,0,0,0,1,0,0,0, 0, 0, 0, 0,
                    0,0,0,0,0,0,1,0,0, 0, 0, 0, 0,
                    0,0,0,0,0,0,0,1,0, 0, 0, 0, 0,
                    0,0,0,0,0,0,0,0,1, 0, 0, 0, 0), 13,13)

# L2             1 2 3 4 5 6 7 8 9 10 11 12 13
M_adj_L2 = matrix(c(0,1,1,1,1,1,1,1,1, 1, 1, 1, 1,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0,
                    1,0,0,0,0,0,0,0,0, 0, 0, 0, 0), 13,13)

# ONE APHID

df_1A_1C_L1 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                               A0=A10, dt, tmax, patch_state=c(1,0,0,0,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", aphid="BB", community="BB", land="L1")
df_1A_1P_L1 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI,
                               A0=A10, dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,0,0,0), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", aphid="BB", community="BB", land="L1")
df_1A_4C_L1 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI,
                               A0=A10, dt, tmax, patch_state=c(1,1,1,1,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", aphid="BB", community="BB", land="L1")
df_1A_4P_L1 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                               A0=A10, dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,1,1,1), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", aphid="BB", community="BB", land="L1")

df_1A_1C_L2 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                               A0=A10, dt, tmax, patch_state=c(1,0,0,0,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", aphid="BB", community="BB", land="L2")
df_1A_1P_L2 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI,
                               A0=A10, dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,0,0,0), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", aphid="BB", community="BB", land="L2")
df_1A_4C_L2 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI,
                               A0=A10, dt, tmax, patch_state=c(1,1,1,1,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", aphid="BB", community="BB", land="L2")
df_1A_4P_L2 = dt_one_aphid_rep(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                               A0=A10, dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,1,1,1), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", aphid="BB", community="BB", land="L2")


# TWO APHIDS

df_2A_1C_L1 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(1,0,0,0,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE", land="L1")
df_2A_1P_L1 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,0,0,0), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE", land="L1")
df_2A_4C_L1 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(1,1,1,1,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE", land="L1")
df_2A_4P_L1 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,1,1,1), 
                               M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE", land="L1")

df_2A_1C_L2 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(1,0,0,0,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE", land="L2")
df_2A_1P_L2 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,0,0,0), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE", land="L2")
df_2A_4C_L2 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(1,1,1,1,0,0,0,0,0,0,0,0,0), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE", land="L2")
df_2A_4P_L2 = dt_two_aphid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                               Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                               dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,1,1,1), 
                               M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE", land="L2")


# TWO APHIDS & PARASITOID

df_2AP_1C_L1 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(1,0,0,0,0,0,0,0,0,0,0,0,0), 
                                      M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE-DR", land="L1")
df_2AP_1P_L1 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,0,0,0), 
                                      M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE-DR", land="L1")
df_2AP_4C_L1 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(1,1,1,1,0,0,0,0,0,0,0,0,0), 
                                      M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE-DR", land="L1")
df_2AP_4P_L1 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,1,1,1), 
                                      M_adj=M_adj_L1, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE-DR", land="L1")

df_2AP_1C_L2 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(1,0,0,0,0,0,0,0,0,0,0,0,0), 
                                      M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE-DR", land="L2")
df_2AP_1P_L2 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,0,0,0), 
                                      M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE-DR", land="L2")
df_2AP_4C_L2 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(1,1,1,1,0,0,0,0,0,0,0,0,0), 
                                      M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE-DR", land="L2")
df_2AP_4P_L2 = dt_two_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                      Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                      beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                      dt, tmax, patch_state=c(0,0,0,0,0,0,0,0,0,1,1,1,1), 
                                      M_adj=M_adj_L2, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE-DR", land="L2")


# COMBINE MODEL RESULTS

# dataframe with population sizes
df_pop_LAND = rbind(rbind(df_1A_1C, df_1A_1P, df_1A_4C, df_1A_4P, 
                          df_2A_1C, df_2A_1P, df_2A_4C, df_2A_4P, 
                          df_2AP_1C, df_2AP_1P, df_2AP_4C, df_2AP_4P) %>%
                      mutate(land="L0"),
                    rbind(df_1A_1C_L1, df_1A_1P_L1, df_1A_4C_L1, df_1A_4P_L1, 
                          df_2A_1C_L1, df_2A_1P_L1, df_2A_4C_L1, df_2A_4P_L1, 
                          df_2AP_1C_L1, df_2AP_1P_L1, df_2AP_4C_L1, df_2AP_4P_L1,
                          df_1A_1C_L2, df_1A_1P_L2, df_1A_4C_L2, df_1A_4P_L2, 
                          df_2A_1C_L2, df_2A_1P_L2, df_2A_4C_L2, df_2A_4P_L2, 
                          df_2AP_1C_L2, df_2AP_1P_L2, df_2AP_4C_L2, df_2AP_4P_L2)) %>%
  mutate(patch=paste0("patch",patch),
         landscape_patches=as.factor(landscape_patches),
         landscape_type=as.factor(landscape_type),
         aphid=as.factor(aphid),
         community=as.factor(community),
         landscape=ifelse(landscape_type=="central", paste0(landscape_patches,"C"), paste0(landscape_patches,"P")))

# dataframe with metapopulation sizes
df_metapop_LAND = df_pop_LAND %>%
  group_by(replica, t, landscape_patches, landscape_type, aphid, community, landscape, land) %>%
  summarise(metapopulation_size=sum(population_size)) %>%
  ungroup()

# dataframe with population recovery credit
df_RC_pop_LAND = df_pop_LAND %>%
  group_by(replica, patch, landscape_patches, landscape_type, aphid, community, land) %>%
  mutate(A_diff=(population_size+lead(population_size))/2*dt,
         A_dN_comm=(dN_comm+lead(dN_comm))/2*dt,
         A_dN_disp=(dN_disp+lead(dN_disp))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE),
            A_dN_comm=sum(A_dN_comm,na.rm=TRUE),
            A_dN_disp=sum(A_dN_disp,na.rm=TRUE)) %>%
  left_join(., df_pop_LAND %>% filter(t==0) %>% select(-t, -dN_comm, -dN_disp) %>%
              rename(pop_t0=population_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-pop_t0*tmax) %>%
  ungroup() %>%
  mutate(patch_type=ifelse(pop_t0!=0, "populated", "empty"))

# dataframe with metapopulation recovery credit
df_RC_metapop_LAND = df_metapop_LAND %>%
  group_by(replica, landscape_patches, landscape_type, aphid, community, landscape, land) %>%
  mutate(A_diff=(metapopulation_size+lead(metapopulation_size))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE)) %>%
  left_join(., df_metapop_LAND %>% filter(t==0) %>% select(-t) %>%
              rename(metapop_t0=metapopulation_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-metapop_t0*tmax) %>%
  ungroup()


# SIMULATIONS - larger communities ----

# THREE APHIDS & PARASITOID

df_3AP_1C = dt_three_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                     Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                     beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                     dt, tmax, patch_state=c(1,0,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE-A3-DR")
df_3AP_1P = dt_three_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                     Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                     beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                     dt, tmax, patch_state=c(0,1,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE-A3-DR")
df_3AP_4C = dt_three_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                     Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                     beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                     dt, tmax, patch_state=c(1,1,1,1,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE-A3-DR")
df_3AP_4P = dt_three_aphid_ptoid_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                                     Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                     beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                     dt, tmax, patch_state=c(0,1,1,1,1), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE-A3-DR")


# TWO APHIDS & PARASITOID & HYPER-PARASITOID

df_2APH_1C = dt_two_aphid_ptoid_hyper_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                          Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                          beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                          dt, tmax, patch_state=c(1,0,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="central", community="BB-LE-DR-HP")
df_2APH_1P = dt_two_aphid_ptoid_hyper_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                          Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                          beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                          dt, tmax, patch_state=c(0,1,0,0,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=1, landscape_type="peripheral", community="BB-LE-DR-HP")
df_2APH_4C = dt_two_aphid_ptoid_hyper_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                          Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                          beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                          dt, tmax, patch_state=c(1,1,1,1,0), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="central", community="BB-LE-DR-HP")
df_2APH_4P = dt_two_aphid_ptoid_hyper_rep(r1_CI, r2_CI, alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI, 
                                          Ae1, Ae2, e1_CI, e2_CI, A10, A20,
                                          beta1_CI, beta2_CI, pmax, f, tau, lambda, Pe, eP,
                                          dt, tmax, patch_state=c(0,1,1,1,1), M_adj, n_rep) %>% 
  mutate(landscape_patches=4, landscape_type="peripheral", community="BB-LE-DR-HP")


# COMBINE MODEL RESULTS

# dataframe with population sizes
df_pop_COMM = rbind(df_1A_1C, df_1A_1P, df_1A_4C, df_1A_4P, 
                    df_2A_1C, df_2A_1P, df_2A_4C, df_2A_4P,
                    df_2AP_1C, df_2AP_1P, df_2AP_4C, df_2AP_4P,
                    df_3AP_1C, df_3AP_1P, df_3AP_4C, df_3AP_4P,
                    df_2APH_1C, df_2APH_1P, df_2APH_4C, df_2APH_4P) %>%
  mutate(patch=paste0("patch",patch),
         landscape_patches=as.factor(landscape_patches),
         landscape_type=as.factor(landscape_type),
         aphid=as.factor(aphid),
         community=as.factor(community),
         landscape=ifelse(landscape_type=="central", paste0(landscape_patches,"C"), paste0(landscape_patches,"P")),
         land="L0")

# dataframe with metapopulation sizes
df_metapop_COMM = df_pop_COMM %>%
  group_by(replica, t, landscape_patches, landscape_type, aphid, community, landscape, land) %>%
  summarise(metapopulation_size=sum(population_size)) %>%
  ungroup()

# dataframe with population recovrey credit
df_RC_pop_COMM = df_pop_COMM %>%
  group_by(replica, patch, landscape_patches, landscape_type, aphid, community, land) %>%
  mutate(A_diff=(population_size+lead(population_size))/2*dt,
         A_dN_comm=(dN_comm+lead(dN_comm))/2*dt,
         A_dN_disp=(dN_disp+lead(dN_disp))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE),
            A_dN_comm=sum(A_dN_comm,na.rm=TRUE),
            A_dN_disp=sum(A_dN_disp,na.rm=TRUE)) %>%
  left_join(., df_pop_COMM %>% filter(t==0) %>% select(-t, -dN_comm, -dN_disp) %>%
              rename(pop_t0=population_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-pop_t0*tmax) %>%
  ungroup() %>%
  mutate(patch_type=ifelse(pop_t0!=0, "populated", "empty"))

# dataframe with metapopulation recovrey credit
df_RC_metapop_COMM = df_metapop_COMM %>%
  group_by(replica, landscape_patches, landscape_type, aphid, community, landscape, land) %>%
  mutate(A_diff=(metapopulation_size+lead(metapopulation_size))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE)) %>%
  left_join(., df_metapop_COMM %>% filter(t==0) %>% select(-t) %>%
              rename(metapop_t0=metapopulation_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-metapop_t0*tmax) %>%
  ungroup()


# PLOTS - experiment ----

# TIME SERIES

ggplot(data=filter(df_pop, aphid=="BB") %>% na.omit(), 
       aes(x=t, y=population_size, group=interaction(replica,community), col=community)) +
  geom_line(alpha=0.1) +
  facet_grid(landscape~patch, scales="free_y") +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  scale_color_manual(values=palette_community) +
  labs(x="time (day)", y="population size") +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=8), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8))

ggplot(data=filter(df_pop, aphid=="LE") %>% na.omit(), 
       aes(x=t, y=population_size, group=interaction(replica,community), col=community)) +
  geom_line(alpha=0.1) +
  facet_grid(landscape~patch, scales="free_y") +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  scale_color_manual(values=palette_community[c(2,3)]) +
  labs(x="time (day)", y="population size") +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=8), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8))

# contributions to population size
df_pop = df_pop %>%
  left_join(., df_pop %>% filter(t==0) %>% 
              mutate(patch_type=ifelse(population_size==0,"empty","populated")) %>%
              select(landscape, community,aphid, patch, replica, patch_type)) %>%
  mutate(dN_disp_fract=dN_disp/(abs(dN_disp)+abs(dN_comm)),
         dN_disp_fract=ifelse(is.nan(dN_disp_fract),0,dN_disp_fract))

ggplot(data=filter(df_pop, aphid=="BB", community=="BB") %>% na.omit(), 
       aes(x=t, y=dN_disp_fract, group=interaction(patch,replica), col=patch_type)) +
  geom_hline(yintercept=0) +
  geom_line(alpha=0.1) +
  facet_grid(landscape~patch, scales="free_y") +
  scale_color_manual(values=c("black",col_green)) +
  labs(x="time (day)", y="relative contribution of dispersal to population change", col="patch type") +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=8), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")


# RECOVERY CREDIT

# combine dataframes for plotting
data_plot <- rbind(df_RC_pop %>%
                     filter(patch_type=="empty") %>%
                     group_by(community, landscape, replica, aphid, land) %>%
                     summarise(recovery=mean(recovery)) %>%
                     ungroup() %>%
                     mutate(scale="empty patches"),
                   df_RC_pop %>%
                     filter(patch_type=="populated") %>%
                     group_by(community, landscape, replica, aphid, land) %>%
                     summarise(recovery=mean(recovery)) %>%
                     ungroup() %>%
                     mutate(scale="populated patches"),
                   df_RC_metapop %>% select(-c(landscape_patches, landscape_type, A_diff, metapop_t0)) %>%
                     mutate(scale="metapopulation"))
data_plot$scale <-  factor(data_plot$scale, c("empty patches","populated patches","metapopulation"))


ggplot(data=filter(data_plot, aphid=="BB"), 
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(alpha=0) +
  facet_grid(scale~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=8), 
        strip.text=element_text(size=8), 
        legend.position="none")

ggplot(data=filter(data_plot, aphid=="LE"), 
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(alpha=0) +
  facet_grid(scale~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=8), 
        strip.text=element_text(size=8), 
        legend.position="none")


# ANOVA PREDICTIONS

plot_model_prediction_fullfig(data_pop=df_RC_pop, data_metapop=df_RC_metapop, 
                              aphid_plot="BB", land_plot="L0", 
                              comms=c("BB","BB-LE","BB-LE-DR"))

plot_model_prediction_fullfig(data_pop=df_RC_pop, data_metapop=df_RC_metapop, 
                              aphid_plot="LE", land_plot="L0", 
                              comms=c("BB","BB-LE","BB-LE-DR"))


# PLOTS - larger landscapes ----

# RECOVERY CREDIT

# combine dataframes for plotting
data_plot_land <- rbind(df_RC_pop_LAND %>%
                          filter(patch_type=="empty") %>%
                          group_by(community, landscape_patches, landscape_type, landscape, replica, aphid, land) %>%
                          summarise(recovery=mean(recovery)) %>%
                          ungroup() %>%
                          mutate(scale="empty patches"),
                        df_RC_pop_LAND %>%
                          filter(patch_type=="populated") %>%
                          group_by(community, landscape_patches, landscape_type, landscape, replica, aphid, land) %>%
                          summarise(recovery=mean(recovery)) %>%
                          ungroup() %>%
                          mutate(scale="populated patches"),
                        df_RC_metapop_LAND %>% select(-c(A_diff, metapop_t0)) %>%
                          mutate(scale="metapopulation"))
data_plot_land$scale <-  factor(data_plot_land$scale, c("empty patches","populated patches","metapopulation"))

ggplot(data=filter(data_plot_land, aphid=="BB", community=="BB"), 
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(alpha=0) +
  facet_grid(scale~land, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none")

ggplot(data=filter(data_plot_land, aphid=="BB", community=="BB", scale=="empty patches"), 
       aes(x=landscape, y=log1p(recovery), col=recovery)) +
  geom_jitter() +
  geom_boxplot(alpha=0) +
  facet_grid(~land, scales="free_y") +
  scale_color_gradient(low="lightgrey", high="black") +
  labs(y="ln(recovery credit +1)", col="recovery\ncredit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(size=8), axis.text.y=element_text(size=8),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        strip.text=element_blank())


# ANOVA PREDICTIONS

plot_model_prediction_fullfig(data_pop=df_RC_pop_LAND, data_metapop=df_RC_metapop_LAND, 
                              aphid_plot="BB", land_plot="L1", 
                              comms=c("BB","BB-LE","BB-LE-DR"))

plot_model_prediction_fullfig(data_pop=df_RC_pop_LAND, data_metapop=df_RC_metapop_LAND, 
                              aphid_plot="BB", land_plot="L2", 
                              comms=c("BB","BB-LE","BB-LE-DR"))


# PLOTS - larger communities ----

# RECOVERY CREDIT

# combine dataframes for plotting
data_plot_comm <- rbind(df_RC_pop_COMM %>%
                          filter(patch_type=="empty") %>%
                          group_by(community, landscape_patches, landscape_type, landscape, replica, aphid, land) %>%
                          summarise(recovery=mean(recovery)) %>%
                          ungroup() %>%
                          mutate(scale="empty patches"),
                        df_RC_pop_COMM %>%
                          filter(patch_type=="populated") %>%
                          group_by(community, landscape_patches, landscape_type, landscape, replica, aphid, land) %>%
                          summarise(recovery=mean(recovery)) %>%
                          ungroup() %>%
                          mutate(scale="populated patches"),
                        df_RC_metapop_COMM %>% select(-c(A_diff, metapop_t0)) %>%
                          mutate(scale="metapopulation"))
data_plot_comm$scale <-  factor(data_plot_comm$scale, 
                                c("empty patches","populated patches","metapopulation"))
data_plot_comm$community <-  factor(data_plot_comm$community, 
                                    c("BB","BB-LE","BB-LE-DR","BB-LE-A3-DR","BB-LE-DR-HP"))

p1 = ggplot(data=filter(data_plot_comm, aphid=="BB", scale=="empty patches"), 
            aes(x=community, y=log1p(recovery), col=recovery)) +
  geom_jitter() +
  geom_boxplot(alpha=0) +
  facet_wrap(~scale, scales="free_y") +
  scale_color_gradient(low="lightgrey", high="black") +
  labs(y="ln(recovery credit +1)") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(size=8, angle=90, hjust=1), axis.text.y=element_text(size=8),
        strip.text=element_text(size=10), axis.title.x=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")

p2 = ggplot(data=filter(data_plot_comm, aphid=="BB", scale=="populated patches"), 
            aes(x=community, y=log1p(recovery), col=recovery)) +
  geom_jitter() +
  geom_boxplot(alpha=0) +
  facet_wrap(~scale, scales="free_y") +
  scale_color_gradient(low="lightgrey", high=col_green) +
  labs(y="ln(recovery credit +1)") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(size=8, angle=90, hjust=1), axis.text.y=element_text(size=8),
        strip.text=element_text(size=10), axis.title=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")

p3 = ggplot(data=filter(data_plot_comm, aphid=="BB", scale=="metapopulation"), 
            aes(x=community, y=log1p(recovery), col=recovery)) +
  geom_jitter() +
  geom_boxplot(alpha=0) +
  facet_wrap(~scale, scales="free_y") +
  scale_color_gradient(low="lightgrey", high=col_blue, breaks=c(5000,20000,35000)) +
  labs(y="ln(recovery credit +1)") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(size=8, angle=90, hjust=1), axis.text.y=element_text(size=8),
        strip.text=element_text(size=10), axis.title=element_blank(),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")

ggarrange(p1,p2,p3, nrow=1)


# ANOVA PREDICTIONS

plot_model_prediction_fullfig(data_pop=df_RC_pop_COMM, data_metapop=df_RC_metapop_COMM, 
                              aphid_plot="BB", land_plot="L0", 
                              comms=c("BB","BB-LE","BB-LE-DR","BB-LE-A3-DR","BB-LE-DR-HP"))
