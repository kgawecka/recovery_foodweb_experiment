rm(list=ls())

# load packages
library(dplyr)
library(ggplot2)
library(marginaleffects)
library(ggpubr)
library(igraph)
library(data.table)
library(forcats)

# define palettes
palette_landscape = c("#4ea33a80","#357029ff","#307ecd80","#1d4d7cff")
palette_community = c("black","#4da43aff","#2f7eceff")
col_blue = "#2f7eceff"
col_green = "#4da43aff"
col_purple = "#8e1558ff"
col_blue_dark = "#1d4d7cff"
col_green_dark = "#357029ff"
col_orange = "#D55E00"
col_orange_dark = "#994400ff"
col_green_light = "#4ea33a80"

# figure widths
width_singlecol = 82
width_twothirds = 110
width_fullpage = 173


# GENERATE LANDSCAPE ----

# generate scale-free network
n_nodes = 50
g = sample_pa(n=n_nodes, power=1, m=1, directed=FALSE)

# check for isolated nodes
isolated_nodes = V(g)[degree(g) == 0]
cat("Number of isolated nodes:", length(isolated_nodes), "\n")

# plot the network
plot(g, vertex.size = 5, vertex.label = NA, layout = layout_with_fr)

# convert to adjacency matrix
M_adj = as_adjacency_matrix(g, sparse=FALSE)

sum(M_adj)/n_nodes # number of links per node
hist(rowSums(M_adj)) # degree distribution

# write out adjacency matrix
write.table(M_adj, paste0("Output/M_land_",n_nodes,".csv"), row.names=FALSE, col.names=FALSE)


# MODEL PARAMETERS ----

# import dataframe with model parameters
# determined from parametrization experiments
species_parameters = read.csv("Data/model_parameters.csv")

# simulation inputs
n_rep = 100 # number of replicas
dt = 1 # timestep size
tmax = 26 # maxmimum number of timesteps
N0 = c(10,10,10,1,1,1) # initial number of aphid1, ahpid2, aphid3, ptoid1, ptoid2, hyper


# FUNCTIONS ----

# run dynamics
f_dynamics = function(r, alpha, Ae, e, beta, gamma, pmax, f, tau, lambda, Pe, eP,
                      N0, n_aphid, n_ptoid, n_hyper,
                      M_land, patch_state, dt, tmax){
  
  # time vector
  time = seq(0,tmax,dt)
  
  # number of time steps
  n_dt = length(time)
  
  # nummber of patches
  n_patch = nrow(M_land)
  
  # number of connections of patches
  c = rowSums(M_land)
  
  # initialise aphid density array
  A = array(0, c(n_dt,n_patch,n_aphid))
  A[1,patch_state==1,] = N0[1:n_aphid]
  
  # initialise matrix of community and dispersal contributions
  dA_comm = array(0, c(n_dt,n_patch,n_aphid))
  dA_disp = array(0, c(n_dt,n_patch,n_aphid))
  
  # initialise female parasitoid/hyperparasitoid density matrix,
  # number of parasitized aphids/parasitoids at each timestep, amd
  # max number of parasitized aphids per timestep
  if(n_ptoid>0){
    P = array(0, c(n_dt,n_patch,n_ptoid))
    P[1:(lambda-1),patch_state==1,] = N0[4:(3+n_ptoid)]
    A_parasit = array(0, c(n_dt,n_patch,n_ptoid))
    pmax_dt = pmax / ((lambda-1)/dt)
  }
  if(n_hyper==1){
    H = matrix(0,n_dt,n_patch)
    H[1:(lambda-1),patch_state==1] = N0[6]
    P_parasit = matrix(0,n_dt,n_patch)
  }
  
  #P_births = array(0, c(n_dt,n_patch,n_ptoid))
  #P_deaths = array(0, c(n_dt,n_patch,n_ptoid))
  #P_emig = array(0, c(n_dt,n_patch,n_ptoid))
  #P_imig = array(0, c(n_dt,n_patch,n_ptoid))
  #H_births = matrix(0,n_dt,n_patch)
  #H_deaths = matrix(0,n_dt,n_patch)
  #H_emig = matrix(0,n_dt,n_patch)
  #H_imig = matrix(0,n_dt,n_patch)
  
  for(t in 1:(n_dt-1)){
    
    # initialise matrix for post-dispersal densities
    A_disp = matrix(0,n_patch,n_aphid)
    if(n_ptoid>0){
      P_disp = matrix(0,n_patch,n_ptoid)
    }
    if(n_hyper==1){
      H_disp = rep(0,n_patch)
    }
    
    # 1st k loop - dispersal dynamics
    for(k in 1:n_patch){
      
      # EMIGRATION - APHID
      # delta for emigration
      deltaA = ifelse(A[t,k,]<Ae, 0, 1)
      # change in A due to emigration
      dA_emig = deltaA * e * (A[t,k,]-Ae) * dt
      
      # EMIGRATION - PARASITOID
      if(n_ptoid>0){
        deltaP = ifelse(P[t,k,]<Pe, 0, 1)
        dP_emig = deltaP * eP * (P[t,k,]-Pe) * dt
      }
      
      # EMIGRATION - HYPERPARASITOID
      if(n_hyper>0){
        deltaH = ifelse(H[t,k]<Pe, 0, 1)
        dH_emig = deltaH * eP * (H[t,k]-Pe) * dt
      }
      
      # IMIGRATION
      # neighbours of patch k
      k_neigh = which(M_land[k,]==1)
      
      # initialise dA due to immigration
      dA_immig = rep(0,n_aphid)
      if(n_ptoid>0){dP_immig = rep(0,n_ptoid)}
      if(n_hyper>0){dH_immig = 0}
      
      # loop through neighbouring patches
      for(n in k_neigh){
        
        # APHID
        # delta of neighbouring patch
        deltaA_neigh = ifelse(A[t,n,]<Ae, 0, 1)
        # update population change due to immigration
        dA_immig = dA_immig + deltaA_neigh * e * (A[t,n,]-Ae)/c[n] * dt
        
        # PARASITOID
        if(n_ptoid>0){
          deltaP_neigh = ifelse(P[t,n,]<Pe, 0, 1)
          dP_immig = dP_immig + deltaP_neigh * eP * (P[t,n,]-Pe)/c[n] * dt
        }
        
        # HYPERPARASITOID
        if(n_hyper>0){
          deltaH_neigh = ifelse(H[t,n]<Pe, 0, 1)
          dH_immig = dH_immig + deltaH_neigh * eP * (H[t,n]-Pe)/c[n] * dt
        }
      }
      
      # NEW POPULATION SIZE
      
      # APHID
      A_disp[k,] = A[t,k,] - dA_emig + dA_immig
      if(any(is.nan(A_disp[k,]))){
        print("Error: A=NaN")
        break}
      A_disp[k,] = ifelse(A_disp[k,]<0,0,A_disp[k,])
      
      # PARASITOID
      if(n_ptoid>0){
        P_disp[k,] = P[t,k,] - dP_emig + dP_immig
        P_disp[k,] = ifelse(P_disp[k,]<0,0,P_disp[k,])
        
        #P_emig[t,k,] = dP_emig
        #P_imig[t,k,] = dP_immig
      }
      
      # HYPERPARASITOID
      if(n_hyper>0){
        H_disp[k] = H[t,k] - dH_emig + dH_immig
        H_disp[k] = ifelse(H_disp[k]<0,0,H_disp[k])
        
        #H_emig[t,k] = dH_emig
        #H_imig[t,k] = dH_immig
      }
      
      # store dispersal contributions
      dA_disp[t+1,k,] = - dA_emig + dA_immig
      
    } # k loop
    
    # 2nd k look - community dynamics
    for(k in 1:n_patch){
      
      # GROWTH & COMPETITION - APHID
      
      # change in A due to growth
      dA_growth_comp = A_disp[k,] * exp((r - alpha %*% A_disp[k,])*dt) - A_disp[k,]
      
      # PREDATION - PARASITOID ON APHID
      # initialise change in A due to parasitism
      dA_parasit = rep(0,n_aphid)
      if(n_ptoid>0){
        
        # potential number of parasitized aphids of each species (type I response)
        A_parasit_max = A_disp[k,] * beta * P_disp[k,] * dt
        
        # proportional parasitoid preference
        beta_prop = beta / apply(beta,1,sum)
        
        # check if number of parasitized aphids > total number of aphids
        # if true, set number of parasitized aphids to total number of aphids
        # according to parasitoid preference
        A_parasit_exc = rowSums(A_parasit_max)>A_disp[k,]
        A_parasit_max[A_parasit_exc,] = A_disp[k,A_parasit_exc] * beta_prop[A_parasit_exc,]
        
        # number of aphids parasitized by each parasitoid
        if(n_ptoid==1){
          n_parasit_P = sum(A_parasit_max)
        }else{
          n_parasit_P = colSums(A_parasit_max)
        }
        
        # potential number of parasitized aphids by each parasitoid
        parasit_max_P = pmax_dt*P_disp[k,]
        
        for(i in 1:n_ptoid){
          # check if maximum number of parasitized aphids exceeded
          if(n_parasit_P[i] > parasit_max_P[i]){
            # total number of parasitized aphids
            A_parasit[t,k,i] = parasit_max_P[i] 
            # change in A due to parasitism
            dA_parasit = dA_parasit + parasit_max_P[i] * A_parasit_max[,i] / n_parasit_P[i]
          } else {
            # total number of parasitized aphids
            A_parasit[t,k,i] = n_parasit_P[i]
            # change in A due to parasitism
            dA_parasit = dA_parasit + A_parasit_max[,i]
          }
        }
        
      }
      
      # PREDATION - HYPERPARASITOIDS ON PARASITOID 
      if(n_hyper>0){
        
        # potential number of parasitized parasitoids (type I response)
        P_parasit_max = A_parasit[t,k,] * gamma * H_disp[k] * dt
        
        # check if number of parasitized parasitoids > total number of mummies
        # if true, set number of parasitized parasitoids to total number of mummies
        P_parasit_exc = P_parasit_max>A_parasit[t,k,]
        P_parasit_max[P_parasit_exc] = A_parasit[t,k,P_parasit_exc]
        
        # number of parasitoids parasitized
        n_parasit = sum(P_parasit_max)
        
        # potential number of parasitized parasitoids
        parasit_max = pmax_dt*H_disp[k]
        
        # check if maximum number of parasitized parasitoids exceeded
        if(n_parasit > parasit_max){
          # total number of parasitized parasitoids
          P_parasit[t,k] = parasit_max
          # change in P due to parasitism
          dP_parasit = parasit_max_P * P_parasit_max / n_parasit
        } else {
          # total number of parasitized parasitoids
          P_parasit[t,k] = n_parasit
          # change in P due to parasitism
          dP_parasit = P_parasit_max
        }
        
        # update number of parasitized aphids (mummies without parasitoid)
        A_parasit[t,k,] = A_parasit[t,k,] - dP_parasit
      }
      
      # BIRTHS & DEATHS - PARASITOID
      if(n_ptoid>0){
        
        # dP due to births
        if(t*dt>tau){
          dP_birth = f*A_parasit[t-(tau)/dt,k,]
        } else {
          dP_birth = rep(0,n_ptoid)
        }
        
        # total number of parasitoids across all patches
        if(n_ptoid==1){
          P_total = sum(P_disp)
        }else{
          P_total = colSums(P_disp)
        }
        
        # dP due to deaths
        if(t*dt==(lambda-1)){
          dP_death = P_disp[k,]
        } else if(t*dt>(tau+lambda-1)) {
          # total emerged at t-(lambda-1)
          # = total parasitised at t-(tau+lambda-1)
          # = total deaths across all patches
          if(n_ptoid==1){
            P_death_total = sum(f*A_parasit[t-(tau+lambda-1)/dt,,])
          }else{
            P_death_total = colSums(f*A_parasit[t-(tau+lambda-1)/dt,,])
          }
          # number of deaths in current patch (proportional to number of ptoids)
          dP_death = P_disp[k,]/P_total * P_death_total
          dP_death[is.nan(dP_death)] = 0
        } else {
          dP_death = rep(0,n_ptoid)
        }
        
      }
      
      # BIRTHS & DEATHS - HYPERPARASITOID
      if(n_hyper>0){
        
        # dH due to births
        if(t*dt>tau){
          dH_birth = f*P_parasit[t-(tau)/dt,k]
        } else {
          dH_birth = 0
        }
        
        # total number of hyperparasitoids across all patches
        H_total = sum(H_disp)
        
        # dH due to deaths
        if(t*dt==(lambda-1)){
          dH_death = H_disp[k]
        } else if(t*dt>(tau+lambda-1) & H_total>0) {
          # total emerged at t-(lambda-1)
          # = total parasitised at t-(tau+lambda-1)
          # = total deaths across all patches
          H_death_total = sum(f*P_parasit[t-(tau+lambda-1)/dt,])
          # number of deaths in current patch (proportional to number of ptoids)
          dH_death = H_disp[k]/H_total * H_death_total
        } else {
          dH_death = 0
        }
        
      }
      
      # NEW POPULATION SIZE
      
      # APHID
      A[t+1,k,] = A_disp[k,] + dA_growth_comp - dA_parasit
      if(any(is.nan(A[t+1,k,]))){
        print("Error: A=NaN")
        break}
      A[t+1,k,] = ifelse(A[t+1,k,]<0,0,A[t+1,k,])
      
      # PARASITOID
      if(n_ptoid>0){
        P[t+1,k,] = P_disp[k,] + dP_birth - dP_death
        P[t+1,k,] = ifelse(P[t+1,k,]<0,0,P[t+1,k,])
        
        #P_births[t,k,] = dP_birth
        #P_deaths[t,k,] = dP_death
      }
      
      # HYPERPARASITOID
      if(n_hyper>0){
        H[t+1,k] = H_disp[k] + dH_birth - dH_death
        H[t+1,k] = ifelse(H[t+1,k]<0,0,H[t+1,k])
        
        #H_births[t,k] = dH_birth
        #H_deaths[t,k] = dH_death
      }
      
      # store community contributions
      dA_comm[t+1,k,] = dA_growth_comp - dA_parasit
      
    } # k loop
    
  } # t loop
  
  # output dataframe
  df_out = data.frame(expand.grid(t=time, patch=1:n_patch, species=paste0("A",1:n_aphid)),
                      population_size=as.vector(A),
                      dN_comm=as.vector(dA_comm),
                      dN_disp=as.vector(dA_disp))
  if(n_ptoid>0){
    df_out = rbind(df_out,
                   data.frame(expand.grid(t=time, patch=1:n_patch, species=paste0("P",1:n_ptoid)),
                              population_size=as.vector(P),
                              dN_comm=NA,
                              dN_disp=NA))
  }
  if(n_hyper>0){
    df_out = rbind(df_out,
                   data.frame(expand.grid(t=time, patch=1:n_patch),
                              species="H",
                              population_size=as.vector(H),
                              dN_comm=NA,
                              dN_disp=NA))
  }
  
  return(df_out)
}

# run replicate simulations (sampling from confidence intervals)
f_reps = function(species_parameters, N0, n_aphid, n_ptoid, n_hyper, 
                  M_land, patch_state, dt, tmax, n_rep){
  
  # aphid1 - BB - sample species parameter values between confidence intervals
  set.seed(1)
  r1_vals      = runif(n_rep, min=species_parameters[1,3], max=species_parameters[1,4]) # growth rate
  alpha11_vals = runif(n_rep, min=species_parameters[3,3], max=species_parameters[3,4]) # intraspecific competition
  alpha12_vals = runif(n_rep, min=species_parameters[5,3], max=species_parameters[5,4]) # interspecific competition
  e1_vals      = runif(n_rep, min=species_parameters[7,3], max=species_parameters[7,4]) # emigration rate
  Ae1 = as.numeric(species_parameters[species_parameters$parameter=="Ae1",2]) # minimum density for emigration
  
  # aphid2 - LE - sample species parameter values between confidence intervals
  set.seed(1)
  r2_vals      = runif(n_rep, min=species_parameters[2,3], max=species_parameters[2,4])
  alpha22_vals = runif(n_rep, min=species_parameters[4,3], max=species_parameters[4,4])
  alpha21_vals = runif(n_rep, min=species_parameters[6,3], max=species_parameters[6,4])
  e2_vals      = runif(n_rep, min=species_parameters[8,3], max=species_parameters[8,4])
  Ae2 = as.numeric(species_parameters[species_parameters$parameter=="Ae2",2])
  
  # aphid3 - MP - estimate
  r3_vals      = r2_vals # the same as LE
  alpha33_vals = alpha22_vals # the same as LE
  alpha3x_vals = (alpha12_vals + alpha21_vals) / 2 # average of alpha12 & alpha21
  e3_vals      = (e1_vals + e2_vals) / 2 # average of BB & LE
  Ae3          = (Ae1 + Ae2) / 2 # average of BB & LE
  
  # aphid - patasitoid
  beta11_vals  = runif(n_rep, min=species_parameters[11,3], max=species_parameters[11,4]) # parasitisation rate on BB
  beta21_vals   = runif(n_rep, min=species_parameters[12,3], max=species_parameters[12,4]) # parasitisation rate on LE
  beta31_vals   = 2/3*beta21_vals # parasitisation rate on MP lower than LE (ratio between beta21 & beta11 = 2/3)
  beta12_vals  =  1/3*beta21_vals # parasitisation rate on BB low
  beta22_vals   = 1/3*beta21_vals # parasitisation rate on LE low
  beta32_vals   = beta11_vals # parasitisation rate on MP high (the same as DR on BB)
  
  # parasitoid & hyperparasitoid - DR - parameters (AC & hyper assumed the same as DR)
  pmax = as.numeric(species_parameters[species_parameters$parameter=="pmax",2]) # maximum parasitization
  tau = round(as.numeric(species_parameters[species_parameters$parameter=="tau",2]),0) # time to emergence
  f = as.numeric(species_parameters[species_parameters$parameter=="f",2]) # fraction of females
  lambda = round(as.numeric(species_parameters[species_parameters$parameter=="lambda",2]),0) # lifespan 
  eP = as.numeric(species_parameters[species_parameters$parameter=="eP",2]) # emigration rate
  Pe = as.numeric(species_parameters[species_parameters$parameter=="Pe",2]) # minimum density for emigration
  
  # parasitoid - hyperparasitoid
  gamma1 = beta11_vals # parasitisation rate on DR (the same as DR on BB)
  gamma2 = beta21_vals # parasitisation rate on DR (the same as DR on LE)
  
  # initialise dataframe for storing results
  df_out = data.frame(replica=integer(),
                      t=double(),
                      patch=integer(),
                      aphid=character(),
                      population_size=double())
  
  # loop through replicas
  for(i in 1:n_rep){
    
    # species parameters for current replica
    if(n_aphid==1){
      r = r1_vals[i]
      alpha = matrix(alpha11_vals[i],
                     1,1, byrow=TRUE)
      e = e1_vals[i]
      Ae = Ae1
      if(n_ptoid==0){
        beta = NULL
      }else{
        beta = matrix(c(beta11_vals[i], beta12_vals[i]),
                      1,2)[,1:n_ptoid, drop=FALSE]
      }
    }
    if(n_aphid==2){
      r = c(r1_vals[i], r2_vals[i])
      alpha = matrix(c(alpha11_vals[i], alpha12_vals[i],
                       alpha21_vals[i], alpha22_vals[i]),
                     2,2, byrow=TRUE)
      e = c(e1_vals[i], e2_vals[i])
      Ae = c(Ae1, Ae2)
      if(n_ptoid==0){
        beta = NULL
      }else{
        beta = matrix(c(beta11_vals[i], beta12_vals[i],
                        beta21_vals[i], beta22_vals[i]),
                      2,2, byrow=TRUE)[,1:n_ptoid, drop=FALSE]
      }
    }
    if(n_aphid==3){
      r = c(r1_vals[i], r2_vals[i], r3_vals[i])
      alpha = matrix(c(alpha11_vals[i], alpha12_vals[i], alpha3x_vals[i],
                       alpha21_vals[i], alpha22_vals[i], alpha3x_vals[i],
                       alpha3x_vals[i], alpha3x_vals[i], alpha33_vals[i]),
                     3,3, byrow=TRUE)
      e = c(e1_vals[i], e2_vals[i], e3_vals[i])
      Ae = c(Ae1, Ae2, Ae3)
      if(n_ptoid==0){
        beta = NULL
      }else{
        beta = matrix(c(beta11_vals[i], beta12_vals[i],
                        beta21_vals[i], beta22_vals[i],
                        beta31_vals[i], beta32_vals[i]),
                      3,2, byrow=TRUE)[,1:n_ptoid, drop=FALSE]
      }
    }
    if(n_hyper==0){
      gamma = NULL
    }else{
      gamma = c(gamma1,gamma2)[1:n_ptoid]
    }
    
    # simulate dynamics
    df_dt = f_dynamics(r, alpha, Ae, e, beta, gamma, pmax, f, tau, lambda, Pe, eP,
                       N0, n_aphid, n_ptoid, n_hyper,
                       M_land, patch_state, dt, tmax) %>%
      mutate(replica=i)
    
    # combine results
    df_out = rbind(df_out, df_dt)
  }
  
  return(df_out)
}

# plotting functions
plot_model_prediction_partfig4 = function(data_test, data_test_exp, df_comms){
  
  # linear anova, recovery log(x+1) transformed (to deal with 0-values)
  anova_m2 = lm(log1p(recovery) ~ landscape_patches*landscape_type*community, data=data_test)
  anova_m2_exp = lm(log1p(recovery) ~ landscape_patches*landscape_type*community, data=data_test_exp)
  
  # predictions
  preds = rbind(avg_predictions(anova_m2, variables="landscape_patches", by="landscape_patches") %>%
                  mutate(sig=ifelse(anova(anova_m2)[1,5]<0.05,"sig","NS"),
                         effect="number of communities") %>%
                  rename(treatment="landscape_patches"),
                avg_predictions(anova_m2, variables="landscape_type", by="landscape_type") %>%
                  mutate(sig=ifelse(anova(anova_m2)[2,5]<0.05,"sig","NS"),
                         effect="location of communities") %>%
                  rename(treatment="landscape_type"),
                avg_predictions(anova_m2, variables="community", by="community") %>%
                  mutate(sig=ifelse(anova(anova_m2)[3,5]<0.05,"sig","NS"),
                         effect="community food web") %>%
                  rename(treatment="community")) %>%
    left_join(., df_comms, by=c("treatment"="community"))
  preds$sig = factor(preds$sig, levels=c("NS", "sig"))
  preds$effect = factor(preds$effect, 
                        levels=c("number of communities","location of communities","community food web"))
  preds[is.na(preds)] = 0
  
  preds_exp = rbind(avg_predictions(anova_m2_exp, variables="landscape_patches", by="landscape_patches") %>%
                      mutate(sig=ifelse(anova(anova_m2_exp)[1,5]<0.05,"sig","NS"),
                             effect="number of communities") %>%
                      rename(treatment="landscape_patches"),
                    avg_predictions(anova_m2_exp, variables="landscape_type", by="landscape_type") %>%
                      mutate(sig=ifelse(anova(anova_m2_exp)[2,5]<0.05,"sig","NS"),
                             effect="location of communities") %>%
                      rename(treatment="landscape_type"),
                    avg_predictions(anova_m2_exp, variables="community", by="community") %>%
                      mutate(sig=ifelse(anova(anova_m2_exp)[3,5]<0.05,"sig","NS"),
                             effect="community food web") %>%
                      rename(treatment="community")) %>%
    left_join(., df_comms, by=c("treatment"="community"))
  preds_exp$sig = factor(preds_exp$sig, levels=c("NS", "sig"))
  preds_exp$effect = factor(preds_exp$effect, 
                            levels=c("number of communities","location of communities","community food web"))
  preds_exp[is.na(preds_exp)] = 0
  
  # transform data for plotting
  data_plot = rbind(data_test %>% select(landscape_patches,recovery,scale) %>%
                      mutate(effect="number of communities") %>%
                      rename(treatment="landscape_patches"),
                    data_test %>% select(landscape_type,recovery,scale) %>%
                      mutate(effect="location of communities") %>%
                      rename(treatment="landscape_type"),
                    data_test %>% select(community,recovery,scale) %>%
                      mutate(effect="community food web") %>%
                      rename(treatment="community")) %>%
    left_join(., df_comms, by=c("treatment"="community"))
  data_plot$effect = factor(data_plot$effect, 
                            levels=c("number of communities","location of communities","community food web"))
  data_plot[is.na(data_plot)] = 0
  data_plot = data_plot %>%
    mutate(treatment=as.factor(fct_reorder(treatment, n_species+0.01*trophic_levels)))
  
  data_plot_exp = rbind(data_test_exp %>% select(landscape_patches,recovery,scale) %>%
                          mutate(effect="number of communities") %>%
                          rename(treatment="landscape_patches"),
                        data_test_exp %>% select(landscape_type,recovery,scale) %>%
                          mutate(effect="location of communities") %>%
                          rename(treatment="landscape_type"),
                        data_test_exp %>% select(community,recovery,scale) %>%
                          mutate(effect="community food web") %>%
                          rename(treatment="community")) %>%
    left_join(., df_comms, by=c("treatment"="community"))
  data_plot_exp$effect = factor(data_plot_exp$effect, 
                                levels=c("number of communities","location of communities","community food web"))
  data_plot_exp[is.na(data_plot_exp)] = 0
  data_plot_exp = data_plot_exp %>%
    mutate(treatment=as.factor(fct_reorder(treatment, n_species+0.01*trophic_levels)))
  
  groups_to_split = c("1A","1A-1P","1A-1P-1H","2A-1P-1H","3A-1P-1H")
  x_positions = as.numeric(factor(groups_to_split, levels=levels(data_plot$treatment))) + 0.5-4
  vline_data = data.frame(xintercept=x_positions,
                          effect="community food web")
  vline_data$effect = factor(vline_data$effect, 
                             levels=c("number of communities","location of communities","community food web"))
  
  # plot
  
  if(unique(data_test$scale)=="empty patches"){
    p = ggplot(data=NULL, aes(x=treatment)) +
      geom_point(data=data_plot, aes(y=log1p(recovery), col="simulation"), 
                 position=position_jitter(width=0.1, height=0), alpha=0.05, size=0.5) +
      geom_errorbar(data=preds, aes(ymin=conf.low, ymax=conf.high, col="simulation"), 
                    width=0) +
      geom_point(data=preds, aes(y=estimate, shape=sig, col="simulation"), 
                 size=2, fill="white") +
      geom_point(data=data_plot_exp, aes(y=log1p(recovery), col="experiment"), 
                 position=position_jitter(width=0.1, height=0), alpha=0.5, size=0.5) +
      geom_errorbar(data=preds_exp, aes(ymin=conf.low, ymax=conf.high, col="experiment"), 
                    width=0) +
      geom_point(data=preds_exp, aes(y=estimate, shape=sig, col="experiment"),
                 size=2, fill="white") +
      geom_vline(data=vline_data, aes(xintercept=xintercept), col="grey", linetype="dashed") +
      facet_grid(scale~effect, scales="free_x", space="free_x",
                 labeller=labeller(effect=label_wrap_gen(width=12))) +
      coord_cartesian(ylim=c(log1p(min(data_plot$recovery)),log1p(max(data_plot$recovery)))) +
      scale_color_manual(name=NULL, 
                         values=c("simulation"=col_blue_dark, "experiment"=col_orange)) +
      scale_shape_manual(values=c(21,19), drop=FALSE, guide="none") +
      scale_x_discrete(labels=c("1"="1\n \n ","4"="4\n \n ", 
                                "central"="central\n \n ","peripheral"=" \nperipheral\n ", 
                                "1A"=" \n \n1A","2A"=" \n \n2A","1A-1P"=" \n1P\n1A",
                                "3A"=" \n \n3A","1A-2P"=" \n2P\n1A","2A-1P"=" \n1P\n2A",
                                "1A-1P-1H"="1H\n1P\n1A",
                                "2A-2P"=" \n2P\n2A","3A-1P"=" \n1P\n3A",
                                "1A-2P-1H"="1H\n2P\n1A","2A-1P-1H"="1H\n1P\n2A",
                                "3A-2P"=" \n2P\n3A","2A-2P-1H"="1H\n2P\n2A",
                                "3A-1P-1H"="1H\n1P\n3A","3A-2P-1H"="1H\n2P\n3A")) +
      labs(y="ln(recovery credit +1)") +
      theme(panel.background=element_rect(fill="white", colour="grey"),
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=8), axis.text.y=element_text(size=6), 
            axis.title=element_text(size=8), axis.title.x=element_blank(),
            strip.background=element_blank(), 
            strip.text=element_text(size=8),
            legend.position="bottom", legend.text=element_text(size=8), legend.key=element_blank())
  }else{
    p = ggplot(data=NULL, aes(x=treatment)) +
      geom_point(data=data_plot, aes(y=log1p(recovery), col="simulation"), 
                 position=position_jitter(width=0.1, height=0), alpha=0.05, size=0.5) +
      geom_errorbar(data=preds, aes(ymin=conf.low, ymax=conf.high, col="simulation"), 
                    width=0) +
      geom_point(data=preds, aes(y=estimate, shape=sig, col="simulation"), 
                 size=2, fill="white") +
      geom_point(data=data_plot_exp, aes(y=log1p(recovery), col="experiment"), 
                 position=position_jitter(width=0.1, height=0), alpha=0.5, size=0.5) +
      geom_errorbar(data=preds_exp, aes(ymin=conf.low, ymax=conf.high, col="experiment"), 
                    width=0) +
      geom_point(data=preds_exp, aes(y=estimate, shape=sig, col="experiment"),
                 size=2, fill="white") +
      geom_vline(data=vline_data, aes(xintercept=xintercept), col="grey", linetype="dashed") +
      facet_grid(scale~effect, scales="free_x", space="free_x",
                 labeller=labeller(effect=label_wrap_gen(width=12))) +
      coord_cartesian(ylim=c(log1p(min(data_plot$recovery)),log1p(max(data_plot$recovery)))) +
      scale_color_manual(name=NULL, 
                         values=c("simulation"=col_blue_dark, "experiment"=col_orange)) +
      scale_shape_manual(values=c(21,19), drop=FALSE, guide="none") +
      scale_x_discrete(labels=c("1"="1\n \n ","4"="4\n \n ", 
                                "central"="central\n \n ","peripheral"=" \nperipheral\n ", 
                                "1A"=" \n \n1A","2A"=" \n \n2A","1A-1P"=" \n1P\n1A",
                                "3A"=" \n \n3A","1A-2P"=" \n2P\n1A","2A-1P"=" \n1P\n2A",
                                "1A-1P-1H"="1H\n1P\n1A",
                                "2A-2P"=" \n2P\n2A","3A-1P"=" \n1P\n3A",
                                "1A-2P-1H"="1H\n2P\n1A","2A-1P-1H"="1H\n1P\n2A",
                                "3A-2P"=" \n2P\n3A","2A-2P-1H"="1H\n2P\n2A",
                                "3A-1P-1H"="1H\n1P\n3A","3A-2P-1H"="1H\n2P\n3A")) +
      labs(y="ln(recovery credit +1)") +
      theme(panel.background=element_rect(fill="white", colour="grey"),
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            axis.text.x=element_text(size=8), axis.text.y=element_text(size=6), 
            axis.title=element_text(size=8), axis.title.x=element_blank(),
            strip.background=element_blank(), 
            strip.text.x=element_blank(), strip.text.y=element_text(size=8),
            legend.position="bottom", legend.text=element_text(size=8), legend.key=element_blank())
  }
  
  return(p)
}
plot_model_prediction_fullfig4 = function(data_pop, data_metapop, data_pop_exp, data_metapop_exp,
                                          df_comms, species_plot){
  
  aphid_plot = ifelse(species_plot=="A1", "BB", "LE")
  
  # EMPTY PATCHES
  
  # average across equivalent patches
  data_test = data_pop %>%
    filter(patch_type=="empty", species==species_plot) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    mutate(scale="empty patches")
  data_test_exp = data_pop_exp %>%
    filter(patch_type=="empty", aphid==aphid_plot) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    mutate(scale="empty patches")
  
  # plot
  plot_empty = plot_model_prediction_partfig4(data_test, data_test_exp, df_comms)
  
  
  # POPULATED PATCHES
  
  # average across equivalent patches
  data_test = data_pop %>%
    filter(patch_type=="populated", species==species_plot) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    mutate(scale="populated patches")
  data_test_exp = data_pop_exp %>%
    filter(patch_type=="populated", aphid==aphid_plot) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    mutate(scale="populated patches")
  
  # plot
  plot_populated = plot_model_prediction_partfig4(data_test, data_test_exp, df_comms)
  
  
  # METAPOPULATION
  
  # subset data for analysis
  data_test = data_metapop %>%
    filter(species==species_plot) %>%
    mutate(scale="metapopulation")
  data_test_exp = data_metapop_exp %>%
    filter(aphid==aphid_plot) %>%
    mutate(scale="metapopulation")
  
  # plot
  plot_meta = plot_model_prediction_partfig4(data_test, data_test_exp, df_comms)
  
  # combine plots
  plot_out = ggarrange(plot_empty, plot_populated, plot_meta,
                       nrow=3, labels=c("A","B","C"), 
                       common.legend=TRUE, legend="bottom", font.label=list(size=10))
  
  return(plot_out)
}

# SIMULATIONS - larger systems ----

# landscape size
n_patches = 50

# landscape adjacency matrix
M_land = as.matrix(read.table(paste0("Output/M_land_",n_patches,".csv"), quote="\"", comment.char=""))

# plot the network
g = graph_from_adjacency_matrix(M_land, mode="undirected")
layout_fixed = layout_with_fr(g)

plot(g, layout=layout_fixed, 
     vertex.label=V(g),
     #vertex.label=NA,
     vertex.size=10, vertex.color="white", 
     vertex.frame.color=col_green_dark, vertex.frame.width=1,
     edge.color=col_green_dark, edge.width=1,
     vertex.label.color="black", vertex.label.family="Arial", vertex.label.cex=0.5)

# patch properties
df_patch = data.frame(patch=1:nrow(M_land),
                      degree=degree(g),
                      closeness_centrality=closeness(g),
                      betweenness_centrality=betweenness(g),
                      eigen_centrality_values=eigen_centrality(g)$vector)

# patch states (initial communities)
patch_states = matrix(0, 4,nrow(M_land))
patch_states[1,2]              = 1 # 1C (50 patches)
patch_states[2,47]             = 1 # 1P (50 patches)
patch_states[3,c(2,3,4,7)]     = 1 # 4C (50 patches)
patch_states[4,c(47,35,41,37)] = 1 # 4P (50 patches)

# dataframe with communities for simulations
df_comms = data.frame(community=c("1A","1A-1P","1A-2P","1A-1P-1H","1A-2P-1H",
                                  "2A","2A-1P","2A-2P","2A-1P-1H","2A-2P-1H",
                                  "3A","3A-1P","3A-2P","3A-1P-1H","3A-2P-1H"),
                      aphids=c(1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3),
                      ptoids=rep(c(0,1,2,1,2),3),
                      hypers=rep(c(0,0,0,1,1),3)) %>%
  mutate(n_species=aphids+ptoids+hypers,
         trophic_levels=ifelse(hypers==0,ifelse(ptoids==0,1,2),3))

# run simulations
for(i in 1:nrow(df_comms)){
  
  n_aphid = df_comms$aphids[i]
  n_ptoid = df_comms$ptoids[i]
  n_hyper = df_comms$hypers[i]
  
  df_sim = rbind(f_reps(species_parameters, N0, n_aphid, n_ptoid, n_hyper, 
                        M_land, patch_state=patch_states[1,], dt, tmax, n_rep) %>% 
                   mutate(landscape_patches=1, landscape_type="central"),
                 f_reps(species_parameters, N0, n_aphid, n_ptoid, n_hyper, 
                        M_land, patch_state=patch_states[2,], dt, tmax, n_rep) %>% 
                   mutate(landscape_patches=1, landscape_type="peripheral"),
                 f_reps(species_parameters, N0, n_aphid, n_ptoid, n_hyper, 
                        M_land, patch_state=patch_states[3,], dt, tmax, n_rep) %>% 
                   mutate(landscape_patches=4, landscape_type="central"),
                 f_reps(species_parameters, N0, n_aphid, n_ptoid, n_hyper, 
                        M_land, patch_state=patch_states[4,], dt, tmax, n_rep) %>% 
                   mutate(landscape_patches=4, landscape_type="peripheral")) %>%
    mutate(community=df_comms$community[i])
  
  write.csv(df_sim, paste0("Output/out_",df_comms$community[i],"_",n_patches,".csv"), row.names = FALSE)
}


# POSTPROCESSING - larger systems ----

# create empty dataframe
df_pop = data.frame(t=integer(),
                    patch=integer(),
                    species=character(),
                    population_size=double(),
                    dN_comm=double(),
                    dN_disp=double(),
                    replica=integer(),
                    landscape_patches=integer(),
                    landscape_type=character(),
                    community=character())
# import and combine results
for(i in 1:nrow(df_comms)){
  df_pop = rbind(df_pop,
                 fread(paste0("Output/out_",df_comms$community[i],"_",n_patches,".csv")))
}
df_pop = df_pop %>%
  mutate(
    patch=paste0("patch",patch),
    landscape_patches=as.factor(landscape_patches),
    landscape_type=as.factor(landscape_type),
    species=as.factor(species),
    community=as.factor(community),
    landscape=ifelse(landscape_type=="central", paste0(landscape_patches,"C"), paste0(landscape_patches,"P")))

df_metapop = df_pop %>%
  group_by(replica, t, landscape_patches, landscape_type, species, community, landscape) %>%
  summarise(metapopulation_size=sum(population_size)) %>%
  ungroup()

df_RC_pop = df_pop %>%
  group_by(replica, patch, landscape_patches, landscape_type, species, community) %>%
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

df_RC_metapop = df_metapop %>%
  group_by(replica, landscape_patches, landscape_type, species, community, landscape) %>%
  mutate(A_diff=(metapopulation_size+lead(metapopulation_size))/2*dt) %>%
  summarise(A_diff=sum(A_diff,na.rm=TRUE)) %>%
  left_join(., df_metapop %>% filter(t==0) %>% select(-t) %>%
              rename(metapop_t0=metapopulation_size)) %>%
  mutate(recovery=A_diff,
         A_diff=A_diff-metapop_t0*tmax) %>%
  ungroup()


# PLOTS - larger systems ----

# TIME SERIES

df_metapop = df_metapop %>%
  left_join(., df_comms) %>%
  mutate(community=as.factor(fct_reorder(community, n_species+0.01*trophic_levels)))

ggplot(data=filter(df_metapop, species=="A1"), 
       aes(x=t, y=metapopulation_size, group=interaction(replica,landscape), col=landscape)) +
  geom_line(alpha=0.1) +
  facet_wrap(~community, scales="free_y") +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  scale_color_manual(values=palette_landscape) +
  labs(x="time (day)", y="metapopulation size") +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")


df_pop = df_pop %>%
  left_join(., df_pop %>% filter(t==0) %>% 
              mutate(patch_type=ifelse(population_size==0,"empty","populated")) %>%
              select(landscape, community, species, patch, replica, patch_type)) %>%
  mutate(dN_disp_fract=dN_disp/(abs(dN_disp)+abs(dN_comm)),
         dN_disp_fract=ifelse(is.nan(dN_disp_fract),0,dN_disp_fract))

ggplot(data=filter(df_pop, species=="A1", community=="1A") %>% 
         na.omit(), 
       aes(x=t, y=dN_disp_fract, group=interaction(patch,replica), col=patch_type)) +
  geom_hline(yintercept=0) +
  geom_line(alpha=0.1) +
  facet_grid(landscape~patch_type) +
  scale_color_manual(values=c("black",col_green)) +
  labs(x="time (day)", y="relative contribution of dispersal to population change", col="patch type") +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")


# RECOVERY CREDIT

data_plot = rbind(df_RC_pop %>%
                    filter(patch_type=="empty") %>%
                    group_by(community, landscape, replica, species) %>%
                    summarise(recovery=mean(recovery)) %>%
                    ungroup() %>%
                    mutate(scale="empty patches"),
                  df_RC_pop %>%
                    filter(patch_type=="populated") %>%
                    group_by(community, landscape, replica, species) %>%
                    summarise(recovery=mean(recovery)) %>%
                    ungroup() %>%
                    mutate(scale="populated patches"),
                  df_RC_metapop %>% select(-c(landscape_patches, landscape_type, A_diff, metapop_t0)) %>%
                    mutate(scale="metapopulation")) %>%
  filter(species=="A1") %>%
  left_join(., df_comms) %>%
  mutate(community=as.factor(fct_reorder(community, n_species+0.01*trophic_levels)))
data_plot$scale =  factor(data_plot$scale, c("empty patches","populated patches","metapopulation"))


ggplot(data=filter(data_plot, scale=="empty patches"),
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(alpha=0) +
  facet_wrap(~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  lims(y=c(0,NA)) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.position="none")

ggplot(data=filter(data_plot, scale=="populated patches"),
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(alpha=0) +
  facet_wrap(~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  lims(y=c(0,NA)) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.position="none")

ggplot(data=filter(data_plot, scale=="metapopulation"),
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(alpha=0) +
  facet_wrap(~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  lims(y=c(0,NA)) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.position="none")


# STATISTICAL ANALYSIS - larger systems ----

# average across equivalent patches
# EMPTY PATCHES
data_test = df_RC_pop %>%
  filter(patch_type=="empty", species=="A1") %>%
  group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
  summarise(recovery=mean(recovery)) %>%
  ungroup() %>%
  mutate(scale="empty patches") %>%
  left_join(., df_comms)
# POPULATED PATCHES
data_test = df_RC_pop %>%
  filter(patch_type=="populated", species=="A1") %>%
  group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
  summarise(recovery=mean(recovery)) %>%
  ungroup() %>%
  mutate(scale="populated patches") %>%
  left_join(., df_comms)
# METAPOPULATION
data_test = df_RC_metapop %>%
  filter(species=="A1") %>%
  mutate(scale="metapopulation") %>%
  left_join(., df_comms)

# linear anova, recovery log(x+1) transformed (to deal with 0-values)
anova_m2 = lm(log1p(recovery) ~ landscape_patches*landscape_type*community, data=data_test)
anova(anova_m2)
performance::check_model(anova_m2) 

# average predictions
pred_landscape_patches = avg_predictions(anova_m2, variables="landscape_patches", by="landscape_patches")
pred_landscape_type = avg_predictions(anova_m2, variables="landscape_type", by="landscape_type")
pred_community = avg_predictions(anova_m2, variables="community", by="community") 

# average comparisons
comp_patches = avg_comparisons(anova_m2, variables="landscape_patches")
comp_type = avg_comparisons(anova_m2, variables="landscape_type")
comp_community_all = t(outer(pred_community$estimate, pred_community$estimate, FUN="-"))
rownames(comp_community_all) = pred_community$community
colnames(comp_community_all) = pred_community$community
comp_community_all[1,-c(2,6)] = NA
comp_community_all[2,-c(3,4,7)] = NA
comp_community_all[3,-c(5,8)] = NA
comp_community_all[4,-c(5,9)] = NA
comp_community_all[5,-c(10)] = NA
comp_community_all[6,-c(7,11)] = NA
comp_community_all[7,-c(8,9,12)] = NA
comp_community_all[8,-c(10,13)] = NA
comp_community_all[9,-c(10,14)] = NA
comp_community_all[10,-c(15)] = NA
comp_community_all[11,-c(12)] = NA
comp_community_all[12,-c(13,14)] = NA
comp_community_all[13,-c(15)] = NA
comp_community_all[14,-c(15)] = NA
comp_community_all[15,] = NA

# % change
comp_patches$estimate / pred_landscape_patches$estimate[1] * 100
comp_type$estimate / pred_landscape_type$estimate[1] * 100
com_community_change = comp_community_all / pred_community$estimate * 100
com_community_change[com_community_change>100] = 0
com_community_change[com_community_change<(-100)] = 0

# dataframe with % change for community effect
df_comm_effect = data.frame(expand.grid(community=row.names(com_community_change),
                                        community2=row.names(com_community_change)),
                            effect=as.vector(com_community_change)) %>%
  na.omit() %>%
  left_join(., df_comms) %>%
  left_join(., df_comms, by=c("community2"="community")) %>%
  mutate(community_change=ifelse((aphids.y-aphids.x)==1,"+A",
                                 ifelse((ptoids.y-ptoids.x)==1,"+P","+H")))
df_comm_effect$community_change = factor(df_comm_effect$community_change, c("+A","+P","+H"))

# plot community change effect
plot_lim = max(abs(df_comm_effect$effect))
ggplot(data=df_comm_effect, 
       aes(x=community_change, y=fct_reorder(community, n_species.x+0.01*trophic_levels.x, .desc=TRUE), 
           fill=effect)) +
  geom_point(size=6, shape=21, col="grey") +
  coord_fixed() +
  scale_fill_distiller(palette="RdBu", limits=c(-plot_lim,plot_lim), direction=1) +
  labs(x="species addition", y="community food web", fill="% change\nin recovery\ncredit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(),
        axis.text=element_text(size=8), axis.title=element_text(size=8),
        strip.text=element_text(size=10),
        legend.text=element_text(size=6), legend.title=element_text(size=8),
        legend.position="right")


# PLOTS - combined experiment & simulations ----

# import experimental results
df_RC_pop_exp = read.csv("Output/data_recovery_pop_exp.csv")
df_RC_metapop_exp = read.csv("Output/data_recovery_metapop_exp.csv")

plot_model_prediction_fullfig4(df_RC_pop, df_RC_metapop, df_RC_pop_exp, df_RC_metapop_exp,
                               df_comms, species_plot="A1")
