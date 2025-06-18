rm(list=ls())

# load packages
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(spatstat.geom)

# ONE APHID SPECIES ----
# - EXPERIMENTAL DATA ANALYSIS ----

# import experimental data

data_BRBR <- read_excel("Data/parameterization_data.xlsx", sheet = "ONE APHID - BRBR") %>% 
  dplyr::select("time", "replica", 
                "patch1_wingless", "patch1_winged", 
                "patch2_wingless", "patch2_winged") %>%
  mutate(replica=as.factor(replica))

data_LIER <- read_excel("Data/parameterization_data.xlsx", sheet = "ONE APHID - LIER") %>% 
  dplyr::select("time", "replica", 
                "patch1_wingless", "patch1_winged", 
                "patch2_wingless", "patch2_winged") %>%
  mutate(replica=as.factor(replica))


# postprocess experimental data

data_BRBR <- data_BRBR %>%
  mutate(t_day = as.numeric(difftime(time, data_BRBR$time[1], units="days")) 
         # time since start in days
  ) %>% 
  group_by(replica) %>%
  mutate(A_t = patch1_wingless,
         # number of aphids in patch 1 at time t
         A_t1 = lead(A_t),
         # number of aphids in patch 1 at time t+1
         E = lead(patch2_winged) + lead(patch2_wingless), 
         # total numer of emigrated aphids (counted at t+1)
         dt = lead(t_day) - t_day
         # delta time
  ) %>%
  ungroup()

data_LIER <- data_LIER %>%
  mutate(t_day = as.numeric(difftime(time, data_BRBR$time[1], units="days")) 
         # time since start in days
  ) %>% 
  group_by(replica) %>%
  mutate(A_t = patch1_wingless,
         # number of aphids in patch 1 at time t
         A_t1 = lead(A_t),
         # number of aphids in patch 1 at time t+1
         E = lead(patch2_winged) + lead(patch2_wingless), 
         # total numer of emigrated aphids (counted at t+1)
         dt = lead(t_day) - t_day
         # delta time
  ) %>%
  ungroup()


# plot time series

ggplot(data=data_BRBR, aes(x=t_day, y=A_t)) +
  geom_point(aes(col=replica)) +
  geom_line(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  ggtitle("BRBR")

ggplot(data=data_LIER, aes(x=t_day, y=A_t)) +
  geom_point(aes(col=replica)) +
  geom_line(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  ggtitle("LIER")


# -- growth rate (r_i) & intraspecific competition (alpha_ii) ----

# fit linear model
m_BRBR <- lm((log(A_t1+E)-log(A_t))/dt ~ A_t, data=data_BRBR)
summary(m_BRBR)
m_LIER <- lm((log(A_t1+E)-log(A_t))/dt ~ A_t, data=data_LIER)
summary(m_LIER)

# extract parameter estimates and confidence intervals
r1 <- as.numeric(m_BRBR$coefficients[1])       # BRBR growth rate
r1_CI <- as.numeric(confint(m_BRBR, level=0.95)[1,])       # growth rate confidence interval
alpha11 <- -as.numeric(m_BRBR$coefficients[2]) # BRBR intraspecific competition
alpha11_CI <- -as.numeric(confint(m_BRBR, level=0.95)[2,]) # intraspecific competition confidence interval

r2 <- as.numeric(m_LIER$coefficients[1])       # LIER growth rate
r2_CI <- as.numeric(confint(m_LIER, level=0.95)[1,])       # growth rate confidence interval
alpha22 <- -as.numeric(m_LIER$coefficients[2]) # LIER intraspecific competition
alpha22_CI <- -as.numeric(confint(m_LIER, level=0.95)[2,]) # intraspecific competition confidence interval

# plot model prediction
BRBR.predict <- cbind(na.omit(data_BRBR), predict(m_BRBR, interval='confidence'))
ggplot(data=BRBR.predict, aes(x=A_t, y=(log(A_t1+E)-log(A_t))/dt)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_line(aes(x=A_t, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  lims(x=c(0,NA)) +
  ggtitle("BRBR")

LIER.predict <- cbind(na.omit(data_LIER), predict(m_LIER, interval='confidence'))
ggplot(data=LIER.predict, aes(x=A_t, y=(log(A_t1+E)-log(A_t))/dt)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_line(aes(x=A_t, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  lims(x=c(0,NA)) +
  ggtitle("LIER")


# -- emigration (e_i) ----

# plot experimental data
ggplot(data=data_BRBR, aes(x=A_t, y=E/dt)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_smooth(col="black") +
  ggtitle("BRBR")
ggplot(data=data_LIER, aes(x=A_t, y=E/dt)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_smooth(col="black") +
  ggtitle("LIER")

# Aphids emigrate to patch 2 once they reach a threshold number in patch 1, Ae_i
# => assume E_i = 0 if A_i<Ae_i
#           E_i = e_i(A_i-Ae_i) if A_i>=Ae_i

# cumulative sums of E
data_BRBR <- data_BRBR %>%
  group_by(replica) %>%
  mutate(E_cumsum=cumsum(E/dt)) %>%
  ungroup()
data_LIER <- data_LIER %>%
  group_by(replica) %>%
  mutate(E_cumsum=cumsum(E/dt)) %>%
  ungroup()

# minimum density for emigration
Ae1 <- data_BRBR %>%
  group_by(replica) %>%
  filter(E_cumsum==0) %>%
  summarise(Ae=max(A_t)) %>%
  ungroup() %>%
  summarise(Ae_mean=mean(Ae)) %>% 
  as.numeric()  # BRBR minimum density for emigration
Ae2 <- data_LIER %>%
  group_by(replica) %>%
  filter(E_cumsum==0) %>%
  summarise(Ae=max(A_t)) %>%
  ungroup() %>%
  summarise(Ae_mean=mean(Ae)) %>% 
  as.numeric() # LIER minimum density for emigration

# create variable such that linear model is fitted through origin
data_BRBR <- data_BRBR %>% mutate(A_Ae=A_t-Ae1)
data_LIER <- data_LIER %>% mutate(A_Ae=A_t-Ae2)

# fit linear model to data where E>0, with x-axis intercept at Ae
m_BRBR <- lm(E/dt~A_Ae+0, data=filter(data_BRBR, E_cumsum>0))
m_LIER <- lm(E/dt~A_Ae+0, data=filter(data_LIER, E_cumsum>0))

# extract parameter estimates and confidence intervals
e1 <- as.numeric(m_BRBR$coefficients) # BRBR per capita emigration rate
e1_CI <- as.numeric(confint(m_BRBR, level=0.95)) # e confidence interval

e2 <- as.numeric(m_LIER$coefficients) # LIER per capita emigration rate
e2_CI <- as.numeric(confint(m_LIER, level=0.95)) # e confidence interval

# plot model prediction
BRBR.predict <- cbind(filter(data_BRBR, E_cumsum>0), predict(m_BRBR, interval='confidence'))
ggplot(data=BRBR.predict, aes(x=A_Ae, y=E/dt)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_line(aes(x=A_Ae, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  lims(x=c(0,NA)) +
  ggtitle("BRBR")

LIER.predict <- cbind(filter(data_LIER, E_cumsum>0), predict(m_LIER, interval='confidence'))
ggplot(data=LIER.predict, aes(x=A_Ae, y=E/dt)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_line(aes(x=A_Ae, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  ggtitle("LIER")


# - MODEL PREDICTION ----
# model prediction using confidence intervals of estimated parameters

# define model function

# inputs:
# r_CI - growth rate confidence interval
# alpha_CI = intraspecific competition confidence interval
# Ae = minimum number of aphids for emigration
# e_CI = emigration rate confidence interval
# A0 = initial number of aphids
# dt = time step size (days)
# tmax = maximum time for simulation (days)
# n_it = number of model iterations
dt_one_aphid_CI = function(r_CI, alpha_CI, Ae, e_CI, A0, dt, tmax, n_it){
  
  time <- seq(0,tmax,dt) # timestep vector
  n_dt <- length(time)   # number of timesteps
  
  # initialise dataframe for storing results
  df_out <- data.frame(expand.grid(iteration=1:n_it, t_day=time),
                       A_t=NA,
                       r=NA, alpha=NA, e=NA)
  
  # assing 0 to negative parameter confidence intervals
  r_CI <- ifelse(r_CI<0,0,r_CI)
  alpha_CI <- ifelse(alpha_CI<0,0,alpha_CI)
  e_CI <- ifelse(e_CI<0,0,e_CI)
  
  # sample parameter values from uniform distributions using confidence intervals
  set.seed(1)
  r_vals <- runif(n_it, min=min(r_CI), max=max(r_CI))
  alpha_vals <- runif(n_it, min=min(alpha_CI), max=max(alpha_CI))
  e_vals <- runif(n_it, min=min(e_CI), max=max(e_CI))
  
  # model iterations loop
  for(n in 1:n_it){
    
    # sample parameter values from uniform distributions using confidence intervals
    r <- r_vals[n]
    alpha <- alpha_vals[n]
    e <- e_vals[n]
    
    # initialise vector for storing number of aphids
    A = rep(NA,n_dt)
    A[1] = A0 # initial number of aphids
    
    # iterate through timesteps
    for(t in 1:(n_dt-1)){
      
      if(is.na(A[t])){next}
      
      # delta for emigration
      delta = ifelse(A[t]<Ae, 0, 1)
      
      # change in A due to emigration
      dA_emig = delta*e*(A[t]-Ae)
      
      # change in A
      dA = A[t]*exp((r-alpha*A[t]) * dt) - A[t] - dA_emig
      
      # new A
      A[t+1] = A[t] + dA
      
      if(is.infinite(A[t+1]) || is.nan(A[t+1])) {A[t+1]=NA}
      if(!is.na(A[t+1]) && A[t+1]<0) {A[t+1]=0}
      
    }
    
    # store result in dataframe
    df_out[df_out$iteration==n,"A_t"] = A
    df_out[df_out$iteration==n,"r"] = r
    df_out[df_out$iteration==n,"alpha"] = alpha
    df_out[df_out$iteration==n,"e"] = e
    
  }
  
  return(df_out)
}


# run model
dt_CI_BRBR = dt_one_aphid_CI(r_CI=r1_CI, alpha_CI=alpha11_CI, Ae=Ae1, e_CI=e1_CI, 
                             A0=5, dt=1, tmax=30, n_it=1000)
dt_CI_LIER = dt_one_aphid_CI(r_CI=r2_CI, alpha_CI=alpha22_CI, Ae=Ae2, e_CI=e2_CI, 
                             A0=5, dt=1, tmax=30, n_it=1000)

# round up time in experimental data
data_BRBR_average <- data_BRBR %>%
  mutate(t_day=round(t_day,0))
data_LIER_average <- data_LIER %>%
  mutate(t_day=round(t_day,0))

# calculate difference between experiment and predictions
dt_CI_BRBR_diff <- dt_CI_BRBR %>%
  left_join(., data_BRBR_average %>% dplyr::select(t_day,replica,A_t) %>% rename(A_data="A_t")) %>%
  mutate(A_diff=abs(A_t-A_data)) %>%
  group_by(iteration) %>%
  summarise(error=sum(A_diff, na.rm=TRUE)) %>%
  ungroup() %>%
  left_join(., unique(dplyr::select(dt_CI_BRBR, c(iteration, r, alpha, e)))) %>%
  mutate(weight=1/(min(error)-max(error))*(error-max(error)),
         top_50=ifelse(error<median(error), "Y", "N"))
dt_CI_LIER_diff <- dt_CI_LIER %>%
  left_join(., data_LIER_average %>% dplyr::select(t_day,replica,A_t) %>% rename(A_data="A_t")) %>%
  mutate(A_diff=abs(A_t-A_data)) %>%
  group_by(iteration) %>%
  summarise(error=sum(A_diff, na.rm=TRUE)) %>%
  ungroup() %>%
  left_join(., unique(dplyr::select(dt_CI_LIER, c(iteration, r, alpha, e)))) %>%
  mutate(weight=1/(min(error)-max(error))*(error-max(error)),
         top_50=ifelse(error<median(error), "Y", "N"))


# TWO APHID SPECIES ----
# - EXPERIMENTAL DATA ANALYSIS ----

# import experimental data
data_2aphids <- read_excel("Data/parameterization_data.xlsx", sheet = "TWO APHIDS")  %>% 
  dplyr::select("time", "replica", "BRBR_wingless", "LIER_wingless") %>%
  mutate(replica=as.factor(replica))

# postprocess experimental data
data_2aphids <- data_2aphids %>%
  mutate(t_day = as.numeric(difftime(time, data_2aphids$time[1], units="days")), 
         # time since start in days
  ) %>% 
  group_by(replica) %>%
  mutate(A1_t = BRBR_wingless,
         # total number of BRBR aphids at time t
         A1_t1 = lead(A1_t),
         # total number of BRBR aphids at time t+1
         A2_t = LIER_wingless,
         # total number of LIER aphids at time t
         A2_t1 = lead(A2_t),
         # total number of LIER aphids at time t+1
         dt = lead(t_day) - t_day  
         # delta time
  ) %>%
  ungroup()


# plot time series

ggplot(data=data_2aphids, aes(x=t_day, y=A1_t)) +
  geom_point(aes(col=replica)) +
  geom_line(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  ggtitle("BRBR")

ggplot(data=data_2aphids, aes(x=t_day, y=A2_t)) +
  geom_point(aes(col=replica)) +
  geom_line(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  ggtitle("LIER")


# -- interspecific competition (alpha_ij) ----

# calculate response variable using extimated r_i and alpha_ii
data_2aphids <- data_2aphids %>%
  mutate(Y1 = r1 - alpha11*A1_t - (log(A1_t1)-log(A1_t))/dt,
         Y2 = r2 - alpha22*A2_t - (log(A2_t1)-log(A2_t))/dt)

# fit linear model
m_BRBR <- lm(Y1 ~ A2_t+0, data=data_2aphids)
summary(m_BRBR)
m_LIER <- lm(Y2 ~ A1_t+0, data=data_2aphids)
summary(m_LIER)

# extract parameter estimates and confidence intervals
alpha12 <- as.numeric(m_BRBR$coefficients) # BRBR interspecific competition
alpha12_CI <- as.numeric(confint(m_BRBR, level=0.95)) # interspecific competition confidence interval
alpha21 <- as.numeric(m_LIER$coefficients) # LIER interspecific competition
alpha21_CI <- as.numeric(confint(m_LIER, level=0.95)) # interspecific competition confidence interval

# plot model prediction
BRBR.predict <- cbind(na.omit(data_2aphids %>% select(replica, A2_t, Y1)), predict(m_BRBR, interval='confidence'))
ggplot(data=BRBR.predict, aes(x=A2_t, y=Y1)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_line(aes(x=A2_t, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  ggtitle("BRBR")

LIER.predict <- cbind(na.omit(data_2aphids %>% select(replica, A1_t, Y2)), predict(m_LIER, interval='confidence'))
ggplot(data=LIER.predict, aes(x=A1_t, y=Y2)) +
  geom_hline(yintercept=0) +
  geom_point(aes(col=replica)) +
  geom_line(aes(x=A1_t, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  ggtitle("LIER")


# - MODEL PREDICTION ----
# model prediction using confidence intervals of estimated parameters

# define model function

# inputs:
# r1_CI - BRBR growth rate confidence interval
# r2_CI - LIER growth rate confidence interval
# alpha11_CI = BRBR intraspecific competition confidence interval
# alpha12_CI = BRBR interspecific competition confidence interval
# alpha22_CI = LIER intraspecific competition confidence interval
# alpha21_CI = LIER interspecific competition confidence interval
# Ae1 = BRBR minimum number of aphids for emigration
# Ae2 = LIER minimum number of aphids for emigration
# e1_CI = BRBR emigration rate confidence interval
# e2_CI = LIER emigration rate confidence interval
# A10 = BRBR initial number of aphids
# A20 = LIER initial number of aphids
# dt = time step size (days)
# tmax = maximum time for simulation (days)
# n_it = number of model iterations
dt_two_aphid_CI = function(r1_CI, r2_CI, 
                           alpha11_CI, alpha22_CI, alpha12_CI, alpha21_CI,
                           Ae1, Ae2, e1_CI, e2_CI,
                           A10, A20, dt, tmax, n_it){
  
  time = seq(0,tmax,dt) # timestep vector
  n_dt = length(time)   # number of timesteps
  
  # initialise dataframe for storing results
  df_out <- data.frame(expand.grid(iteration=1:n_it, t_day=time),
                       A1_t=NA, A2_t=NA,
                       r1=NA, r2=NA, 
                       alpha11=NA, alpha22=NA, alpha12=NA, alpha21=NA,
                       e1=NA, e2=NA)
  
  # assign 0 to negative parameter confidence intervals
  r1_CI <- ifelse(r1_CI<0,0,r1_CI)
  r2_CI <- ifelse(r2_CI<0,0,r2_CI)
  alpha11_CI <- ifelse(alpha11_CI<0,0,alpha11_CI)
  alpha22_CI <- ifelse(alpha22_CI<0,0,alpha22_CI)
  alpha12_CI <- ifelse(alpha12_CI<0,0,alpha12_CI)
  alpha21_CI <- ifelse(alpha21_CI<0,0,alpha21_CI)
  e1_CI <- ifelse(e1_CI<0,0,e1_CI)
  e2_CI <- ifelse(e2_CI<0,0,e2_CI)
  
  # sample parameter values from uniform distributions using confidence intervals
  set.seed(1)
  r1_vals <- runif(n_it, min=min(r1_CI), max=max(r1_CI))
  alpha11_vals <- runif(n_it, min=min(alpha11_CI), max=max(alpha11_CI))
  e1_vals <- runif(n_it, min=min(e1_CI), max=max(e1_CI))
  
  set.seed(1)
  r2_vals <- runif(n_it, min=min(r2_CI), max=max(r2_CI))
  alpha22_vals <- runif(n_it, min=min(alpha22_CI), max=max(alpha22_CI))
  e2_vals <- runif(n_it, min=min(e2_CI), max=max(e2_CI))
  
  set.seed(1)
  alpha12_vals <- runif(n_it, min=min(alpha12_CI), max=max(alpha12_CI))
  alpha21_vals <- runif(n_it, min=min(alpha21_CI), max=max(alpha21_CI))
  
  # model iterations loop
  for(n in 1:n_it){
    
    # sample parameter values from uniform distributions using confidence intervals
    r1 <- r1_vals[n]
    r2 <- r2_vals[n]
    alpha11 <- alpha11_vals[n]
    alpha22 <- alpha22_vals[n]
    alpha12 <- alpha12_vals[n]
    alpha21 <- alpha21_vals[n]
    e1 <- e1_vals[n]
    e2 <- e2_vals[n]
    
    # initialise vector for storing number of aphids
    A1 = rep(NA,n_dt)
    A1[1] = A10 # initial number of aphids
    A2 = rep(NA,n_dt)
    A2[1] = A20 # initial number of aphids
    
    # iterate through timesteps
    for(t in 1:(n_dt-1)){
      
      if(is.na(A1[t]) || is.na(A2[t])){next}
      
      # delta for emigration
      delta1 = ifelse(A1[t]<Ae1, 0, 1)
      delta2 = ifelse(A2[t]<Ae2, 0, 1)
      
      # change in A due to emigration
      dA_emig1 = delta1*e1*(A1[t]-Ae1)*dt
      dA_emig2 = delta2*e2*(A2[t]-Ae2)*dt
      
      # change in A
      dA1 = A1[t]*exp((r1-alpha11*A1[t]-alpha12*A2[t])*dt) - A1[t] - dA_emig1
      dA2 = A2[t]*exp((r2-alpha22*A2[t]-alpha21*A1[t])*dt) - A2[t] - dA_emig2
      
      # new A
      A1[t+1] = A1[t] + dA1
      A2[t+1] = A2[t] + dA2
      
      if(is.infinite(A1[t+1]) || is.nan(A1[t+1])) {A1[t+1]=NA}
      if(is.infinite(A2[t+1]) || is.nan(A2[t+1])) {A2[t+1]=NA}
      if(!is.na(A1[t+1]) && A1[t+1]<0) {A1[t+1]=0}
      if(!is.na(A2[t+1]) && A2[t+1]<0) {A2[t+1]=0}
      
    }
    
    # store result in dataframe
    df_out[df_out$iteration==n,"A1_t"] = A1
    df_out[df_out$iteration==n,"A2_t"] = A2
    df_out[df_out$iteration==n,"r1"] = r1
    df_out[df_out$iteration==n,"r2"] = r2
    df_out[df_out$iteration==n,"alpha11"] = alpha11
    df_out[df_out$iteration==n,"alpha22"] = alpha22
    df_out[df_out$iteration==n,"alpha12"] = alpha12
    df_out[df_out$iteration==n,"alpha21"] = alpha21
    df_out[df_out$iteration==n,"e1"] = e1
    df_out[df_out$iteration==n,"e2"] = e2
    
  }
  
  return(df_out)
}

# run model
# e = 0 - no emigration in experiment
dt_CI_2aphids = dt_two_aphid_CI(r1=r1_CI, r2=r2_CI, 
                                alpha11_CI=alpha11_CI, alpha22_CI=alpha22_CI,
                                alpha12_CI=alpha12_CI, alpha21_CI=alpha21_CI,
                                Ae1=0, Ae2=0, e1=c(0,0), e2=c(0,0),
                                A10=5, A20=5, dt=1, tmax=30, n_it=1000)

# round up time in experimental data
data_2aphids_average <- data_2aphids %>%
  mutate(t_day=round(t_day,0))

# calculate difference between experiment and predictions
dt_CI_2aphids_diff <- dt_CI_2aphids %>%
  left_join(., data_2aphids_average %>% dplyr::select(t_day,replica,A1_t,A2_t) %>% rename(A1_data="A1_t", A2_data="A2_t")) %>%
  mutate(A1_diff=abs(A1_t-A1_data),
         A2_diff=abs(A2_t-A2_data)) %>%
  group_by(iteration) %>%
  summarise(error1=sum(A1_diff, na.rm=TRUE),
            error2=sum(A2_diff, na.rm=TRUE)) %>%
  mutate(error=error1+error2) %>%
  ungroup() %>%
  left_join(unique(dplyr::select(dt_CI_2aphids, c(iteration, r1, r2, alpha11, alpha22, alpha12, alpha21)))) %>%
  mutate(weight=1/(min(error)-max(error))*(error-max(error)),
         top_50=ifelse(error<median(error), "Y", "N"))


# ONE & TWO APHID SPECIES ----
# - ITERATED PARAMETERS ----

# combine errors from one aphid and two aphid simulations
df_error <- rbind(dt_CI_2aphids_diff %>% dplyr::select(-error1,-error2) %>%
                    mutate(e1=NA,e2=NA, community="BRBR-LIER"),
                  dt_CI_BRBR_diff %>% 
                    rename(r1="r",alpha11="alpha",e1="e") %>%
                    mutate(r2=NA,alpha22=NA,alpha12=NA,alpha21=NA,e2=NA, community="BRBR"),
                  dt_CI_LIER_diff %>% 
                    rename(r2="r",alpha22="alpha",e2="e") %>%
                    mutate(r1=NA,alpha11=NA,alpha12=NA,alpha21=NA,e1=NA, community="LIER"))

# weighted distribution of parameters
df_error_plot <- df_error %>%
  pivot_longer(cols=c("r1","r2","alpha11","alpha22","alpha12","alpha21","e1","e2"), 
               names_to="parameter", values_to="value")

ggplot(data=filter(df_error_plot, top_50=="Y"),
       aes(x=parameter, y=value)) +
  geom_jitter(aes(size=weight), alpha=0.1) +
  geom_boxplot(aes(weight=weight), alpha=0) +
  facet_wrap(~parameter, scales="free")

# summarise parameters in top 50% simulations
df_param_sims <- df_error_plot %>%
  filter(top_50=="Y") %>%
  group_by(parameter) %>%
  summarise(min=min(value, na.rm=TRUE),
            Q1=weighted.quantile(value, w=weight, probs=0.25, na.rm=TRUE),
            median=weighted.median(value, w=weight, na.rm=TRUE),
            Q3=weighted.quantile(value, w=weight, probs=0.75, na.rm=TRUE),
            max=max(value, na.rm=TRUE)) %>%
  ungroup()

r1_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="r1")$Q1,
                          filter(df_param_sims, parameter=="r1")$Q3))
r2_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="r2")$Q1,
                          filter(df_param_sims, parameter=="r2")$Q3))
alpha11_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="alpha11")$Q1,
                               filter(df_param_sims, parameter=="alpha11")$Q3))
alpha22_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="alpha22")$Q1,
                               filter(df_param_sims, parameter=="alpha22")$Q3))
alpha12_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="alpha12")$Q1,
                               filter(df_param_sims, parameter=="alpha12")$Q3))
alpha21_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="alpha21")$Q1,
                               filter(df_param_sims, parameter=="alpha21")$Q3))
e1_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="e1")$Q1,
                          filter(df_param_sims, parameter=="e1")$Q3))
e2_CI_itr <- as.numeric(c(filter(df_param_sims, parameter=="e2")$Q1,
                          filter(df_param_sims, parameter=="e2")$Q3))


# - MODEL PREDICTION ----
# model prediction using confidence intervals of iterated parameters

# run model
dt_CI_BRBR = dt_one_aphid_CI(r_CI=r1_CI_itr, alpha_CI=alpha11_CI_itr,
                             Ae=Ae1, e_CI=e1_CI_itr, 
                             A0=5, dt=1, tmax=30, n_it=1000)

dt_CI_LIER = dt_one_aphid_CI(r_CI=r2_CI_itr, alpha_CI=alpha22_CI_itr, 
                             Ae=Ae2, e_CI=e2_CI_itr, 
                             A0=5, dt=1, tmax=30, n_it=1000)

dt_CI_2aphids = dt_two_aphid_CI(r1=r1_CI_itr, r2=r2_CI_itr, 
                                alpha11_CI=alpha11_CI_itr, 
                                alpha22_CI=alpha22_CI_itr,
                                alpha12_CI=alpha12_CI_itr,
                                alpha21_CI=alpha21_CI_itr,
                                Ae1=0, Ae2=0, e1=c(0,0), e2=c(0,0),
                                A10=5, A20=5, dt=1, tmax=30, n_it=1000)


# plot experimental data and model predictions

p1 <- ggplot(data=data_BRBR, aes(x=t_day, y=A_t)) +
  geom_line(data=dt_CI_BRBR, aes(group=iteration), alpha=0.1) +
  geom_point(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,max(data_BRBR$A_t, na.rm=TRUE)+100)) +
  scale_color_brewer(palette="Dark2") +
  labs(x="time (day)", y="population size", subtitle="BB") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), 
        axis.text=element_text(size=8), axis.title=element_text(size=8),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        plot.subtitle=element_text(size=10), legend.position="bottom")

p2 <- ggplot(data=data_LIER, aes(x=t_day, y=A_t)) +
  geom_line(data=dt_CI_LIER, aes(group=iteration), alpha=0.1) +
  geom_point(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,max(data_LIER$A_t, na.rm=TRUE)+100)) +
  scale_color_brewer(palette="Dark2") +
  labs(x="time (day)", y="population size", subtitle="LE") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), 
        axis.text=element_text(size=8), axis.title=element_text(size=8),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        plot.subtitle=element_text(size=10), legend.position="bottom")

p3 <- ggplot(data=data_2aphids, aes(x=t_day, y=A1_t)) +
  geom_line(data=dt_CI_2aphids, aes(group=iteration), alpha=0.1) +
  geom_point(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,max(data_2aphids$A1_t, na.rm=TRUE)+100)) +
  scale_color_brewer(palette="Dark2") +
  labs(x="time (day)", y="population size", subtitle="BB-LE: BB") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), 
        axis.text=element_text(size=8), axis.title=element_text(size=8),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        plot.subtitle=element_text(size=10), legend.position="bottom")

p4 <- ggplot(data=data_2aphids, aes(x=t_day, y=A2_t)) +
  geom_line(data=dt_CI_2aphids, aes(group=iteration), alpha=0.1) +
  geom_point(aes(col=replica)) +
  coord_cartesian(x=c(0,NA), y=c(0,max(data_2aphids$A2_t, na.rm=TRUE)+100)) +
  scale_color_brewer(palette="Dark2") +
  labs(x="time (day)", y="population size", subtitle="BB-LE: LE") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background=element_blank(), 
        axis.text=element_text(size=8), axis.title=element_text(size=8),
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        plot.subtitle=element_text(size=10), legend.position="bottom")

ggarrange(p1,p2,p3,p4, nrow=2, ncol=2, common.legend=TRUE, legend="bottom")


# PARASITOID ----
# - EXPERIMENTAL DATA ANALYSIS ----
# -- maximum parasitization ----

# import experimental data
data_ptoid_max_parasit <- read_excel("Data/parameterization_data.xlsx", sheet = "PTOID - MAX PARASIT") %>%
  dplyr::select("aphid_species", "replica", "no_mummies") %>%
  mutate(replica=as.factor(replica))

# plot experimental data
ggplot(data=data_ptoid_max_parasit,
       aes(x=aphid_species, y=no_mummies, col=replica)) +
  geom_jitter()

# maximum parasitization
pmax <- mean(data_ptoid_max_parasit$no_mummies)


# -- functional response ----

# import experimental data
data_ptoid_function <- read_excel("Data/parameterization_data.xlsx", sheet = "PTOID - FUNCTION") %>%   
  dplyr::select("aphid_species", "initial_no_aphids", "replica", "no_mummies") %>%   
  mutate(replica=as.factor(replica))

# plot experimental data
ggplot(data=data_ptoid_function, aes(x=initial_no_aphids, y=no_mummies)) +
  geom_jitter(aes(col=replica)) +   
  geom_smooth(col="black") +   
  facet_wrap(~aphid_species, scales="free_y") +   
  coord_cartesian(y=c(0,NA))

# fit linear model, with x-axis intercept at 0
m_ptoid_BRBR <- lm(formula = no_mummies ~ initial_no_aphids + 0, 
                   data = filter(data_ptoid_function, aphid_species=="BRBR"))
summary(m_ptoid_BRBR)
m_ptoid_LIER <- lm(formula = no_mummies ~ initial_no_aphids + 0, 
                   data = filter(data_ptoid_function, aphid_species=="LIER"))
summary(m_ptoid_LIER)

# extract fitted parameters
beta1 <- as.numeric(m_ptoid_BRBR$coefficients) # BRBR parasitization rate
beta1_CI <- as.numeric(confint(m_ptoid_BRBR, level=0.95)) # parasitization rate confidence interval
beta2 <- as.numeric(m_ptoid_LIER$coefficients) # LIER parasitization rate
beta2_CI <- as.numeric(confint(m_ptoid_LIER, level=0.95)) # parasitization rate confidence interval

# plot model prediction
ptoid.predict <- cbind(filter(data_ptoid_function, aphid_species=="BRBR"), predict(m_ptoid_BRBR, interval='confidence'))
ggplot(data=ptoid.predict, aes(x=initial_no_aphids, y=no_mummies)) +
  geom_hline(yintercept=0) +
  geom_jitter(aes(col=as.factor(replica)), width=0.1, height=0.1) + 
  geom_line(aes(x=initial_no_aphids, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  ggtitle("BRBR")

ptoid.predict <- cbind(filter(data_ptoid_function, aphid_species=="LIER"), predict(m_ptoid_LIER, interval='confidence'))
ggplot(data=ptoid.predict, aes(x=initial_no_aphids, y=no_mummies)) +
  geom_hline(yintercept=0) +
  geom_jitter(aes(col=as.factor(replica)), width=0.1, height=0.1) + 
  geom_line(aes(x=initial_no_aphids, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  ggtitle("LIER")


# -- time to emergence ----

# import experimental data
data_ptoid_emergence <- read_excel("Data/parameterization_data.xlsx", sheet = "PTOID - EMERGENCE") %>%
  dplyr::select("date_setup", "date_emerged_count", "aphid_species", "no_emerged") %>%
  mutate(date_setup=as.Date(date_setup), date_emerged_count=as.Date(date_emerged_count)) %>%
  mutate(tau = as.numeric(date_emerged_count-date_setup),
         # temporal delay between parasitoid attack and new parasitoid emergence (days)
  ) %>% 
  filter(no_emerged>0)

# plot experimental data
ggplot(data=data_ptoid_emergence, aes(x=tau, y=no_emerged)) +
  geom_bar(stat="identity") +
  facet_wrap(~aphid_species)

# time to emergence
tau <- weighted.mean(data_ptoid_emergence$tau, data_ptoid_emergence$no_emerged)


# -- fraction of females ----

# import experimental data
data_ptoid_females <- read_excel("Data/parameterization_data.xlsx", sheet = "PTOID - FEMALES") %>%
  dplyr::select("date", "batch", "no_ptoids", "no_females") %>%
  mutate(date=as.Date(date), batch=as.factor(batch))

# fraction of females
f <- sum(data_ptoid_females$no_females) / sum(data_ptoid_females$no_ptoids)


# -- mortality ----

# import experimental data
data_ptoid_mortality <- read_excel("Data/parameterization_data.xlsx", sheet = "PTOID - MORTALITY")  %>% 
  dplyr::select("date", "replica", "no_ptoids") %>%
  mutate(date=as.Date(date), replica=as.factor(replica))
data_ptoid_mortality <- data_ptoid_mortality %>%
  mutate(age = as.numeric(date-data_ptoid_mortality$date[1]+1)
         # ptoid age (days); note: 1 day old at setup
  ) %>% 
  group_by(replica) %>%
  mutate(no_dead = first(no_ptoids) - no_ptoids
         # no of dead parasitoids
  ) %>%
  ungroup()

# plot experimental data
ggplot(data=data_ptoid_mortality, aes(x=age, y=no_dead, fill=replica)) +
  geom_bar(stat="identity")

# lifespan
lambda <- weighted.mean(data_ptoid_mortality$age, data_ptoid_mortality$no_dead)


# -- dispersal ----

# import experimental data
data_ptoid_dispersal <- read_excel("Data/parameterization_data.xlsx", sheet = "PTOID - DISPERSAL")  %>% 
  dplyr::select("date", "initial_no_ptoids", "replica", "patch1_ptoids", "patch2_ptoids") %>%
  mutate(date=as.Date(date), replica=as.factor(replica))
data_ptoid_dispersal <- data_ptoid_dispersal %>%
  mutate(t_day = as.numeric(date-data_ptoid_dispersal$date[1]),
         # time since start in days
         P = patch1_ptoids
         # number of parasitoids in patch 1
  ) %>% 
  group_by(initial_no_ptoids, replica) %>%
  mutate(dt = lead(t_day) - t_day,  # delta time
         E = (lead(patch2_ptoids) - patch2_ptoids)/dt
         # emigration rate to patch 2
  ) %>%
  ungroup()

# plot experimental data
ggplot(data=na.omit(data_ptoid_dispersal), aes(x=P, y=E)) +   
  geom_hline(yintercept=0) +   
  geom_jitter(aes(col=as.factor(t_day)), width=0.1, height=0.1) +   
  geom_smooth(col="black") +
  labs(col="t_day")

# Parasitoids emigrate to patch 2 once they reach a threshold number in patch 1, Pe
# => assume E = 0 if P<Pe
#           E = e_P(P-Pe) if P>=Pe

# fit linear model to determine minimum density for emigration
m_ptoid <- lm(E~P, data=data_ptoid_dispersal)
eP <- as.numeric(m_ptoid$coefficients[2]) # slope
e_I <- as.numeric(m_ptoid$coefficients[1]) # intercept with E-axis
Pe <- -e_I/eP # intercept with P-axis - minimum density for emigration

# create variable such that linear model is fitted through origin
data_ptoid_dispersal <- mutate(data_ptoid_dispersal, 
                               P_Pe=P-Pe) 

# fit linear model with x-axis intercept at Pe
m_ptoid <- lm(E~P_Pe+0, data=data_ptoid_dispersal)
summary(m_ptoid)

# extract parameter estimates
eP <- as.numeric(m_ptoid$coefficients) # emigration rate
eP_CI <- as.numeric(confint(m_ptoid)) # emigration rate confidence interval

# plot model prediction
ptoid.predict <- cbind(na.omit(data_ptoid_dispersal), predict(m_ptoid, interval='confidence'))
ggplot(data=ptoid.predict, aes(x=P_Pe, y=E)) +
  geom_hline(yintercept=0) +
  geom_jitter(aes(col=as.factor(t_day)), width=0.1, height=0.1) + 
  geom_line(aes(x=P_Pe, y=fit), col="black", size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr), alpha=0.3) +
  labs(col="t_day")


# SUMMARY OF PARAMETERS ----

df_parameters <- data.frame(parameter=c("r1","r2","alpha11","alpha22","alpha12","alpha21",
                                        "e1","e2","Ae1","Ae2","beta1","beta2",
                                        "pmax","tau","f","lambda","eP","Pe"),
                            value=c(NA,NA,NA,NA,NA,NA,
                                    NA,NA,Ae1,Ae2,NA,NA,
                                    pmax,tau,f,lambda,eP,Pe),
                            CI_low=c(r1_CI_itr[1],r2_CI_itr[1],
                                     alpha11_CI_itr[1],alpha22_CI_itr[1],
                                     alpha12_CI_itr[1],alpha21_CI_itr[1],
                                     e1_CI_itr[1],e2_CI_itr[1],NA,NA,
                                     beta1_CI[1],beta2_CI[1],
                                     NA,NA,NA,NA,NA,NA),
                            CI_high=c(r1_CI_itr[2],r2_CI_itr[2],
                                      alpha11_CI_itr[2],alpha22_CI_itr[2],
                                      alpha12_CI_itr[2],alpha21_CI_itr[2],
                                      e1_CI_itr[2],e2_CI_itr[2],NA,NA,
                                      beta1_CI[2],beta2_CI[2],
                                      NA,NA,NA,NA,NA,NA))

write.csv(df_parameters, paste0("Data/model_parameters.csv", sep=""), row.names=FALSE)
