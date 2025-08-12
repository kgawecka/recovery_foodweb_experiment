rm(list=ls())

# load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(performance)
library(marginaleffects)
library(ggpubr)
library(lubridate)

# define palettes
palette_landscape <- c("#4ea33a80","#357029ff","#307ecd80","#1d4d7cff")
palette_community <- c("black","#4da43aff","#2f7eceff")
col_blue <- "#2f7eceff"
col_green <- "#4da43aff"
col_purple <- "#8e1558ff"

# figure widths
width_singlecol = 82
width_twothirds = 110
width_fullpage = 173


# DATA IMPORT ----

# import data
data <- read.csv("Data/experiment_data.csv")

# calculate time since start
data <- data %>%
  mutate(community=as.factor(community), 
         # insect community: BB / BB&LE / BB&LE&DR
         landscape=as.factor(landscape), 
         # inital community placement: 5 = 5-patches
         #                             4C = 4-patches-central
         #                             4P = 4-patches-peripheral
         #                             1C = 1-patch-central
         #                             1P = 1-patch-peripheral
         replica=as.factor(replica)) %>%
  group_by(community, landscape, replica) %>%
  mutate(datetime = dmy_hms(paste(date, time)),
         # combine date and time into one string
         t_day = as.numeric(difftime(datetime, dplyr::first(datetime), units="days")) 
         # time since start in days
  ) %>% 
  ungroup()

# transform into longer dataframe, remove landscape "5" and rename insects
data_pop <- data %>%
  pivot_longer(cols=c(7:16), names_to="patch_aphid", values_to="count") %>%
  mutate(patch=as.factor(sub("_.*", "", patch_aphid)),
         aphid=as.factor(sub(".*_", "", patch_aphid))) %>%
  dplyr::select(-c("patch_aphid")) %>%
  dplyr::select(community, landscape, replica, t_day, count, patch, aphid) %>%
  filter(landscape!="5") %>%
  droplevels() %>%
  mutate(aphid=as.factor(ifelse(aphid=="BRBR", "BB", "LE")), 
         community=as.factor(ifelse(community=="BRBR","1A",ifelse(community=="BRBR_LIER","2A","2A-1P"))))

# calculate metapopulation size - total count across all patches
data_metapop <- data_pop %>%
  group_by(t_day, community, landscape, replica, aphid) %>%
  summarise(metapopulation_size=sum(count)) %>%
  ungroup() %>%
  group_by(community, landscape, replica, aphid) %>%
  mutate(metapopulation_size_rel=metapopulation_size/metapopulation_size[t_day==0]) %>%
  ungroup()


# RECOVERY CREDIT ----

# crop and interpolate data to ensure equal length time series

# shortest experiment duration
tmin <- data_pop %>%
  group_by(community, landscape, replica, patch, aphid) %>%
  slice_max(t_day) %>%
  ungroup() %>%
  summarise(tmin=min(t_day)) %>% as.numeric

# population (patch) time series
data_pop_rec <- rbind(data_pop, data_pop %>%
                        dplyr::select(community, landscape, replica, patch, aphid) %>%
                        unique() %>%
                        mutate(t_day=tmin, count=NA)) %>%
  group_by(community, landscape, replica, patch, aphid) %>%
  arrange(t_day) %>%
  mutate(count_interpol=(lead(count)-lag(count))/(lead(t_day)-lag(t_day))*(t_day-lead(t_day))+lead(count),
         count=ifelse(is.na(count), count_interpol, count)) %>%
  ungroup() %>%
  dplyr::select(-count_interpol) %>%
  filter(t_day<=tmin) %>%
  na.omit()

# metapopulation time series
data_metapop_rec <- rbind(data_metapop, data_metapop %>%
                            dplyr::select(community, landscape, replica, aphid) %>%
                            unique() %>%
                            mutate(t_day=tmin, metapopulation_size=NA, metapopulation_size_rel=NA)) %>%
  group_by(community, landscape, replica, aphid) %>%
  arrange(t_day) %>%
  mutate(metapop_interpol=(lead(metapopulation_size)-lag(metapopulation_size))/(lead(t_day)-lag(t_day))*(t_day-lead(t_day))+lead(metapopulation_size),
         metapopulation_size=ifelse(is.na(metapopulation_size), metapop_interpol, metapopulation_size),
         metapoprel_interpol=(lead(metapopulation_size_rel)-lag(metapopulation_size_rel))/(lead(t_day)-lag(t_day))*(t_day-lead(t_day))+lead(metapopulation_size_rel),
         metapopulation_size_rel=ifelse(is.na(metapopulation_size_rel), metapoprel_interpol, metapopulation_size_rel)) %>%
  ungroup() %>%
  dplyr::select(-metapop_interpol, -metapoprel_interpol) %>%
  filter(t_day<=tmin) %>%
  na.omit()


# calculate recovery credit

# recovery of population
data_recovery_pop <- data_pop %>%
  group_by(community, landscape, replica, patch, aphid) %>%
  mutate(recovery_total=(count+lead(count))/2*(lead(t_day)-t_day)) %>%
  summarise(recovery_total=sum(recovery_total,na.rm=TRUE)) %>%
  left_join(., data_pop_rec %>% filter(t_day==0) %>% dplyr::select(-t_day) %>%
              rename(count_t0=count)) %>%
  left_join(., data_pop_rec %>% 
              group_by(community, landscape, replica, patch, aphid) %>%
              slice_max(t_day) %>% dplyr::select(-count) %>%
              rename(tmax=t_day)) %>%
  ungroup() %>%
  mutate(recovery=recovery_total,
         patch_type=ifelse(count_t0!=0, "populated", "empty"),
         patch_position=ifelse(patch=="patch1","centre","periphery"),
         landscape_patches=substring(landscape, 1,1),
         landscape_type=substring(landscape, 2,2),
         landscape_type=as.factor(ifelse(landscape_type=="C","central","peripheral"))) %>%
  dplyr::select(-c(recovery_total,count_t0,tmax))

# recovery of metapopulation
data_recovery_metapop <- data_metapop_rec %>%
  group_by(community, landscape, replica, aphid) %>%
  mutate(recovery_total=(metapopulation_size+lead(metapopulation_size))/2*(lead(t_day)-t_day)) %>%
  summarise(recovery_total=sum(recovery_total,na.rm=TRUE)) %>%
  left_join(., data_metapop_rec %>% filter(t_day==0) %>% dplyr::select(-t_day, -metapopulation_size_rel) %>%
              rename(metapop_t0=metapopulation_size)) %>%
  left_join(., data_metapop_rec %>% 
              group_by(community, landscape, replica, aphid) %>%
              slice_max(t_day) %>% dplyr::select(-metapopulation_size, -metapopulation_size_rel) %>%
              rename(tmax=t_day)) %>%
  ungroup() %>%
  mutate(recovery=recovery_total,
         #recovery=recovery_total-metapop_t0*tmax,
         landscape_patches=substring(landscape, 1,1),
         landscape_type=substring(landscape, 2,2),
         landscape_type=as.factor(ifelse(landscape_type=="C","central","peripheral"))) %>%
  dplyr::select(-c(recovery_total,metapop_t0,tmax))

# write out dataframes
write.csv(data_recovery_pop, "Output/data_recovery_pop_exp.csv", row.names=FALSE)
write.csv(data_recovery_metapop, "Output/data_recovery_metapop_exp.csv", row.names=FALSE)


# STATISTICAL ANALYSIS ----

# function for plotting data and anova prediction for given scale and aphid
plot_model_prediction_partfig <- function(data_test, fig_col){
  
  # linear anova, recovery log(x+1) transformed (to deal with 0-values)
  anova_mod <- lm(log1p(recovery) ~ landscape_patches*landscape_type*community, data=data_test)
  
  # effect of number of patches
  pred_landscape_patches <- avg_predictions(anova_mod, variables="landscape_patches", by="landscape_patches") %>%
    mutate(sig=ifelse(anova(anova_mod)[1,5]<0.05,"sig","NS"))
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
  pred_landscape_type <- avg_predictions(anova_mod, variables="landscape_type", by="landscape_type") %>%
    mutate(sig=ifelse(anova(anova_mod)[2,5]<0.05,"sig","NS"))
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
  pred_community <- avg_predictions(anova_mod, variables="community", by="community") %>%
    mutate(sig=ifelse(anova(anova_mod)[3,5]<0.05,"sig","NS"))
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

# function for plotting data and anova prediction for all scales and given aphid
plot_model_prediction_fullfig <- function(data_pop, data_metapop, aphid_plot){
  
  # EMPTY PATCHES
  
  # average across equivalent patches
  data_test <- data_pop %>%
    filter(patch_type=="empty", aphid==aphid_plot) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    droplevels()
  
  # plot
  plot_empty = plot_model_prediction_partfig(data_test, fig_col=col_blue)
  
  
  # POPULATED PATCHES
  
  # average across equivalent patches
  data_test <- data_pop %>%
    filter(patch_type=="populated", aphid==aphid_plot) %>%
    group_by(community, landscape_patches, landscape_type, landscape, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    droplevels()
  
  # plot
  plot_populated = plot_model_prediction_partfig(data_test, fig_col=col_green)
  
  
  # METAPOPULATION
  
  # subset data for analysis
  data_test <- data_metapop %>%
    filter(aphid==aphid_plot) %>%
    droplevels()
  
  # plot
  plot_meta = plot_model_prediction_partfig(data_test, fig_col=col_purple)
  
  # combine plots
  plot_out = ggarrange(plot_empty, plot_populated, plot_meta,
                       nrow=3, labels=c("A","B","C"))
  
  return(plot_out)
}


# --- SPECIFY SCALE & APHID ---

# specify scale for analysis:
# "population_empty" / "population_populated" / "metapopulation"
scale = "metapopulation"

# specify aphid species for analysis:
# "BB" / "LE"
aphid_sp = "LE"


# --- ANALYSIS FOR SPECIFIED SCALE & APHID ---

# subset data for analysis
if(scale=="population_empty"){
  data_test <- data_recovery_pop %>%
    filter(patch_type=="empty", aphid==aphid_sp) %>%
    group_by(community, landscape, landscape_patches, landscape_type, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    droplevels()
} else if(scale=="population_populated"){
  data_test <- data_recovery_pop %>%
    filter(patch_type=="populated", aphid==aphid_sp) %>%
    group_by(community, landscape, landscape_patches, landscape_type, replica) %>%
    summarise(recovery=mean(recovery)) %>%
    ungroup() %>%
    droplevels()
} else {
  data_test <- data_recovery_metapop %>%
    filter(aphid==aphid_sp) %>%
    droplevels()
}

# examine data
hist(data_test$recovery)
summary(data_test$recovery)

# linear anova, recovery log(x+1) transformed (to deal with 0-values)
anova_mod <- lm(log1p(recovery) ~ landscape_patches*landscape_type*community, data=data_test)
anova(anova_mod)

# check model performance
performance::check_model(anova_mod) 

# plot: all predictions
marginaleffects::plot_predictions(anova_mod, by=c("landscape_patches","landscape_type","community"), points=1) +
  geom_point(data=data_test, aes(x=landscape_patches, y=log1p(recovery), col=landscape_type), 
             position=position_jitter(width=0.1, height=0), alpha=0.2) +
  coord_cartesian(ylim=c(log1p(min(data_test$recovery)),log1p(max(data_test$recovery))))

# plot: average predictions
plot_model_prediction_partfig(data_test, fig_col="black")

# average predictions
pred_landscape_patches <- avg_predictions(anova_mod, variables="landscape_patches", by="landscape_patches")
pred_landscape_type <- avg_predictions(anova_mod, variables="landscape_type", by="landscape_type")
pred_community <- avg_predictions(anova_mod, variables="community", by="community")

# average comparisons
comp_patches = avg_comparisons(anova_mod, variables="landscape_patches")
comp_patches$estimate / pred_landscape_patches$estimate[1] * 100 # % change
comp_type = avg_comparisons(anova_mod, variables="landscape_type")
comp_type$estimate / pred_landscape_type$estimate[1] * 100 # % change
if(aphid_sp=="BB"){
  comp_community = avg_predictions(anova_mod, variables="community", hypothesis=c("b2-b1=0","b3-b2=0"))
  c(comp_community$estimate[1] / pred_community$estimate[1] * 100,
    comp_community$estimate[2] / pred_community$estimate[2] * 100) # % change
} else {
  comp_community = avg_comparisons(anova_mod, variables="community")
  comp_community$estimate / pred_community$estimate[1] * 100 # % change
}

# pariswise comparisons
pwc_patches <- marginaleffects::comparisons(anova_mod, variables="landscape_patches", by=c("landscape_type")) %>%
  mutate(sig=ifelse(p.value<=0.001,"***",ifelse(p.value<=0.01,"**",ifelse(p.value<=0.05,"*",ifelse(p.value<=0.1,".","ns")))))
pwc_type <- marginaleffects::comparisons(anova_mod, variables="landscape_type", by=c("landscape_patches")) %>%
  mutate(sig=ifelse(p.value<=0.001,"***",ifelse(p.value<=0.01,"**",ifelse(p.value<=0.05,"*",ifelse(p.value<=0.1,".","ns")))))


# PLOTS ----

# ANOVA PREDICTIONS

plot_model_prediction_fullfig(data_pop=data_recovery_pop, 
                              data_metapop=data_recovery_metapop, 
                              aphid_plot="BB")

plot_model_prediction_fullfig(data_pop=data_recovery_pop, 
                              data_metapop=data_recovery_metapop, 
                              aphid_plot="LE")


# TIME SERIES

ggplot(data=filter(data_pop, aphid=="BB") %>% na.omit(), 
       aes(x=t_day, y=count, group=interaction(replica,community), col=community)) +
  geom_point(alpha=0.3) +
  geom_line() +
  facet_grid(landscape~patch, scales="free_y") +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  scale_color_manual(values=palette_community) +
  labs(x="time (day)", y="population size") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")

ggplot(data=filter(data_pop, aphid=="LE") %>% na.omit(), 
       aes(x=t_day, y=count, group=interaction(replica,community), col=community)) +
  geom_point(alpha=0.3) +
  geom_line() +
  facet_grid(landscape~patch, scales="free_y") +
  coord_cartesian(x=c(0,NA), y=c(0,NA)) +
  scale_color_manual(values=palette_community[c(2,3)]) +
  labs(x="time (day)", y="population size") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")

df_foodweb = data_pop  %>%
  mutate(count=ifelse(is.na(count),0,count))
df_foodweb = df_foodweb %>%
  left_join(., df_foodweb %>% filter(t_day==0) %>% 
              mutate(patch_type=ifelse(count==0,"empty","populated")) %>%
              select(landscape, community, aphid, patch, replica, patch_type)) %>%
  mutate(aphids_present=ifelse(count>0 & aphid=="BB",1,
                               ifelse(count>0 & aphid=="LE",2,0))) %>%
  group_by(community, landscape, replica, patch, patch_type, t_day) %>%
  summarise(total=sum(aphids_present)) %>%
  mutate(aphids_present=ifelse(total==1,"BB only",ifelse(total==2,"LE only",ifelse(total==3,"BB&LE","none"))))
df_foodweb$aphids_present = factor(df_foodweb$aphids_present, levels=c("none","BB only","LE only","BB&LE"))

ggplot(data=filter(df_foodweb, patch_type=="empty"), 
       aes(x=t_day, y=aphids_present, group=interaction(patch,replica), col=aphids_present)) +
  geom_jitter(alpha=0.5) +
  facet_grid(community~landscape, scales="free_y") +
  scale_color_brewer(palette="Dark2") +
  labs(x="time (day)", y="aphid species present", col="aphid species present") +
  guides(colour=guide_legend(override.aes=list(alpha=1))) +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.text=element_text(size=8), legend.title=element_text(size=8),
        legend.position="bottom")


# RECOVERY CREDIT

# combine dataframes for plotting
data_plot <- rbind(data_recovery_pop %>%
                     filter(patch_type=="empty") %>%
                     group_by(community, landscape, replica, aphid) %>%
                     summarise(recovery=mean(recovery)) %>%
                     ungroup() %>%
                     mutate(scale="empty patches"),
                   data_recovery_pop %>%
                     filter(patch_type=="populated") %>%
                     group_by(community, landscape, replica, aphid) %>%
                     summarise(recovery=mean(recovery)) %>%
                     ungroup() %>%
                     mutate(scale="populated patches"),
                   data_recovery_metapop %>% select(-landscape_patches, -landscape_type) %>%
                     mutate(scale="metapopulation"))
data_plot$scale <-  factor(data_plot$scale, c("empty patches","populated patches","metapopulation"))

ggplot(data=filter(data_plot, aphid=="BB"), 
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter() +
  geom_boxplot(alpha=0) +
  facet_grid(scale~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.position="none")

ggplot(data=filter(data_plot, aphid=="LE"), 
       aes(x=landscape, y=recovery, col=landscape)) +
  geom_jitter() +
  geom_boxplot(alpha=0) +
  facet_grid(scale~community, scales="free_y") +
  scale_color_manual(values=palette_landscape) +
  labs(y="recovery credit") +
  theme(panel.background=element_rect(fill="white", colour="grey"),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        axis.title=element_text(size=8), axis.text=element_text(size=6), 
        strip.text=element_text(size=8), 
        legend.position="none")
