rm(list=ls())
source("0_Functions.R")
# --------------------------- Design final figures ----------------------------------------


# Figure 2

library(ggraph)
d=readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model.rds")
p1=Plot_variance_partitioning(d$SPEI06$Coeff_variance,alpha1 = 1,alpha2 = 1)
p2=Plot_variance_partitioning(d$SPEI06$Coeff_variance,alpha1 = 1,alpha2 = 1,
                              stability_var = "Resilience_Isbell_abs")

d=readRDS("./Results/Moving_window_gradient_all_results_SPEI.rds")
p3=Plot_variance_partitioning(d$SPEI06$Coeff_variance,alpha1 = 1,alpha2 = 1,negative_x = T)
p4=Plot_variance_partitioning(d$SPEI06$Coeff_variance,alpha1 = 1,alpha2 = 1,negative_x = T,
                              stability_var = "Resilience_Isbell_abs")

figure2=(ggarrange(p1$p2+labs(y="Explained variance (%)",x="Land-use intensity")+
                     geom_node_label(
                       data = NULL,
                       aes(x = 1.7, y = 1.05, label = "Drought resistance"),
                       color = "white",
                       fill="black",
                       label.size = 1,   
                       family = "NewCenturySchoolbook",
                       label.padding = unit(.5, "lines"),
                       label.r = unit(0.3, "lines"),  
                       size = 4.5             
                     )+ guides(fill = guide_legend(nrow = 1)),
                   p2$p2+labs(y="Explained variance (%)",x="Land-use intensity")+
                     geom_node_label(
                       data = NULL,
                       aes(x = 1.7, y = 1.05, label = "Drought resilience"),
                       color = "white",
                       fill="black",
                       label.size = 1,   
                       family = "NewCenturySchoolbook",
                       label.padding = unit(.5, "lines"),
                       label.r = unit(0.3, "lines"),  
                       size = 4.5             
                     ),
                   p3$p2+labs(y="Explained variance (%)",x="Drought intensity \n <- Moderate drought      Extreme drought ->"),
                   p4$p2+labs(y="Explained variance (%)",x="Drought intensity \n <- Moderate drought      Extreme drought ->"),nrow=2,ncol=2,common.legend = T,legend = "bottom"))

ggsave("./Figures/Figure_2.pdf",figure2,width = 8,height = 8)



# Figure 3

d=readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model.rds")

p_CWM = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI06$Partial_residuals,
  d_Rl = d$SPEI06$Partial_residuals,
  grepl_character = "CWM"
)

p_FD = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI06$Partial_residuals,
  d_Rl = d$SPEI06$Partial_residuals,
  grepl_character = "FD"
)
ggsave("./Figures/Figure_3.pdf",ggarrange(p_CWM+ylab("Effect size of mean community trait"),
                                          p_FD+ylab("Effect size of trait diversity"),
                                          nrow = 2,labels = letters[1:2],common.legend = T,legend = "bottom"),
       width = 8,height = 8)

# Figure 4

d=readRDS("./Results/Moving_window_gradient_all_results_SPEI.rds")

p_FD = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI06$Partial_residuals,
  d_Rl = d$SPEI06$Partial_residuals,negative_x = T,
  grepl_character = "FD"
)
ggsave("./Figures/Figure_4.pdf",p_FD+ylab("Effect size of trait diversity")+
         xlab("Drought intensity \n <- Moderate drought      Extreme drought ->"),
       width = 8,height = 4)


# Figure 5

d = readRDS("./Results/Paired_results.rds")

join_responses_boot = function(res_obj, rel_obj) {
  
  res = res_obj$Partial_residuals %>%
    dplyr::filter(Stability_var == "Resistance_Isbell") %>%
    dplyr::select(-Stability_var) %>%
    tidyr::pivot_longer(-Variable_gradient,
                        names_to = "Predictor",
                        values_to = "Effect_size_Rs") %>%
    dplyr::filter(grepl("CWM|FD", Predictor), !is.na(Effect_size_Rs)) %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::mutate(N_boot = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  rel = rel_obj$Partial_residuals %>%
    dplyr::filter(Stability_var == "Resilience_Isbell_abs") %>%
    dplyr::select(-Stability_var) %>%
    tidyr::pivot_longer(-Variable_gradient,
                        names_to = "Predictor",
                        values_to = "Effect_size_Rl") %>%
    dplyr::filter(grepl("CWM|FD", Predictor), !is.na(Effect_size_Rl)) %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::mutate(N_boot = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  dplyr::inner_join(res, rel,
                    by = c("Predictor", "Variable_gradient", "N_boot"))
}

boot_LUI  = join_responses_boot(d$Res_LUI,  d$Rel_LUI)
boot_SPEI = join_responses_boot(d$Res_SPEI, d$Rel_SPEI)

# per-draw correlations for panel b

boot_cors = function(boot_df, gradient_label) {
  boot_df %>%
    dplyr::group_by(Predictor, N_boot) %>%
    dplyr::summarise(
      r = cor(Effect_size_Rs, Effect_size_Rl,
              method = "spearman", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Gradient = gradient_label)
}

cors_all = dplyr::bind_rows(
  boot_cors(boot_LUI,  "LUI"),
  boot_cors(boot_SPEI, "SPEI")
)

# averaged per-window table for panel a

avg_pairs = function(boot_df) {
  boot_df %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::summarise(
      Effect_size_Rs = mean(Effect_size_Rs, na.rm = TRUE),
      Effect_size_Rl = mean(Effect_size_Rl, na.rm = TRUE),
      .groups = "drop"
    )
}
avg_LUI  = avg_pairs(boot_LUI)
avg_SPEI = avg_pairs(boot_SPEI)


# relabel predictors

relabel = function(x) {
  dplyr::recode(x,
                "CWM_PC1" = "CWM PC1: Economics (slow - fast)",
                "CWM_PC2" = "CWM PC2: Absorption (AMF - hairs)",
                "CWM_PC3" = "CWM PC3: Exploration (AD - SRL)",
                "FD_PC1"  = "FD PC1: Economics (slow - fast)",
                "FD_PC2"  = "FD PC2: Absorption (AMF - hairs)",
                "FD_PC3"  = "FD PC3: Exploration (AD - SRL)"
  )
}

cors_all = cors_all %>% dplyr::mutate(Predictor = relabel(Predictor))
avg_LUI  = avg_LUI  %>% dplyr::mutate(Predictor = relabel(Predictor))
avg_SPEI = avg_SPEI %>% dplyr::mutate(Predictor = relabel(Predictor))
boot_LUI = boot_LUI %>% dplyr::mutate(Predictor = relabel(Predictor))

predictor_colors = c(
  "CWM PC1: Economics (slow - fast)"      = "#3D405B",
  "CWM PC2: Absorption (AMF - hairs)"     = "#E07A5F",
  "CWM PC3: Exploration (AD - SRL)"  = "#81B29A",
  "FD PC1: Economics (slow - fast)"       = "#F2CC8F",
  "FD PC2: Absorption (AMF - hairs)"      = "#9A8C98",
  "FD PC3: Exploration (AD - SRL)"   = "#22223B"
)


#panela
d_A_avg = avg_LUI %>% dplyr::filter(Predictor %in% c("CWM PC3: Exploration (AD - SRL)",
                                                     "FD PC3: Exploration (AD - SRL)"))
# bootstrapped slopes fitted on the averaged curve, one per draw
d_A_slopes = boot_LUI %>%
  dplyr::filter(Predictor %in% c("CWM PC3: Exploration (AD - SRL)",
                                 "FD PC3: Exploration (AD - SRL)")) %>%
  dplyr::group_by(N_boot,Predictor) %>%
  dplyr::summarise(
    slope     = coef(lm(Effect_size_Rl ~ Effect_size_Rs))[2],
    intercept = coef(lm(Effect_size_Rl ~ Effect_size_Rs))[1],
    .groups   = "drop"
  )

p_A = ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_abline(data = d_A_slopes,
              aes(slope = slope, intercept = intercept,group=Predictor),
              color = "#B5EAD7", alpha = 0.15, linewidth = 0.3) +
  facet_wrap(.~Predictor,nrow=2,scales="free")+
  geom_point(data = d_A_avg,
             aes(x = Effect_size_Rs, y = Effect_size_Rl),
             color = "black", size = 2) +
  labs(x = "Effect size on drought resistance",
       y = "Effect size on drought resilience") +
  the_theme2

cors_all = cors_all %>%
  dplyr::mutate(
    PC         = sub("^\\S+\\s+(PC\\d).*$", "\\1", Predictor),
    Trait_type = sub("^(\\S+).*$", "\\1", Predictor)
  )%>%
  dplyr::mutate(., Gradient=recode_factor(Gradient,"LUI"="Along LUI gradient",
                                          "SPEI"="Along SPEI gradient"))

p_B = ggplot(cors_all,
             aes(x = r, fill = PC, color = PC)) +
  geom_density(alpha = 0.7, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(Gradient ~ Trait_type, scales = "fixed") +
  scale_fill_manual(values  = c("PC1" = "#C7CEEA",
                                "PC2" = "#FFDAC1",
                                "PC3" = "#B5EAD7"),
                    labels = c("PC1: Economics (slow - fast)",
                               "PC2: Absorption (AMF - hairs)",
                               "PC3: Exploration (AD - SRL)")) +
  scale_color_manual(values = c("PC1" = "#C7CEEA",
                                "PC2" = "#FFDAC1",
                                "PC3" = "#B5EAD7"),
                     labels = c("PC1: Economics (slow - fast)",
                                "PC2: Absorption (AMF - hairs)",
                                "PC3: Exploration (AD - SRL)")) +
  labs(x = "Correlation of effect size on resistance and resilience",
       y = "Density", fill = "", color = "") +
  the_theme2 +
  theme(    
    strip.text.x.top = element_text(colour = "white", size = 13),
    strip.text.y.right = element_text(colour = "white", size = 13),
    legend.position = "bottom",
    strip.background = element_rect(fill = "black", color = NA),
    panel.spacing = unit(0.6, "lines"))

ggsave("./Figures/Figure_5_tradeoff.pdf",
       ggarrange(p_A, p_B,ncol = 2, nrow = 1,
                 widths = c(1, 1.6),labels = c("a", "b")),
       width = 13, height = 7)


# --------------------------- Design Extended ----------------------------------------

d=readRDS("./Results/Moving_window_gradient_all_results_SPEI.rds")

p_CWM = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI06$Partial_residuals,
  d_Rl = d$SPEI06$Partial_residuals,negative_x = T,
  grepl_character = "CWM"
)
ggsave("./Figures/Figure_extended_SPEI_CWM.pdf",
       p_CWM+ylab("Effect size of mean community trait")+
         xlab("Drought intensity \n <- Moderate drought      Extreme drought ->"),
       width = 8,height = 4)



d = readRDS("./Results/Paired_results.rds")

pairs = list(
  LUI  = join_responses(d$Res_LUI,  d$Rel_LUI),
  SPEI = join_responses(d$Res_SPEI, d$Rel_SPEI)
)

pairs = lapply(pairs, relabel_predictors)

quad_colors = c("Trade-off (+ on resilience, - on resistance)"= "pink",
                "Trade-off (- on resilience, + on resistance)"="#DA6A6A",
  "Synergy (+ on both)"           = "lightgreen",
  "Synergy (- on both)"           = "forestgreen"
)
quad_LUI  = quadrant_tab(pairs$LUI)
quad_SPEI = quadrant_tab(pairs$SPEI)


p_tradeoff_LUI  = plot_tradeoff(pairs$LUI,  "Land-use intensity")
p_tradeoff_SPEI = plot_tradeoff(pairs$SPEI, "SPEI")

ggsave("./Figures/Extended_Tradeoff_Rs_Rl_LUI.pdf",
       p_tradeoff_LUI,  width = 10, height = 6)
ggsave("./Figures/Extended_Tradeoff_Rs_Rl_SPEI.pdf",
       p_tradeoff_SPEI, width = 10, height = 6)


ggsave("./Figures/Extended_tradeoff_Rs_Rl.pdf",
       ggarrange(plot_quad(quad_LUI, "Land use intensity (LUI)"),
                 plot_quad(quad_SPEI, "Drought intensity (SPEI)")+
                   theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()),
                 common.legend = T,legend = "bottom",widths = c(1,.6),
                 ncol=2), 
       width = 10, height = 5)



#NULL change effects SPEI

d = readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model_shuffle.rds")
bin_width = 0.03     

d_plot = d %>%
  dplyr::filter(grepl("CWM|FD", Predictor)) %>%
  Rename_PCA_axes(.,melt = T)%>%
  dplyr::mutate(
    Variable_gradient_bin = round(Variable_gradient / bin_width) * bin_width
  )

alpha_ = .3

p = ggplot(d_plot) +
  stat_summary(aes(
    x = Variable_gradient_bin,
    y = Effect_size,
    group = Predictor),fill = "grey",
    fun.data = function(x) {
      data.frame(ymin = quantile(x, (1 - 95/100) / 2),
                 ymax = quantile(x, 95/100 + (1 - 95/100) / 2))
    },
    geom = "ribbon", alpha = alpha_) +
  # stat_summary(aes(
  #   x = Variable_gradient_bin,
  #   y = Effect_size,
  #   group = Predictor),fill = "grey",
  #   fun.data = function(x) {
  #     data.frame(ymin = quantile(x, (1 - 95/100) / 2),
  #                ymax = quantile(x, 95/100 + (1 - 95/100) / 2))
  #   },
  #   geom = "smooth", alpha = alpha_) +
  geom_smooth(data=d_plot%>%
                dplyr::group_by(., Variable_gradient_bin,Predictor)%>%
                dplyr::summarise(., .groups = "keep",Effect_size=mean(Effect_size,na.rm=T)),
              aes(x = Variable_gradient_bin,
                  y = Effect_size),
              color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Predictor, scales = "free_y", nrow = 2) +
  labs(x = "Land-use intensity", y = "Effect size") +
  the_theme2+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "black"),
        strip.text       = element_text(colour = "white"))


#NULL change effects SPEI

d = readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model_shuffle_SPEI.rds")
bin_width = 0.03     

d_plot = d %>%
  dplyr::filter(grepl("CWM|FD", Predictor)) %>%
  Rename_PCA_axes(.,melt = T)%>%
  dplyr::mutate(
    Variable_gradient_bin = round(Variable_gradient / bin_width) * bin_width
  )

alpha_ = .3

p2 = ggplot(d_plot) +
  stat_summary(aes(
    x = Variable_gradient_bin,
    y = Effect_size,
    group = Predictor),fill = "grey",
    fun.data = function(x) {
      data.frame(ymin = quantile(x, (1 - 95/100) / 2),
                 ymax = quantile(x, 95/100 + (1 - 95/100) / 2))
    },
    geom = "ribbon", alpha = alpha_) +
  # stat_summary(aes(
  #   x = Variable_gradient_bin,
  #   y = Effect_size,
  #   group = Predictor),fill = "grey",
  #   fun.data = function(x) {
  #     data.frame(ymin = quantile(x, (1 - 95/100) / 2),
  #                ymax = quantile(x, 95/100 + (1 - 95/100) / 2))
  #   },
  #   geom = "smooth", alpha = alpha_) +
  geom_smooth(data=d_plot%>%
                dplyr::group_by(., Variable_gradient_bin,Predictor)%>%
                dplyr::summarise(., .groups = "keep",Effect_size=mean(Effect_size,na.rm=T)),
              aes(x = Variable_gradient_bin,
                  y = Effect_size),
              color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Predictor, scales = "free_y", nrow = 2) +
  labs(x = "Drought intensity \n <- Moderate drought      Extreme drought ->", y = "Effect size") +
  the_theme2+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "black"),
        strip.text       = element_text(colour = "white"))

ggsave("./Figures/Moving_window_shuffle_bootstraps_binned.pdf",
       ggarrange(p+ylim(-.5,.5),
                 p2+ylim(-.5,.5),nrow=2,labels = letters[1:2]),width = 11, height = 12)


#NULL Distributions  correlaitons 
obs = readRDS("./Results/Paired_results.rds")

obs_from_paired = function(res_obj, rel_obj) {
  
  d_Rs = res_obj$Partial_residuals %>%
    dplyr::select(-Stability_var) %>%
    tidyr::pivot_longer(
      cols      = -Variable_gradient,
      names_to  = "Predictor",
      values_to = "Effect_size_Rs"
    ) %>%
    dplyr::filter(!is.na(Effect_size_Rs)) %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::summarise(
      Effect_size_Rs = mean(Effect_size_Rs, na.rm = TRUE),
      .groups        = "drop"
    )
  
  d_Rl = rel_obj$Partial_residuals %>%
    dplyr::select(-Stability_var) %>%
    tidyr::pivot_longer(
      cols      = -Variable_gradient,
      names_to  = "Predictor",
      values_to = "Effect_size_Rl"
    ) %>%
    dplyr::filter(!is.na(Effect_size_Rl)) %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::summarise(
      Effect_size_Rl = mean(Effect_size_Rl, na.rm = TRUE),
      .groups        = "drop"
    )
  
  dplyr::inner_join(d_Rs, d_Rl,
                    by = c("Predictor", "Variable_gradient")) %>%
    dplyr::filter(grepl("CWM|FD|Species|PSV|MNTD", Predictor)) %>%
    dplyr::group_by(Predictor) %>%
    dplyr::summarise(
      slope_obs = coef(lm(Effect_size_Rl ~ Effect_size_Rs))[2],
      r_obs     = cor(Effect_size_Rs, Effect_size_Rl,
                      method = "spearman", use = "complete.obs"),
      .groups   = "drop"
    )
}

obs_LUI  = obs_from_paired(obs$Res_LUI,  obs$Rel_LUI)  %>% dplyr::mutate(Gradient = "Along LUI gradient")
obs_all = dplyr::bind_rows(obs_LUI)
null_LUI  = readRDS("./Results/Null_paired_LUI.rds")  %>% dplyr::mutate(Gradient = "Along LUI gradient")
null_all = dplyr::bind_rows(null_LUI)


p_vals = null_all %>%
  dplyr::left_join(obs_all, by = c("Gradient", "Predictor")) %>%
  dplyr::group_by(Gradient, Predictor) %>%
  dplyr::summarise(
    r_null_mean = mean(r, na.rm = TRUE),
    r_null_lo   = quantile(r, 0.025, na.rm = TRUE),
    r_null_hi   = quantile(r, 0.975, na.rm = TRUE),
    r_obs       = dplyr::first(r_obs),
    p_two_sided = {
      med = median(r, na.rm = TRUE)
      mean(abs(r - med) >= abs(dplyr::first(r_obs) - med), na.rm = TRUE)
    },
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    label = sprintf("%s\np = %.3f", Predictor, p_two_sided)
  )

null_plot = null_all %>%
  dplyr::left_join(p_vals %>% dplyr::select(Gradient, Predictor, label),
                   by = c("Gradient", "Predictor"))


null_plot = null_plot %>% dplyr::mutate(label = recode_factor(label,
                                                                  "CWM_PC1\np = 0.040" = "CWM PC1: Economics\np = 0.040",
                                                                  "CWM_PC2\np = 0.030" = "CWM PC2: Absorption\np = 0.030",
                                                                  "CWM_PC3\np = 0.450" = "CWM PC3: Exploration\np = 0.450",
                                                                  "FD_PC1\np = 0.340"  = "FD PC1: Economics\np = 0.340",
                                                                  "FD_PC2\np = 0.310"  = "FD PC2: Absorption\np = 0.310",
                                                                  "FD_PC3\np = 0.030"  = "FD PC3: Exploration\np = 0.030",
                                                              ))%>%
  dplyr::filter(,label %!in% c("MNTD_weighted\np = 0.470","PSV\np = 0.750","Species_richness\np = 0.190"))

p_vals    = p_vals    %>% dplyr::mutate(label = recode_factor(label,
                                                                  "CWM_PC1\np = 0.040" = "CWM PC1: Economics\np = 0.040",
                                                                  "CWM_PC2\np = 0.030" = "CWM PC2: Absorption\np = 0.030",
                                                                  "CWM_PC3\np = 0.450" = "CWM PC3: Exploration\np = 0.450",
                                                                  "FD_PC1\np = 0.340"  = "FD PC1: Economics\np = 0.340",
                                                                  "FD_PC2\np = 0.310"  = "FD PC2: Absorption\np = 0.310",
                                                                  "FD_PC3\np = 0.030"  = "FD PC3: Exploration\np = 0.030",
                                                              ))%>%
  dplyr::filter(,label %!in% c("MNTD_weighted\np = 0.470","PSV\np = 0.750","Species_richness\np = 0.190"))

p_null = ggplot(null_plot %>% dplyr::filter(Gradient == "Along LUI gradient"),
                aes(x = r)) +
  geom_density(fill = "grey75", color = "grey40",
               alpha = 0.7, linewidth = 0.4) +
  geom_vline(data = p_vals %>% dplyr::filter(Gradient == "Along LUI gradient"),
             aes(xintercept = r_obs),
             color = "#D7263D", linewidth = 0.9) +
  geom_vline(data = p_vals %>% dplyr::filter(Gradient == "Along LUI gradient"),
             aes(xintercept = r_null_lo),
             color = "grey30", linetype = "dashed", linewidth = 0.4) +
  geom_vline(data = p_vals %>% dplyr::filter(Gradient == "Along LUI gradient"),
             aes(xintercept = r_null_hi),
             color = "grey30", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~ label, ncol = 3, scales = "fixed") +
  labs(x = "Spearman correlation of effects between resistance and resilience",
       y = "Density") +
  the_theme2 +
  theme(strip.background = element_rect(fill = "black"),
        strip.text       = element_text(colour = "white", size = 8))

ggsave("./Figures/Null_distributions_correlation_effects.pdf",p_null,width = 9,height = 5)

# --------------------------- Supplementary ------------------------------------
## >> Map of the studied sites ----

#import kml of vegetation plots
Kml_obs.dir = list.files(path = "./Data/Geo/", pattern = ".kml", full.names = T)

Kml_obs = lapply(Kml_obs.dir, st_read)

names(Kml_obs) = c("A", "H", "S")

str(Kml_obs$A)

Kml_obs = lapply(Kml_obs,function(.) {
  .$MyDescr = ifelse(.$Description == "VIP - Forest" | .$Description == "EP - Forest", "Forest", "Grassland")
  .
})


Kml_obs.gr = lapply(Kml_obs,function(.) {
  df.gr = .[which(.$MyDescr == "Grassland"), ]
  return(df.gr)
}) 

Centroids.gr = lapply(Kml_obs.gr,function(.) {
  Pl_obs = st_transform(., crs = 25832)
  Cntrd = st_centroid(st_union(Pl_obs))
  Cntrd.geo = st_transform(Cntrd, crs = 4326)
})

Centroids.gr = c(Centroids.gr$A, Centroids.gr$H, Centroids.gr$S)



##----------------------------------------------------------------Map of the study area

#get DEM for the 3 regions
Bbox_ALB = st_as_sf(data.frame(ID = seq_len(4), x = c(9.2, 9.2, 9.6, 9.6), y = c(48.35, 48.525, 48.35, 48.525)), coords = c("x", "y"))
Bbox_HAI = st_as_sf(data.frame(ID = seq_len(4), x = c(10.12, 10.12, 10.8, 10.8), y = c(50.92, 51.36, 50.92, 51.36)), coords = c("x", "y"))
Bbox_SCH = st_as_sf(data.frame(ID = seq_len(4), x = c(13.55, 13.55, 14.09, 14.09), y = c(52.83, 53.18, 52.83, 53.18)), coords = c("x", "y"))

st_crs(Bbox_ALB) = st_crs(Bbox_HAI) = st_crs(Bbox_SCH) = 4326

DEM_ALB = elevatr::get_elev_raster(locations = Bbox_ALB, z = 10, clip = "locations")
DEM_ALB = as.data.frame(DEM_ALB, na.rm = T, xy = T)

DEM_HAI = elevatr::get_elev_raster(locations = Bbox_HAI, z = 10, clip = "locations")
DEM_HAI = as.data.frame(DEM_HAI, na.rm = T, xy = T)

DEM_SCH = elevatr::get_elev_raster(locations = Bbox_SCH, z = 10, clip = "locations")
DEM_SCH = as.data.frame(DEM_SCH, na.rm = T, xy = T)

#rename elevation column
colnames(DEM_ALB)[3] = colnames(DEM_HAI)[3] = colnames(DEM_SCH)[3] = "Elevation"


#ALB_map, HAI_map and SCH_map updated to use terrain.colors as palette for elevation on 14/11/2023
p1=ggplot() +
  geom_raster(data = DEM_ALB, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_viridis_c(option = "F")+
  geom_sf(data = Kml_obs.gr$A, col = "black", size = 3) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.1, "cm"), height = unit(0.5, "cm"), text_cex = 1) +
  ylab(NULL) + xlab(NULL) + labs(fill="Altitude (m)")+
  the_theme2+ggtitle("South-West")

p2 = ggplot() +
  geom_raster(data = DEM_HAI, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_viridis_c(option = "F")+
  geom_sf(data = Kml_obs.gr$H, col = "black", size = 3) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.1, "cm"), height = unit(0.5, "cm"), text_cex = 1) +
  ylab(NULL) + xlab(NULL) + labs(fill="Altitude (m)")+
  the_theme2+ggtitle("Center")

p3 =  ggplot() +
  geom_raster(data = DEM_SCH, aes(x = x, y = y, fill = Elevation)) +
  scale_fill_viridis_c(option = "F")+
  geom_sf(data = Kml_obs.gr$S, col = "black", size = 3) +
  annotation_scale(pad_x = unit(1, "cm"), pad_y = unit(0.1, "cm"), height = unit(0.5, "cm"), text_cex = 1) +
  ylab(NULL) + xlab(NULL) + labs(fill="Altitude (m)")+
  the_theme2+ggtitle("North-East")

p1_tot=ggarrange(p1,p2,p3,ncol=3,labels = c("a","",""),widths = c(1,1.3,1.3),align = "hv")

clim_data=read.table("./Data/Geo/SM_10_Ta_200_precipitation_radolan_bc8f4bdf98468e3c/plots.csv",sep=",",header = T)%>%
  dplyr::mutate(., 
                Year = as.numeric(substr(datetime, 1, 4)),
                Month = as.numeric(substr(datetime, 6, 7))
  )%>%
  dplyr::group_by(., plotID,Year)%>%
  dplyr::summarise(.,.groups = "drop",MT=mean(Ta_200,na.rm=T),MP=sum(precipitation_radolan,na.rm=T))%>%
  dplyr::group_by(., plotID)%>%
  dplyr::summarise(.,.groups = "drop",MT=mean(MT,na.rm=T),MP=mean(MP))%>%
  merge(., read.table("./Data/Geo/Geographical_climatic.csv",sep=";")%>%
          dplyr::select(., site_name,Elevation),by.x="plotID",by.y="site_name")


p2=ggplot(clim_data %>% 
            melt(., id.vars = c("plotID")) %>%
            dplyr::mutate(., 
                          plotID = gsub("[0-9]", "", plotID),
                          variable = dplyr::recode(variable,
                                                   "MT" = "Mean annual temperature",
                                                   "MP" = "Mean annual precipitation"
                          ))%>%
            dplyr::mutate(., 
                          plotID = dplyr::recode(plotID,"AEG"="South-West","SEG"="North-East","HEG"="Center"))) +
  geom_histogram(aes(x = value, fill = plotID)) +
  the_theme2 +
  facet_wrap(. ~ variable, 
             scales = "free",
             strip.position = "bottom") +
  scale_fill_manual(values = c("#C7CEEA", "#FFDAC1", "#B5EAD7")) +
  labs(x = "", y = "Count", fill = "") +
  theme(
    strip.placement = "outside"  )

ggsave("./Figures/Maps_sites.pdf",ggarrange(p1_tot,p2,nrow=2,labels=c("","b")),width = 12,height = 10)

## >> Distribution traits ----

trait_sp=read.table(paste0("./Data/Global_trait_dataset_Hairphae.csv"),sep=";")%>%
  dplyr::select(.,-PC1,-PC2,-PC3)
colnames(trait_sp)=c("SRL","RTD","AD","RN","RHL",'AMF','RHI',"SM","SLA","Height","LDMC","LNC","LPC")

p=Plot_distributions_variables(trait_sp)+labs(x="Trait value",y="Count")

ggsave("./Figures/Trait_distribution.pdf",p,width = 8,height = 5)

## >> PCA traits ----

global_traits=read.table("./Data/Global_trait_dataset.csv",sep=";")%>%
  Closer_to_normal_traits(.)%>%
  dplyr::select(., all_of(c(colnames(.)[grep("_HP",colnames(.))],"Seed_mass",
                            # "Myco_intensity",
                            "SLA_all","Height","LDMC","LeafN","LeafP")))%>%
  drop_na(.)

colnames(global_traits)=c("SRL","RTD","AD","RN","RHL",'AMF','RHI',"SM","SLA","Height","LDMC","LNC","LPC")

paranResult = paran::paran(global_traits,
                           iterations = 5000, centile = 95, quietly = TRUE, status = FALSE)

paranResult$Retained
fitAll = psych::principal(r = global_traits,nfactors = 3,rotate = "varimax")
p1=plotPCA(fitAll$scores,
           fitAll$loadings,
           fitAll$Vaccounted,annotateFactor = 3.2,
           xIndex = 1,yIndex = 2,
           xLim = c(-3,3),yLim = c(-3,3))

p2=plotPCA(fitAll$scores,
           fitAll$loadings,
           fitAll$Vaccounted,annotateFactor = 3.2,
           xIndex = 1,yIndex = 3,
           xLim = c(-3,3),yLim = c(-3,3))

ggsave("./Figures/PCA_traits_above_below_Hairphae.pdf",ggarrange(p1,p2),width = 10,height = 4)

## >> Correlation functional diversity ----

d_mod=readRDS("./Data/Merged_datasets_Hairphae.rds")[[1]]%>%
  dplyr::filter(., 
                Trait_representativity>70)

p_pc1 = Plot_correlation_variables(d_mod[, c("FDis_PC1", "FD_PC1", "RaoQ_PC1")])
p_pc2 = Plot_correlation_variables(d_mod[, c("FDis_PC2", "FD_PC2", "RaoQ_PC2")])
p_pc3 = Plot_correlation_variables(d_mod[, c("FDis_PC3", "FD_PC3", "RaoQ_PC3")])

ggarrange(p_pc1,p_pc2,p_pc3,ncol=3)

## >> Correlation predictors LUI ----

d=readRDS("./Data/Merged_datasets_Hairphae.rds")[[1]]%>%
  dplyr::filter(., 
                Trait_representativity>70)%>%
  reshape2::melt(., measure.vars = c("CWM_PC1","CWM_PC2","CWM_PC3",
                                     paste0("FD_PC",1:3),"Species_richness",
                                     "PSV","MNTD_weighted"))%>%
  dplyr::mutate(
    variable = dplyr::recode(variable,
                              "CWM_PC1" = "CWM PC1: Economics (slow - fast)",
                              "CWM_PC2" = "CWM PC2: Absorption (AMF - hairs)",
                              "CWM_PC3" = "CWM PC3: Exploration (AD - SRL)",
                              "FD_PC1"  = "FD PC1: Economics (slow - fast)",
                              "FD_PC2"  = "FD PC2: Absorption (AMF - hairs)",
                              "FD_PC3"  = "FD PC3: Exploration (AD - SRL)"
    )
  )


variables_list = levels(d$variable)
plot_list = list()

for(i in seq_along(variables_list)) {
  var_name = variables_list[i]
  
  d_subset = d %>% filter(variable == var_name)
  
  significance = coef(summary(lm(value ~ LUI_global, data = d_subset)))[2, 4]
  linetype_val = ifelse(significance < 0.05, "solid", "dashed")
  
  assign(paste0("p_", i), 
         ggplot(d_subset, aes(x = LUI_global, y = value)) +
           geom_point(fill = "white", color = "black", shape = 21) +
           geom_smooth(method = "lm", formula = y ~ poly(x, 1), 
                       color = "red", fill = "pink", linetype = linetype_val) +
           the_theme2 +
           labs(x = "LUI", y = var_name) +
           theme(plot.title = element_text(hjust = 0.5, size = 10)))
  
  plot_list[[i]] = get(paste0("p_", i))
}


p_all = ggarrange(plotlist = plot_list, ncol = 3, nrow = 3, align = "hv")

ggsave("./Figures/Correlation_LUI_predictors.pdf",p_all,width = 12,height = 10)



## >> Timeseries ----

Merged_data=readRDS("./Data/Vegetation_SPEI.rds")
pdf("./Figures/All_vegetation_series.pdf",width = 20,height = 5)
par(mfrow=c(1,3))

data_k_detrended=Merged_data$spei01%>%
  dplyr::mutate(., Total_biomass=as.numeric(Total_biomass))%>%
  dplyr::group_by(Site) %>%
  dplyr::arrange(Year, .by_group = TRUE) %>%
  dplyr::mutate(time = dplyr::row_number(),  
                biomass_detrended_loess = Total_biomass - predict(loess(Total_biomass ~ time, na.action = na.exclude)),
                biomass_detrended_lm = residuals(lm(Total_biomass ~ time, na.action = na.exclude)))

for (k in unique(Merged_data$spei01$Site)){
  
  plot(as.numeric(data_k_detrended$Year[which(data_k_detrended$Site==k)]),
       as.numeric(data_k_detrended$Total_biomass[which(data_k_detrended$Site==k)]),
       type="b",xlab="Year",ylab="Total biomass")
  plot(as.numeric(data_k_detrended$Year[which(data_k_detrended$Site==k)]),
       as.numeric(data_k_detrended$biomass_detrended_loess[which(data_k_detrended$Site==k)]),
       type="b",xlab="Year",ylab="Total biomass (detrended loess)")
  plot(as.numeric(data_k_detrended$Year[which(data_k_detrended$Site==k)]),
       as.numeric(data_k_detrended$biomass_detrended_lm[which(data_k_detrended$Site==k)]),
       type="b",xlab="Year",ylab="Total biomass (detrended lm)")
  mtext(k)
}
dev.off()





Merged_data=readRDS("./Data/Vegetation_SPEI.rds")
pdf("./Figures/All_vegetation_series_Cover.pdf",width = 20,height = 5)
par(mfrow=c(1,3))

data_k_detrended=Merged_data$spei01%>%
  dplyr::mutate(., Total_biomass=as.numeric(Cover))%>%
  dplyr::group_by(Site) %>%
  dplyr::arrange(Year, .by_group = TRUE) %>%
  dplyr::mutate(time = dplyr::row_number(),  
                biomass_detrended_loess = Cover - predict(loess(Cover ~ time, na.action = na.exclude)),
                biomass_detrended_lm = residuals(lm(Cover ~ time, na.action = na.exclude)))

for (k in unique(Merged_data$spei01$Site)){
  
  plot(as.numeric(data_k_detrended$Year[which(data_k_detrended$Site==k)]),
       as.numeric(data_k_detrended$Cover[which(data_k_detrended$Site==k)]),
       type="b",xlab="Year",ylab="Total biomass")
  plot(as.numeric(data_k_detrended$Year[which(data_k_detrended$Site==k)]),
       as.numeric(data_k_detrended$biomass_detrended_loess[which(data_k_detrended$Site==k)]),
       type="b",xlab="Year",ylab="Total biomass (detrended loess)")
  plot(as.numeric(data_k_detrended$Year[which(data_k_detrended$Site==k)]),
       as.numeric(data_k_detrended$biomass_detrended_lm[which(data_k_detrended$Site==k)]),
       type="b",xlab="Year",ylab="Total biomass (detrended lm)")
  mtext(k)
}
dev.off()


## >> Correlations predictors ----

d_mod=readRDS("./Data/Merged_datasets_Hairphae.rds")[["spei06"]]%>%
  dplyr::filter(., 
                Trait_representativity>70,
                Resistance_Isbell>0,
                Drought_category %in% c("Moderate drought","Extreme drought"))%>%
  dplyr::mutate(.,Species_richness=sqrt(Species_richness))%>%
  dplyr::select(., 
                Species_richness,
                FD_PC1,FD_PC2,FD_PC3,
                CWM_PC1,CWM_PC2,CWM_PC3,
                MNTD_weighted,PSV,Temperature,Soil_moisture,Precipitation,
                Species_richness,LUI_global)%>%
  drop_na(.)

p=Plot_correlation_variables(d_mod)+the_theme2+theme(axis.text.x = element_text(angle=60,hjust=1))

ggsave("./Figures/Correlation_predictors.pdf",p,width = 6,height = 6)


d_mod=readRDS("./Data/Merged_datasets_Hairphae.rds")[["spei06"]]%>%
  dplyr::filter(., 
                Trait_representativity>70,
                Resistance_Isbell>0,
                Drought_category %in% c("Moderate drought","Extreme drought"))%>%
  dplyr::mutate(.,Species_richness=sqrt(Species_richness))%>%
  dplyr::select(., 
                Species_richness,
                FD_PC1,FD_PC2,FD_PC3,
                CWM_PC1,CWM_PC2,CWM_PC3,
                MNTD_weighted,PSV,
                Species_richness,LUI_global)%>%
  drop_na(.)

p=Plot_distributions_variables(d_mod)+labs(x="Value",y="Count")+facet_wrap(.~variable,scales="free",ncol=3)

ggsave("./Figures/Distribution_predictors.pdf",p,width = 7,height = 7)

readRDS("./Data/Merged_datasets_Underplot.rds")[["spei06"]]%>%
  dplyr::filter(., 
                Trait_representativity>70,
                Resistance_Isbell>0,
                Drought_category %in% c("Moderate drought","Extreme drought"))%>%
  dplyr::mutate(.,Species_richness=sqrt(Species_richness))%>%
  dplyr::select(., 
                SPEI_value,
                LUI_global,
                LUI_regional,
                Logbiomass_t_median,
                Logbiomass_t_minus_1,
                Resistance_Isbell,
                Resilience_Isbell,
                Total_biomass,
                Species_richness)%>%
  drop_na(.)%>%
  Plot_correlation_variables(.)


## >> Correlation stability, LUI, Clim ----


data_hairphae = readRDS("./Data/Merged_datasets_Hairphae.rds")[["spei06"]]
d_res = data_hairphae %>%
  dplyr::filter(., 
                Trait_representativity > 70,
                Resistance_Isbell > 0,
                Drought_category %in% c("Moderate drought", "Extreme drought")) %>%
  dplyr::select(Site, Year, Resistance_Isbell,Precipitation,Temperature,Soil_moisture,LUI_global) %>%
  dplyr::mutate(., Resistance_Isbell = bcPower(Resistance_Isbell, -0.7))

d_resil = data_hairphae %>%
  dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
  dplyr::filter(., 
                Trait_representativity > 70,
                Resilience_Isbell > 0,
                Drought_category %in% c("Moderate drought", "Extreme drought"),
                !is.na(Resilience_Isbell),
                is.na(Recovery_context)) %>%
  dplyr::select(Site, Year, Resilience_Isbell,Precipitation,Temperature,Soil_moisture,LUI_global) %>%
  dplyr::mutate(., Resilience_Isbell = log(Resilience_Isbell))
d_all=d_resil%>%
  merge(., d_res, by = "Site", all.x = TRUE, suffixes = c("", ".res"))

d_all = d_all %>%
  dplyr::rename(
    Resilience = Resilience_Isbell,
    Resistance = Resistance_Isbell
  )


d_all_melt = d_all %>%
  melt(., measure.vars = c("Resistance", "Resilience"),
       id.vars = c("Site", "Year", "LUI_global", "Temperature", "Precipitation", "Soil_moisture"))

d_all_melt = d_all_melt %>%
  dplyr::filter(!is.na(value))

nice_names = c("Resistance" = "Drought resistance","Resilience" = "Drought resilience")
drivers = c("LUI_global", "Temperature", "Precipitation", "Soil_moisture")
nice_driver_names = c(
  "LUI_global" = "Land-use intensity",
  "Temperature" = "Temperature",
  "Precipitation" = "Precipitation",
  "Soil_moisture" = "Soil moisture")

variables_list = unique(d_all_melt$variable)
get_lmer_p =function(model, driver) {
  p_value =coef(summary(model))[driver, "Pr(>|t|)"]
  return(p_value)
}

plot_list =list()

for (d in seq_along(drivers)) {
  driver =drivers[d]
  driver_nice =nice_driver_names[driver]
  
  for (v in seq_along(variables_list)) {
    var_name =variables_list[v]
    
    d_subset =d_all_melt %>%
      dplyr::filter(variable == var_name) %>%
      dplyr::filter(!is.na(.data[[driver]]))
    
    if (nrow(d_subset) < 3) next
    
    
    d_subset[[driver]] =scale(d_subset[[driver]])[, 1]
    d_subset[["value"]] =scale(d_subset[["value"]])[, 1]
    
    
    model =lmer(as.formula(paste0("value ~ ", driver, " + (1|Site)+ (1|Year)")), 
                data = d_subset)
    p_value =coef(summary(model))[2, "Pr(>|t|)"]
    
    linetype_val =ifelse(p_value < 0.05, "solid", "dashed")
    
    p =ggplot(d_subset, aes_string(x = driver, y = "value")) +
      geom_point(color="black",fill = "white", shape = 21, alpha = 0.6, size = 2, show.legend = FALSE) +
      geom_smooth(method = "lm", formula = y ~ poly(x, 1), 
                  color = "red", fill = "pink", linetype = linetype_val, alpha = 0.3) +
      the_theme2 +
      labs(x = driver_nice, y = nice_names[var_name]) +
      the_theme2
    
    plot_list[[length(plot_list) + 1]] =p
  }
}

p_all =ggarrange(plotlist = plot_list, 
                 ncol = length(variables_list), 
                 nrow = length(drivers), 
                 align = "hv")

ggsave("./Figures/Correlation_stability_clim_land_use.pdf", p_all, width = 7, height = 8)



d_cor = readRDS("./Data/Merged_datasets_Hairphae.rds")[["spei06"]] %>%
  dplyr::filter(., 
                Trait_representativity > 70,
                Resistance_Isbell > 0,
                Drought_category %in% c("Moderate drought", "Extreme drought")) %>%
  dplyr::select(., LUI_global, Temperature, Precipitation, Soil_moisture) %>%
  drop_na(.)

p = Plot_correlation_variables(d_cor) + the_theme2 + theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("./Figures/Correlation_drivers.pdf", p, width = 5, height = 5)

## >> Spatial and temporal autocorrelation models ----

d=readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model.rds")
temporal=d$SPEI06$Temporal_autocorrelation%>%
  dplyr::mutate(., Stability_var=recode_factor(Stability_var,
                                               "Resilience_Isbell_abs"="Resilience",
                                               "Resistance_Isbell"="Resistance"))%>%
  dplyr::filter(., Stability_var %in% c("Resistance","Resilience"))
spatial=d$SPEI06$Spatial_autocorrelation%>%
  dplyr::mutate(., Stability_var=recode_factor(Stability_var,
                                               "Resilience_Isbell_abs"="Resilience",
                                               "Resistance_Isbell"="Resistance"))%>%
  dplyr::filter(., Stability_var %in% c("Resistance","Resilience"))

make_autocorr_plot = function(stability_var_name, spatial_df, temporal_df, 
                              scale_factor = .15, scale_factor2 = 3) {
  
  sp = spatial_df %>% dplyr::filter(Stability_var == stability_var_name)
  te = temporal_df %>% dplyr::filter(Stability_var == stability_var_name)
  
  p_spatial = ggplot(sp, aes(x = Variable_gradient)) +
    geom_line(aes(y = Statistic), linewidth = .5, color = "grey") +
    geom_line(aes(y = Pvalue * scale_factor), linewidth = .5, color = "grey") +
    geom_point(aes(y = Pvalue * scale_factor, color = Pvalue < .05), size = 3) +
    geom_point(aes(y = Statistic), size = 3, color = "#767FC7") +
    scale_y_continuous(name = "Moran's I statistic",
                       sec.axis = sec_axis(~ . * scale_factor, name = paste0("P-value (scaled factor = ",scale_factor,")"))) +
    geom_hline(yintercept = .05 * scale_factor, color = "#BA70BB") +
    scale_color_manual(values = c("#BA70BB", "red")) +
    labs(x = "Land-use intensity") +
    ggtitle("Spatial autocorrelation") +
    the_theme2 + guides(color = "none", fill = "none") +
    theme(axis.title.y.right = element_text(color = "#BA70BB"))
  
  p_temporal = ggplot(te, aes(x = Variable_gradient)) +
    geom_line(aes(y = Statistic), linewidth = .5, color = "grey") +
    geom_line(aes(y = Pvalue * scale_factor2), linewidth = .5, color = "grey") +
    geom_point(aes(y = Pvalue * scale_factor2, color = Pvalue < .05), size = 3) +
    geom_point(aes(y = Statistic), size = 3, color = "#767FC7") +
    scale_y_continuous(name = "Durbin-Watson statistic",
                       sec.axis = sec_axis(~ . * scale_factor, name = paste0("P-value (scaled factor = ",scale_factor2,")"))) +
    geom_hline(yintercept = .05 * scale_factor2, color = "#BA70BB") +
    scale_color_manual(values = c("#BA70BB", "red")) +
    labs(x = "Land-use intensity") +
    ggtitle("Temporal autocorrelation") +
    the_theme2 + guides(color = "none", fill = "none") +
    theme(axis.title.y.right = element_text(color = "#BA70BB"))
  
  return(ggarrange(p_spatial, p_temporal, ncol = 2))
}

p1 = make_autocorr_plot("Resistance", spatial, temporal,scale_factor2 = 1.5)
p2 = make_autocorr_plot("Resilience", spatial, temporal,scale_factor2 = 2)

p_tot = ggarrange(p1, p2, nrow = 2, labels = letters[1:2])

ggsave("./Figures/Test_autocorrelation_space_time.pdf", p_tot, width = 10, height = 7)



## >> Fertilization, Mowing, Grazing ----

d=readRDS("./Results/Moving_window_gradient_all_results_fertilization.rds")

Plot_MV_binary_by_axis = function(d_res, d_rel,
                                  grepl_character = "CWM",
                                  CI_inner = 90,
                                  CI_outer = 95,
                                  gradientname = "fertilisation",
                                  color_Rs = "#22223B",
                                  color_Rl = "#EC7692",
                                  alpha_non_sig = 0.25) {
  
  prep = function(df, resp_label, stab) {
    df %>%
      dplyr::filter(Stability_var == stab) %>%
      dplyr::select(-Stability_var) %>%
      tidyr::pivot_longer(-Variable_gradient,
                          names_to  = "Predictor",
                          values_to = "Effect_size") %>%
      dplyr::filter(grepl(grepl_character, Predictor),
                    !is.na(Effect_size)) %>%
      dplyr::mutate(
        Axis = dplyr::case_when(
          grepl("PC1", Predictor) ~ "PC1: Economics (slow - fast)",
          grepl("PC2", Predictor) ~ "PC2: Absorption (AMF - hairs)",
          grepl("PC3", Predictor) ~ "PC3: Exploration (AD - SRL)",
          TRUE ~ "Other"
        ),
        Response = resp_label,
        Variable_gradient = factor(Variable_gradient,
                                   levels = c("NO", "YES"),
                                   labels = c(paste0("No ",   gradientname),
                                              paste0("With ", gradientname)))
      )
  }
  
  d_plot = dplyr::bind_rows(
    prep(d_res, "Drought resistance", "Resistance_Isbell"),
    prep(d_rel, "Drought resilience", "Resilience_Isbell_abs")
  ) %>%
    dplyr::filter(Axis != "Other") %>%
    dplyr::mutate(Response = factor(Response,
                                    levels = c("Drought resistance",
                                               "Drought resilience")))
  
  d_sum = d_plot %>%
    dplyr::group_by(Axis, Response, Variable_gradient) %>%
    dplyr::summarise(
      Stat   = median(Effect_size, na.rm = TRUE),
      q1_in  = quantile(Effect_size, (1 - CI_inner/100)/2, na.rm = TRUE),
      q3_in  = quantile(Effect_size, (1 - (1 - CI_inner/100)/2), na.rm = TRUE),
      q1_out = quantile(Effect_size, (1 - CI_outer/100)/2, na.rm = TRUE),
      q3_out = quantile(Effect_size, (1 - (1 - CI_outer/100)/2), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(significant = (q1_out > 0) | (q3_out < 0))
  
  pos = position_dodge(width = 0.6)
  
  p = ggplot(d_sum, aes(x = Variable_gradient, y = Stat,
                        color = Response, fill = Response,
                        alpha = significant,group=Response)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_linerange(aes(ymin = q1_out, ymax = q3_out),
                   lwd = .5, position = pos) +
    geom_linerange(aes(ymin = q1_in, ymax = q3_in),
                   lwd = 1.5, position = pos) +
    geom_point(shape = 21, color = "black", fill = "white",
               size = 4, position = pos) +
    scale_alpha_manual(values = c(`FALSE` = alpha_non_sig, `TRUE` = 1),
                       guide = "none") +
    scale_color_manual(values = c("Drought resistance" = color_Rs,
                                  "Drought resilience" = color_Rl)) +
    scale_fill_manual(values  = c("Drought resistance" = color_Rs,
                                  "Drought resilience" = color_Rl)) +
    facet_wrap(~Axis, nrow = 1, scales = "fixed") +
    labs(x = "", y = "Effect size", color = "", fill = "") +
    the_theme2 +
    theme(strip.background = element_rect(fill = "black"),
          strip.text       = element_text(colour = "white"),
          legend.position  = "bottom")  
  return(p)
}
d_fert = readRDS("./Results/Moving_window_gradient_all_results_fertilization.rds")

p_fert_CWM = Plot_MV_binary_by_axis(
  d_res = d_fert$SPEI06$Partial_residuals,
  d_rel = d_fert$SPEI06$Partial_residuals,
  grepl_character = "CWM",
  gradientname    = "fertilisation"
)

p_fert_FD = Plot_MV_binary_by_axis(
  d_res = d_fert$SPEI06$Partial_residuals,
  d_rel = d_fert$SPEI06$Partial_residuals,
  grepl_character = "FD",
  gradientname    = "fertilisation"
)

# Mowing
d_mow = readRDS("./Results/Moving_window_gradient_all_results_mowing.rds")

p_mow_CWM = Plot_MV_binary_by_axis(
  d_res = d_mow$SPEI06$Partial_residuals,
  d_rel = d_mow$SPEI06$Partial_residuals,
  grepl_character = "CWM",
  gradientname    = "mowing"
)

p_mow_FD = Plot_MV_binary_by_axis(
  d_res = d_mow$SPEI06$Partial_residuals,
  d_rel = d_mow$SPEI06$Partial_residuals,
  grepl_character = "FD",
  gradientname    = "mowing"
)

# Grazing
d_graz = readRDS("./Results/Moving_window_gradient_all_results_grazing.rds")

p_graz_CWM = Plot_MV_binary_by_axis(
  d_res = d_graz$SPEI06$Partial_residuals,
  d_rel = d_graz$SPEI06$Partial_residuals,
  grepl_character = "CWM",
  gradientname    = "grazing"
)

p_graz_FD = Plot_MV_binary_by_axis(
  d_res = d_graz$SPEI06$Partial_residuals,
  d_rel = d_graz$SPEI06$Partial_residuals,
  grepl_character = "FD",
  gradientname    = "grazing"
)
ggsave("./Figures/Trait_effects_resistance_fertilization.pdf",
       ggarrange(p_fert_CWM, p_fert_FD,nrow = 2, ncol = 1,common.legend = TRUE, legend = "bottom",
         labels = c("a", "b")),width = 9,height = 7)



ggsave("./Figures/Trait_effects_resistance_mowing.pdf",
       ggarrange(p_mow_CWM, p_mow_FD,nrow = 2, ncol = 1,
         common.legend = TRUE, legend = "bottom",
         labels = c("a", "b")),width = 9,height = 7)

ggsave("./Figures/Trait_effects_resistance_grazing.pdf",
       ggarrange(
         p_graz_CWM, p_graz_FD,
         nrow = 2, ncol = 1,
         common.legend = TRUE, legend = "bottom",
         labels = c("a", "b")),width = 9,height = 7)




## >> SPEI-3, 9  ----

d=readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model.rds")

p_CWM3 = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI03$Partial_residuals,
  d_Rl = d$SPEI03$Partial_residuals,
  grepl_character = "CWM"
)

p_FD3 = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI03$Partial_residuals,
  d_Rl = d$SPEI03$Partial_residuals,
  grepl_character = "FD"
)

ggsave("./Figures/LUI_drought_timescale_resistance_resilience_SPEI3.pdf",
       ggarrange(p_CWM3+ylab("Effect size of mean community trait (spei 3 months)"),
                 p_FD3+ylab("Effect size of community trait diversity (spei 3 months)"),
                 nrow = 2,labels = letters[1:2],common.legend = T,legend = "bottom"),
       width = 8,height = 8)


p_CWM9 = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI09$Partial_residuals,
  d_Rl = d$SPEI09$Partial_residuals,
  grepl_character = "CWM"
)

p_FD9 = Plot_MV_paired_by_axis_join(
  d_Rs = d$SPEI09$Partial_residuals,
  d_Rl = d$SPEI09$Partial_residuals,
  grepl_character = "FD"
)


ggsave("./Figures/LUI_drought_timescale_resistance_resilience_SPEI9.pdf",
       ggarrange(p_CWM9+ylab("Effect size of mean community trait (spei 9 months)"),
                 p_FD9+ylab("Effect size of community trait diversity (spei 9 months)"),
                 nrow = 2,labels = letters[1:2],common.legend = T,legend = "bottom"),
       width = 8,height = 8)

## >> Individual traits ----

d=readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model_all_individual_traits.rds")
pdf("./Test.pdf",width = 9,height = 4)

for (k in names(d$spei06)){
  p1=Plot_MV_figure(Rename_PCA_axes(d$spei06[[k]]$Partial_residuals),grepl_character = "CWM",alpha_ = .6)
  p2=Plot_MV_figure(Rename_PCA_axes(d$spei06[[k]]$Partial_residuals),grepl_character = "FD",alpha_ = .6)
  print(ggarrange(p1,p2,ncol=2,labels = letters[1:2]))
}
dev.off()

## >> Correlation resistance resilience ----

Moving_window_Rs_Rl_correlation = function(d,
                                           X_var,
                                           resistance_col = "Resistance_Isbell",
                                           resilience_col = "Resilience_Isbell_abs",
                                           n_min = 110,
                                           n_step = 2,
                                           method = "spearman") {
  
  d = d %>%
    dplyr::filter(!is.na(.data[[X_var]]),
                  !is.na(.data[[resistance_col]]),
                  !is.na(.data[[resilience_col]])) %>%
    dplyr::arrange(.data[[X_var]]) %>%
    dplyr::mutate(row_id = dplyr::row_number())
  
  results = list()
  index = 0
  
  repeat {
    index = index + 1
    
    id_max = min(n_min + n_step * (index - 1), nrow(d))
    id_min = max(id_max - n_min, 1)
    win = d[id_min:id_max, ]
    
    r_obs = cor(win[[resistance_col]], win[[resilience_col]],
                use = "complete.obs", method = method)
    
    results[[index]] = tibble(
      X_mean   = mean(win[[X_var]], na.rm = TRUE),
      n_events = nrow(win),
      r_obs    = r_obs
    )
    
    if (id_max == nrow(d)) break
  }
  
  bind_rows(results) %>% dplyr::mutate(X_var = X_var)
}
cor_LUI  = Moving_window_Rs_Rl_correlation(d_events, "LUI_global")
cor_SPEI = Moving_window_Rs_Rl_correlation(d_events, "SPEI_value")
ggplot(cor_LUI, aes(x = X_mean, y = r_obs)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Land-use intensity", y = "Correlation (Rs vs. Rl)") +
  theme_bw()
ggplot(cor_SPEI, aes(x = X_mean, y = r_obs)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "Drought intensity", y = "Correlation (Rs vs. Rl)") +
  theme_bw()

## >> R2 models ----

d_lui  = readRDS("./Results/Moving_window_gradient_all_results_Geo_clim_model.rds")
d_spei = readRDS("./Results/Moving_window_gradient_all_results_SPEI.rds")

r2_all = dplyr::bind_rows(
  d_lui$SPEI06$R2_model  %>% dplyr::filter(Stability_var == "Resistance_Isbell") %>%
    dplyr::transmute(Gradient = Variable_gradient, R2m = R2, R2c = R2c,
                     Stability_var = "Resistance", Gradient_type = "LUI"),
  d_lui$SPEI06$R2_model  %>% dplyr::filter(Stability_var == "Resilience_Isbell_abs") %>%
    dplyr::transmute(Gradient = Variable_gradient, R2m = R2, R2c = R2c,
                     Stability_var = "Resilience", Gradient_type = "LUI"),
  d_spei$SPEI06$R2_model %>% dplyr::filter(Stability_var == "Resistance_Isbell") %>%
    dplyr::transmute(Gradient = Variable_gradient, R2m = R2, R2c = R2c,
                     Stability_var = "Resistance", Gradient_type = "SPEI"),
  d_spei$SPEI06$R2_model %>% dplyr::filter(Stability_var == "Resilience_Isbell_abs") %>%
    dplyr::transmute(Gradient = Variable_gradient, R2m = R2, R2c = R2c,
                     Stability_var = "Resilience", Gradient_type = "SPEI")
)

r2_lui_plot = r2_all %>% dplyr::filter(Gradient_type == "LUI") %>%
  tidyr::pivot_longer(c(R2m, R2c), names_to = "R2_type", values_to = "R2_value") %>%
  ggplot(aes(Gradient, R2_value, colour = R2_type, shape = Stability_var)) +
  geom_line(lwd = 1) + geom_point(size = 3, fill = "white") +
  scale_colour_manual(values = c(R2m = "black", R2c = "brown"),
                      labels = c(R2m = expression(R[m]^2), R2c = expression(R[c]^2))) +
  scale_shape_manual(values = c(Resistance = 21, Resilience = 23)) +
  labs(x = "Land-use intensity", y = expression(R^2), colour = NULL, linetype = NULL) +
  the_theme2 + facet_wrap(. ~ Stability_var, ncol = 2) +
  theme(legend.position = "bottom", legend.box = "horizontal")+guides(shape="none")

r2_spei_plot = r2_all %>% dplyr::filter(Gradient_type == "SPEI") %>%
  tidyr::pivot_longer(c(R2m, R2c), names_to = "R2_type", values_to = "R2_value") %>%
  ggplot(aes(-Gradient, R2_value, colour = R2_type, shape = Stability_var)) +
  geom_line(lwd = 1) + geom_point(size = 3, fill = "white") +
  scale_colour_manual(values = c(R2m = "black", R2c = "brown"),
                      labels = c(R2m = expression(R[m]^2), R2c = expression(R[c]^2))) +
  scale_shape_manual(values = c(Resistance = 21, Resilience = 23)) +
  labs(x = "Drought intensity \n <- Moderate drought      Extreme drought ->", y = expression(R^2),
       colour = NULL, linetype = NULL) +
  the_theme2 + facet_wrap(. ~ Stability_var, ncol = 2) +
  theme(legend.position = "bottom", legend.box = "horizontal")+guides(shape="none")

R2_plot = ggarrange(r2_lui_plot, r2_spei_plot, nrow =1, common.legend = TRUE, legend = "bottom")

ggsave("./Figures/R2_models.pdf",R2_plot,width = 11,height = 4)

