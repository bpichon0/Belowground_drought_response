rm(list=ls())
source("0_Functions.R")

d=readRDS("./Data/Merged_datasets_Hairphae.rds")

## 1) Moving windows LUI gradient ----

SPEI_timescale = c("spei03", "spei06", "spei09")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T, F)
datasets = c("HairPhae")
diversity_types = c("FD")
climate_params = c("SPEI_value", "longitude", "latitude","Temperature",
                   "Precipitation","Soil_moisture")
name_functions = c("FD") 


all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 2,
            variable_gradient = "LUI_global",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}


saveRDS(list(SPEI03 =all_results$spei03$HairPhae_FD_phylo,
             SPEI06 =all_results$spei06$HairPhae_FD_phylo,
             SPEI09 = all_results$spei09$HairPhae_FD_phylo,
             SPEI03_no_phylo =all_results$spei03$HairPhae_FD_no_phylo,
             SPEI06_no_phylo =all_results$spei06$HairPhae_FD_no_phylo,
             SPEI09_no_phylo = all_results$spei09$HairPhae_FD_no_phylo),
        "./Results/Moving_window_gradient_all_results_Geo_clim_model.rds")





SPEI_timescale = c("spei06")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T)
N_boot=100
datasets = c("HairPhae")
diversity_types = c("FD")
climate_params = c("SPEI_value", "longitude", "latitude","Temperature",
                   "Precipitation","Soil_moisture")
name_functions = c("FD") 


all_results = tibble()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          for (n_boot in 1:N_boot){
            stability_var = stability_vars[f]
            print(n_boot)
            # Calculate n_min_sites based on number of predictors
            result = Moving_window_gradient_diversity(
              dataset = dataset,
              stability_var = stability_var,
              n_min_sites = NULL,  
              n_site_increase = 6,
              variable_gradient = "LUI_global",
              N_bootstrap = 100,
              shuffle_gradient=T,
              SPEI = SPEI_timescale[timescale_k],
              climate_params = climate_params,
              random = "both",
              correct_clim_residuals = T,
              phylo = phylo,
              drought_category = c("Moderate drought", "Extreme drought"),
              diversity_type = diversity_type
            )
            all_results=rbind(all_results,result$Coeff_variance%>%
                                dplyr::mutate(., N_boot=n_boot))
            
          }
        }
      }
    }
  }
}


saveRDS(all_results,"./Results/Moving_window_gradient_all_results_Geo_clim_model_shuffle.rds")

## 2) Moving windows SPEI gradient ----

SPEI_timescale = c("spei03", "spei06", "spei09")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T,F)
datasets = c( "HairPhae")
diversity_types = c("FD")
climate_params = c("longitude", "latitude","LUI_global",
                   "Temperature","Precipitation","Soil_moisture")
name_functions = c("FD") 

all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 2,
            variable_gradient = "SPEI_value",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}

saveRDS(list(SPEI03 =all_results$spei03$HairPhae_FD_phylo,
             SPEI06 =all_results$spei06$HairPhae_FD_phylo,
             SPEI09 = all_results$spei09$HairPhae_FD_phylo,
             SPEI03_no_phylo =all_results$spei03$HairPhae_FD_no_phylo,
             SPEI06_no_phylo =all_results$spei06$HairPhae_FD_no_phylo,
             SPEI09_no_phylo = all_results$spei09$HairPhae_FD_no_phylo),
        "./Results/Moving_window_gradient_all_results_SPEI.rds")




SPEI_timescale = c("spei06")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T)
N_boot=50
datasets = c("HairPhae")
diversity_types = c("FD")
climate_params = c("LUI_global", "longitude", "latitude","Temperature",
                   "Precipitation","Soil_moisture")
name_functions = c("FD") 


all_results = tibble()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          for (n_boot in 1:N_boot){
            stability_var = stability_vars[f]
            print(n_boot)
            # Calculate n_min_sites based on number of predictors
            result = Moving_window_gradient_diversity(
              dataset = dataset,
              stability_var = stability_var,
              n_min_sites = NULL,  
              n_site_increase = 6,
              variable_gradient = "SPEI_value",
              N_bootstrap = 100,
              shuffle_gradient=T,
              SPEI = SPEI_timescale[timescale_k],
              climate_params = climate_params,
              random = "both",
              correct_clim_residuals = T,
              phylo = phylo,
              drought_category = c("Moderate drought", "Extreme drought"),
              diversity_type = diversity_type
            )
            all_results=rbind(all_results,result$Coeff_variance%>%
                                dplyr::mutate(., N_boot=n_boot))
            
          }
        }
      }
    }
  }
}


saveRDS(all_results,"./Results/Moving_window_gradient_all_results_Geo_clim_model_shuffle_SPEI.rds")

## 3) Moving windows fertilization, mowing, grazing gradients ----

SPEI_timescale = c("spei03", "spei06", "spei09")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T,F)
datasets = c( "HairPhae")
diversity_types = c("FD")
climate_params = c("longitude", "latitude","SPEI_value",
                   "Temperature","Precipitation","Soil_moisture")
name_functions = c("FD") 

all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity_binary(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 2,
            variable_gradient = "TotalFertilization",
            variable_gradient_type ="binned",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}

saveRDS(list(SPEI03 =all_results$spei03$HairPhae_FD_phylo,
             SPEI06 =all_results$spei06$HairPhae_FD_phylo,
             SPEI09 = all_results$spei09$HairPhae_FD_phylo,
             SPEI03_no_phylo =all_results$spei03$HairPhae_FD_no_phylo,
             SPEI06_no_phylo =all_results$spei06$HairPhae_FD_no_phylo,
             SPEI09_no_phylo = all_results$spei09$HairPhae_FD_no_phylo),
        "./Results/Moving_window_gradient_all_results_fertilization.rds")



SPEI_timescale = c("spei03", "spei06", "spei09")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T,F)
datasets = c( "HairPhae")
diversity_types = c("FD")
climate_params = c("longitude", "latitude","SPEI_value",
                   "Temperature","Precipitation","Soil_moisture")
name_functions = c("FD") 

all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity_binary(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 2,
            variable_gradient = "TotalGrazing",
            variable_gradient_type ="binned",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}

saveRDS(list(SPEI03 =all_results$spei03$HairPhae_FD_phylo,
             SPEI06 =all_results$spei06$HairPhae_FD_phylo,
             SPEI09 = all_results$spei09$HairPhae_FD_phylo,
             SPEI03_no_phylo =all_results$spei03$HairPhae_FD_no_phylo,
             SPEI06_no_phylo =all_results$spei06$HairPhae_FD_no_phylo,
             SPEI09_no_phylo = all_results$spei09$HairPhae_FD_no_phylo),
        "./Results/Moving_window_gradient_all_results_grazing.rds")



SPEI_timescale = c("spei03", "spei06", "spei09")
stability_vars = c("Resilience_Isbell_abs", "Resistance_Isbell")
phylo_options = c(T,F)
datasets = c( "HairPhae")
diversity_types = c("FD")
climate_params = c("longitude", "latitude","SPEI_value",
                   "Temperature","Precipitation","Soil_moisture")
name_functions = c("FD") 

all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity_binary(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 2,
            variable_gradient = "TotalMowing",
            variable_gradient_type ="binned",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}

saveRDS(list(SPEI03 =all_results$spei03$HairPhae_FD_phylo,
             SPEI06 =all_results$spei06$HairPhae_FD_phylo,
             SPEI09 = all_results$spei09$HairPhae_FD_phylo,
             SPEI03_no_phylo =all_results$spei03$HairPhae_FD_no_phylo,
             SPEI06_no_phylo =all_results$spei06$HairPhae_FD_no_phylo,
             SPEI09_no_phylo = all_results$spei09$HairPhae_FD_no_phylo),
        "./Results/Moving_window_gradient_all_results_mowing.rds")

## 4) Checking axes of traits ----

veg_data = read.table("./Data/Vegetation/vegetation_species.csv", sep=",", header = T) %>%
  dplyr::rename(Site = Useful_EP_PlotID) %>%
  dplyr::mutate(
    Species = case_when(
      Species == "cf_Crepis_biennis" ~ "Crepis_biennis",
      Species == "Epilobium_cf_tetragonum" ~ "Epilobium_tetragonum",
      Species == "Euphorbia_cf_helioscopia" ~ "Euphorbia_helioscopia",
      Species == "Pinus_cf_mugo" ~ "Pinus_mugo",
      Species == "Ploygonum_cf_lapathifolium" ~ "Ploygonum_lapathifolium",
      Species == "Vicia_cf_lathyroides" ~ "Vicia_lathyroides",
      Species == "Persicaria_lapathifolia" ~ "Ploygonum_lapathifolium",
      Species == "Tripleurospermum_perforatum" ~ "Tripleurospermum_inodorum",
      TRUE ~ Species)) %>%
  dplyr::group_by(Year, Site, Species) %>%
  dplyr::summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop")%>%
  dplyr::mutate(Species = gsub("_aggr.", "", Species))

set.seed(1)
global_traits=read.table("./Data/Global_trait_dataset.csv",sep=";")%>%
  Closer_to_normal_traits(.)

colnames(global_traits)=gsub("_HP","",colnames(global_traits))

traits_all    = c("SRL","RTD","AD","RN","RHL","AMF","RHI",
                  "Seed_mass","SLA_all","Height","LDMC","LeafN","LeafP")
traits_below  = c("SRL","RTD","AD","RN","AMF","RHI")
traits_mixed  = c("SRL","RTD","AD","RN","AMF","SLA_all","LDMC","LeafN","LeafP")
traits_above    = c("SLA_all","Height","LDMC","LeafN","LeafP")

data_all   = global_traits[, traits_all]%>%drop_na(.)
data_below = global_traits[, traits_below]%>%drop_na(.)
data_mixed = global_traits[, traits_mixed]%>%drop_na(.)
data_above = global_traits[, traits_above]%>%drop_na(.)


run_pca = function(df, label) {
  pr = paran::paran(df, iterations = 5000, centile = 95,
                    quietly = TRUE, status = FALSE)
  nfac = max(pr$Retained, 1)
  fit = psych::principal(r = df, nfactors = nfac, rotate = "varimax")
  
  scores = as.data.frame(fit$scores)
  colnames(scores) = paste0(label, "_PC", seq_len(ncol(scores)))
  scores$Species = rownames(df)
  
  list(label = label, nfac = nfac, fit = fit, scores = scores,Vaccounted=fit$Vaccounted,loadings=fit$loadings)
}

res_all   = run_pca(data_all,   "All")
res_below = run_pca(data_below, "Below")
res_mixed = run_pca(data_mixed, "Mixed")
res_above = run_pca(data_above, "Above")

scores_merged = res_all$scores %>%
  dplyr::inner_join(res_below$scores, by = "Species") %>%
  dplyr::inner_join(res_above$scores, by = "Species") %>%
  dplyr::inner_join(res_mixed$scores, by = "Species")

pc_cols = setdiff(colnames(scores_merged), "Species")
cor_mat = cor(scores_merged[, pc_cols], use = "pairwise.complete.obs")

groups = sub("_PC.*", "", pc_cols)
cross_only = outer(groups, groups, FUN = "!=")
cor_mat_cross = cor_mat
cor_mat_cross[!cross_only] = NA


melt_df = reshape2::melt(cor_mat_cross, na.rm = TRUE,
                         varnames = c("Axis1", "Axis2"),
                         value.name = "r")%>%
  dplyr::filter(., Axis1 %!in% c("Mixed_PC1","Mixed_PC2"),
                Axis2 %!in% c("Mixed_PC1","Mixed_PC2"))

p_heatmap = ggplot(melt_df, aes(Axis1, Axis2, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1)) +
  the_theme2+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Pearson r")



p11=plotPCA(res_all$scores,
            res_all$loadings,
            res_all$Vaccounted,annotateFactor = 3,
            xIndex = 1,yIndex = 2,
            xLim = c(-3,3),yLim = c(-3,3))

p12=plotPCA(res_all$scores,
            res_all$loadings,
            res_all$Vaccounted,annotateFactor = 3,
            xIndex = 1,yIndex = 3,
            xLim = c(-3,3),yLim = c(-3,3))


p2=plotPCA(res_below$scores,
           res_below$loadings,
           res_below$Vaccounted,annotateFactor = 3,
           xIndex = 1,yIndex = 2,
           xLim = c(-3,3),yLim = c(-3,3))

p3=plotPCA(res_mixed$scores,
           res_mixed$loadings,
           res_mixed$Vaccounted,annotateFactor = 3,
           xIndex = 1,yIndex = 2,
           xLim = c(-3,3),yLim = c(-3,3))

p4=plotPCA(res_above$scores,
           res_above$loadings,
           res_above$Vaccounted,annotateFactor = 3,
           xIndex = 1,yIndex = 2,
           xLim = c(-3,3),yLim = c(-3,3))


p_tot=ggarrange(p_heatmap,ggarrange(p2+ggtitle("Belowground traits"),p4+ggtitle("Aboveground traits"),nrow=2,ncol=1,labels = c("b","c")),ncol=2,labels = c("a",""))
ggsave("./Figures/PCA_traits_above_below_comparison.pdf",p_tot,width = 8,height = 6)



# Trait CWM
rm(list=ls())
source("0_Functions.R")

traits_all   = c("SRL","RTD","AD","RN","RHL","AMF","RHI",
                 "Seed_mass","SLA_all","Height","LDMC","LeafN","LeafP")

traits_below = c("SRL","RTD","AD","RN","AMF")

traits_mixed = c("SRL","RTD","AD","RN","AMF","SLA_all","LDMC","LeafN","LeafP")


traits_above   = c("SLA_all","Height","LDMC","LeafN","LeafP")

global_traits=read.table("./Data/Global_trait_dataset.csv",sep=";")%>%
  Closer_to_normal_traits(.)

colnames(global_traits)=gsub("_HP","",colnames(global_traits))
trait_list = list(
  All   = run_pca_scores(global_traits[, traits_all]%>%drop_na(.)),
  Below = run_pca_scores(global_traits[, traits_below]%>%drop_na(.)),
  Above = run_pca_scores(global_traits[, traits_above]%>%drop_na(.)),
  Mixed = run_pca_scores(global_traits[, traits_mixed]%>%drop_na(.))
)
compute_FD_CWM = function(Cover_k, Species_k, trait_sp) {
  
  n_pc_cols = grep("^PC[0-9]+$", colnames(trait_sp), value = TRUE)
  
  # ---- empty / degenerate community ----
  if (!any(Cover_k) | length(Cover_k) == 0) {
    empty = tibble(
      Species_richness = NA_real_,
      FEve_all = NA, FDis_all = NA, FDiv_all = NA, RaoQ_all = NA,
      RaoQ_all_FD = NA,
      Trait_representativity = NA
    )
    for (k in seq_along(n_pc_cols)) {
      empty[[paste0("FEve_PC", k)]]   = NA
      empty[[paste0("FDis_PC", k)]]   = NA
      empty[[paste0("FDiv_PC", k)]]   = NA
      empty[[paste0("RaoQ_PC", k)]]   = NA
      empty[[paste0("RaoQ_PC", k, "_FD")]] = NA
      empty[[paste0("CWM_PC", k)]]    = NA
      empty[[paste0("CWvar_PC", k)]]  = NA
      empty[[paste0("CWskew_PC", k)]] = NA
      empty[[paste0("CWkurt_PC", k)]] = NA
      empty[[paste0("FD_PC", k)]]     = NA
    }
    return(empty)
  }
  
  if (any(Cover_k == 0)) {
    which_sp = which(Cover_k == 0)
    Cover_k  = Cover_k[-which_sp]
    Species_k = Species_k[-which_sp]
  }
  
  Trait_representativity = 100
  if (any(Species_k %!in% trait_sp$scientificName)) {
    which_sp = which(Species_k %!in% trait_sp$scientificName)
    Trait_representativity = 100 * (1 - sum(Cover_k[which_sp]) / sum(Cover_k))
    Cover_k   = Cover_k[-which_sp]
    Species_k = Species_k[-which_sp]
  }
  
  trait_sp_k = dplyr::filter(trait_sp, scientificName %in% Species_k)
  
  Cover_k = as.matrix(t(Cover_k))
  colnames(Cover_k) = Species_k
  
  Trait_matrix = as.matrix(trait_sp_k[, n_pc_cols, drop = FALSE])
  rownames(Trait_matrix) = trait_sp_k$scientificName
  Trait_matrix = Trait_matrix[colnames(Cover_k), , drop = FALSE]  # align order
  
  Dist_matrix = compute_dist_matrix(Trait_matrix, metric = "euclidean")
  Dist_matrix = Dist_matrix / max(Dist_matrix, na.rm = TRUE)
  FD_trait_weighted = FD::dbFD(Dist_matrix, Cover_k)
  
  RaoQ_traits_all = RaoQ_traits(Cover_k, Dist_matrix)
  RaoQ_traits_all[RaoQ_traits_all == 1] = NA
  
  out = tibble(
    Species_richness = length(Species_k),
    FEve_all = FD_trait_weighted$FEve,
    FDis_all = FD_trait_weighted$FDis,
    FDiv_all = FD_trait_weighted$FDiv,
    RaoQ_all = RaoQ_traits_all,
    RaoQ_all_FD = FD_trait_weighted$RaoQ,
    Trait_representativity = Trait_representativity
  )
  
  for (k in seq_along(n_pc_cols)) {
    pc_name = n_pc_cols[k]
    
    fd_k = tryCatch({
      Dist_k = compute_dist_matrix(Trait_matrix[, pc_name], metric = "euclidean")
      Dist_k = Dist_k / max(Dist_k, na.rm = TRUE)
      FD::dbFD(Dist_k, Cover_k)
    }, error = function(e) list(FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA))
    
    raoq_k = tryCatch({
      Dist_k = compute_dist_matrix(Trait_matrix[, pc_name], metric = "euclidean")
      Dist_k = Dist_k / max(Dist_k, na.rm = TRUE)
      rq = RaoQ_traits(Cover_k, Dist_k)
      rq[rq == 1] = NA
      rq
    }, error = function(e) NA)
    
    cover_norm = as.numeric(Cover_k) / sum(Cover_k)
    val = Trait_matrix[, pc_name]
    cwm  = sum(cover_norm * val, na.rm = TRUE)
    cvar = sum(cover_norm * (val - cwm)^2, na.rm = TRUE)
    cskew = sum(cover_norm * (val - cwm)^3, na.rm = TRUE) / cvar^(3/2)
    ckurt = sum(cover_norm * (val - cwm)^4, na.rm = TRUE) / cvar^2
    fd_disp = sum(cover_norm * (abs(val - cwm) / sum(abs(val - cwm))))
    
    out[[paste0("FEve_PC", k)]]        = fd_k$FEve
    out[[paste0("FDis_PC", k)]]        = fd_k$FDis
    out[[paste0("FDiv_PC", k)]]        = fd_k$FDiv
    out[[paste0("RaoQ_PC", k)]]        = raoq_k
    out[[paste0("RaoQ_PC", k, "_FD")]] = fd_k$RaoQ
    out[[paste0("CWM_PC", k)]]         = cwm
    out[[paste0("CWvar_PC", k)]]       = cvar
    out[[paste0("CWskew_PC", k)]]      = cskew
    out[[paste0("CWkurt_PC", k)]]      = ckurt
    out[[paste0("FD_PC", k)]]          = fd_disp
  }
  
  out
}


compute_FD_CWM_all_subsets = function(Cover_k, Species_k, trait_list) {
  shared_cols = c("Species_richness", "Trait_representativity")
  
  results = purrr::imap(trait_list, function(trait_sp, label) {
    res = compute_FD_CWM(Cover_k, Species_k, trait_sp)
    to_prefix = setdiff(colnames(res), shared_cols)
    colnames(res)[colnames(res) %in% to_prefix] = paste0(label, "_", to_prefix)
    res
  })
  
  # keep shared cols from the first subset only, drop from the rest
  for (i in seq_along(results)) {
    if (i > 1) results[[i]] = dplyr::select(results[[i]], -any_of(shared_cols))
  }
  
  dplyr::bind_cols(results)
}

Trait_communities = read.table("./Data/Vegetation/vegetation_species.csv",
                               sep = ",", header = TRUE) %>%
  dplyr::rename(Site = Useful_EP_PlotID) %>%
  dplyr::mutate(
    Species = case_when(
      Species == "cf_Crepis_biennis" ~ "Crepis_biennis",
      Species == "Epilobium_cf_tetragonum" ~ "Epilobium_tetragonum",
      Species == "Euphorbia_cf_helioscopia" ~ "Euphorbia_helioscopia",
      Species == "Pinus_cf_mugo" ~ "Pinus_mugo",
      Species == "Ploygonum_cf_lapathifolium" ~ "Ploygonum_lapathifolium",
      Species == "Vicia_cf_lathyroides" ~ "Vicia_lathyroides",
      Species == "Persicaria_lapathifolia" ~ "Ploygonum_lapathifolium",
      Species == "Tripleurospermum_perforatum" ~ "Tripleurospermum_inodorum",
      TRUE ~ Species)) %>%
  dplyr::group_by(Year, Site, Species) %>%
  dplyr::summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Species = gsub("_aggr.", "", Species)) %>%
  dplyr::filter(!is.na(Cover)) %>%
  dplyr::group_by(Year, Site) %>%
  dplyr::group_modify(~ compute_FD_CWM_all_subsets(.x$Cover, .x$Species, trait_list))


saveRDS(object = Trait_communities,file = "./Data/Community_diversity_HairPhae_above_below_mixed.rds")


# 
# 
# TEST_traits = readRDS("./TEST_traits.rds")%>%
#   dplyr::filter(., Trait_representativity>70)
# 
# Plot_correlation_variables(TEST_traits[,colnames(TEST_traits)[grep("CWM",colnames(TEST_traits))]])
# 










# Model results

Moving_window_gradient_diversity = function(dataset = "Underplot",
                                            stability_var,
                                            n_min_sites = NULL,
                                            n_site_increase = 10,
                                            variable_gradient = "LUI_global",
                                            N_bootstrap = 100,
                                            SPEI = "spei06",
                                            climate_params = c(""),
                                            correct_clim_residuals=T,
                                            random = "both",
                                            phylo=T,
                                            drought_category = c("Moderate drought", "Extreme drought"),
                                            diversity_type = c("FD")) {
  
  if (diversity_type == "RaoQ") {
    diversity_predictors = c("Species_richness",
                             "RaoQ_PC1", "RaoQ_PC2", "RaoQ_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type=="FD_below") {  
    diversity_predictors = c("Species_richness",
                             "Below_FD_PC1", "Below_FD_PC2",
                             "Below_CWM_PC1", "Below_CWM_PC2")
  } else if (diversity_type=="FD_above") {  
    diversity_predictors = c("Species_richness",
                             "Above_FD_PC1", "Above_FD_PC2",
                             "Above_CWM_PC1", "Above_CWM_PC2")
  }else if (diversity_type=="FD_mixed") {  
    diversity_predictors = c("Species_richness",
                             "Mixed_FD_PC1", "Mixed_FD_PC2",
                             "Mixed_CWM_PC1", "Mixed_CWM_PC2")
  } else if (diversity_type=="FDis") {  
    diversity_predictors = c("Species_richness",
                             "FDis_PC1", "FDis_PC2", "FDis_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  }else if (diversity_type=="all") {  
    diversity_predictors = c("Species_richness",
                             "FDis_PC1", "FDis_PC2", "FDis_PC3",
                             "FEve_PC1", "FEve_PC2", "FEve_PC3",
                             "FDiv_PC1", "FDiv_PC2", "FDiv_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  }
  
  if (phylo){
    diversity_predictors=c(diversity_predictors,"PSV", "MNTD_weighted")
  }
  
  if (is.null(n_min_sites)){
    n_min_sites=(length(diversity_predictors)+
                   ifelse(length(climate_params)==1 | correct_clim_residuals,0,length(climate_params))+3)*10
  }
  
  data_path = paste0("./Data/Merged_datasets_", dataset, ".rds")
  
  if (stability_var == "Temporal_stability") {  
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., Trait_representativity > 70) %>%
      dplyr::mutate(., Species_richness = sqrt(Species_richness)) %>%
      dplyr::distinct(., .keep_all = T, Site) 
    
  } else if (stability_var == "Resilience_Isbell_abs") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Drought_category %in% drought_category,
                    !is.na(Resilience_Isbell_abs), 
                    is.na(Recovery_context)) %>%
      dplyr::mutate(., 
                    Resilience_Isbell_abs = log(Resilience_Isbell_abs),
                    Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Resilience_Isbell") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Resilience_Isbell > 0,
                    Drought_category %in% drought_category,
                    !is.na(Resilience_Isbell), 
                    is.na(Recovery_context)) %>%
      dplyr::mutate(., 
                    Resilience_Isbell = log(Resilience_Isbell),
                    Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Recovery") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Drought_category %in% drought_category,
                    !is.na(Recovery), 
                    is.na(Recovery_context)) %>%
      dplyr::mutate(., 
                    Recovery = bcPower(Recovery,-.2),
                    Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Resistance_Isbell") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Resistance_Isbell > 0,
                    Drought_category %in% drought_category) %>%
      dplyr::mutate(., 
                    Resistance_Isbell = bcPower(Resistance_Isbell, -0.7),
                    Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Logbiomass_t_median") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Drought_category %in% drought_category,
                    !is.na(Logbiomass_t_median)) %>%
      dplyr::mutate(., Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Sensitivity_water_biomass") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    !is.na(Sensitivity_water_biomass)) %>%
      dplyr::mutate(., Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Sensitivity_water_cover") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    !is.na(Sensitivity_water_cover)) %>%
      dplyr::mutate(., Species_richness = sqrt(Species_richness))
    
  } else if (stability_var == "Logbiomass_t_minus_1") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Drought_category_previous %in% drought_category,
                    !is.na(Logbiomass_t_minus_1)) %>%
      dplyr::mutate(., Species_richness = sqrt(Species_richness))
  }else if (stability_var == "Resilience_wet") {
    
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., 
                    Trait_representativity > 70,
                    Resistance_Isbell < 0,
                    Drought_category %in% drought_category) %>%
      dplyr::mutate(., 
                    Resistance_Isbell = bcPower(Resistance_Isbell, -0.7),
                    Species_richness = sqrt(Species_richness))
    
  }
  
  d2=readRDS("./Data/Community_diversity_HairPhae_above_below_mixed.rds")
  
  d_file=merge(d_file,d2%>%dplyr::select(.,-Species_richness,-Trait_representativity),by=c("Site","Year"))
  
  # select data
  d_mod = d_file %>%
    dplyr::mutate(., 
                  Variable_gradient = as.numeric(dplyr::pull(., dplyr::all_of(variable_gradient))),
                  Response_var = as.numeric(dplyr::pull(., dplyr::all_of(stability_var)))) %>%
    dplyr::select(., 
                  Response_var,
                  Variable_gradient,
                  Site, Region, any_of(climate_params),
                  Year,
                  all_of(diversity_predictors)) %>%
    dplyr::arrange(., Variable_gradient) %>%
    dplyr::mutate(., Variable_gradient = as.character(Variable_gradient)) %>%
    drop_na(.)
  
  # correcting with residuals
  if (climate_params[1] != "" & correct_clim_residuals) {
    d_mod$Response_var = residuals(lm(as.formula(paste0("Response_var ~ ", paste0(climate_params, collapse = "+"))), 
                                      data = d_mod))
    d_mod = d_mod %>%
      dplyr::select(., -any_of(climate_params))
  }
  
  
  # tibbles
  all_results = all_variance = all_R2 = 
    all_partial_residuals = all_spatial_autocorrelation =
    all_temporal_autocorrelation = tibble()
  
  print(dim(d_mod))
  
  # Main loop
  for (index in 1:326) {
    
    if (index %% 20 == 0) print(index)
    
    id_max = ifelse((n_min_sites + n_site_increase * (index - 1)) > nrow(d_mod),
                    nrow(d_mod), 
                    (n_min_sites + n_site_increase * (index - 1)))
    id_min = ifelse((id_max - n_min_sites) < 1, 1, id_max - n_min_sites)
    save_d = d_mod[id_min:id_max, ]
    
    # Scale data
    data_model = d_mod[id_min:id_max, ] %>%
      dplyr::mutate(., Variable_gradient=as.numeric(Variable_gradient))%>%
      Numeric_scaling_df(.)
    
    # Fit LMM with random effects
    if (random == "both") {
      mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Site) + (1|Year) + ",
                                           paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                         data = data_model, REML = T, control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
    } else if (random %in% c("Site", "site")) {
      mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Site) + ",
                                           paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                         data = data_model, REML = T, control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
    } else if (random %in% c("Year", "year")) {
      mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Year) + ",
                                           paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                         data = data_model, REML = T, control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
    } else if (random %in% c("Region", "region")) {
      mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Region) + ",
                                           paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                         data = data_model, REML = T, control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
    } else if (random=="lme"){
      # mod_drought=lme(as.formula(paste0("Response_var ~ ",paste0(colnames(data_model %>%
      #                                                                       dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
      #                 random =~ 1|Site/Year,
      #          data = data_model, method = "ML", na.action = na.fail,
      #          control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
      mod_drought=gls(as.formula(paste0("Response_var ~ ",paste0(colnames(data_model %>%
                                                                            dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                      data = data_model, method = "ML", na.action = na.fail,
                      control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
      mod_drought=update(mod_drought,correlation = corCAR1(form = as.formula(paste0("~ Year|Site"))))
      
    } else if (random == "GAM"){
      
      predictors = colnames(data_model %>% 
                              dplyr::select(-Response_var, -Site, -Region))
      
      spatial_predictors = c("longitude", "latitude") 
      linear_predictors = predictors[!predictors %in% spatial_predictors]
      
      
      mod_drought = gam(
        as.formula(paste0(
          "Response_var ~ ",
          paste0(linear_predictors, collapse = " + "),  
          " + s(longitude, latitude)"
        )),
        data = data_model,
        method = "REML",
        na.action = na.omit)      
      
      mod_drought = gamm(
        as.formula(paste0(
          "Response_var ~ ",
          paste0(linear_predictors, collapse = " + "),  # Linear terms
          " + s(longitude, latitude) + s(Site, bs='re')"
        )),
        data = data_model,
        method = "REML",
        na.action = na.omit)
      
    } 
    
    if (random !="lme"){
      # Extract R2
      all_R2 = rbind(all_R2, 
                     tibble(R2 = r.squaredGLMM(mod_drought)[1], 
                            R2c = r.squaredGLMM(mod_drought)[2]) %>%
                       add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      # Bootstrap coefficients
      if (random!="GAM"){
        boot_strapmod = fixef(arm::sim(mod_drought))[, -1]
        boot_mod = list()
        boot_mod$t = boot_strapmod
      }else{
        boot_mod = list()
        boot_mod$t = tibble()
      }
      
      # Partial residuals for each predictor
      df = tibble()
      for (pred in diversity_predictors) {
        
        partial_res = visreg::visreg(mod_drought, pred, plot = F)
        mod_cov = boot(data = partial_res$res, 
                       statistic = boot_function_lm_partialres, 
                       R = 99, 
                       formula = as.formula(paste0("visregRes ~", pred)))
        
        mod_cov$t = as.data.frame(mod_cov$t)
        colnames(mod_cov$t) = pred
        if (nrow(df) == 0) {
          df = mod_cov$t
        } else {
          df = bind_cols(df, mod_cov$t)
        }
      }
      
      all_partial_residuals = rbind(all_partial_residuals, 
                                    df %>% add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      # Store results
      all_variance = rbind(all_variance, 
                           data.frame(Effect_size = colMeans(boot_mod$t),
                                      Predictor = colnames(boot_mod$t)) %>%
                             add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      all_results = rbind(all_results, 
                          data.frame(boot_mod$t) %>%
                            add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      # Test for spatial autocorrelation
      #because we have multiple obs per Plot (years of events), we need to recalculate the residuals per plot
      data_model$location_id = paste(d_file$longitude[match(data_model$Site, d_file$Site)], 
                                     d_file$latitude[match(data_model$Site, d_file$Site)],
                                     sep = "_")
      
      # we recompute residuals
      sim_res_agg = recalculateResiduals(simulateResiduals(mod_drought), 
                                         group = data_model$location_id)
      
      unique_coords = data_model %>%
        dplyr::mutate(.,
                      longitude=d_file$longitude[match(data_model$Site, d_file$Site)],
                      latitude=d_file$latitude[match(data_model$Site, d_file$Site)])%>%
        dplyr::distinct(longitude, latitude) %>%
        dplyr::arrange(longitude, latitude)
      
      spatial_test = testSpatialAutocorrelation(sim_res_agg,
                                                x = unique_coords$longitude,
                                                y = unique_coords$latitude,plot = F)
      
      all_spatial_autocorrelation=rbind(all_spatial_autocorrelation,
                                        tibble(Statistic=spatial_test$statistic[1],Pvalue=spatial_test$p.value)%>%
                                          add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      
      #Same with temporal autocorrelation
      sim_res_agg_time = recalculateResiduals(simulateResiduals(mod_drought), 
                                              group = data_model$Year)
      
      unique_years = data_model %>%
        dplyr::distinct(Year) %>%
        dplyr::arrange(Year)
      
      temporal_test = testTemporalAutocorrelation(sim_res_agg_time,
                                                  time = unique_years$Year,plot = F)
      
      all_temporal_autocorrelation = rbind(all_temporal_autocorrelation,
                                           tibble(Statistic = temporal_test$statistic[1],
                                                  Pvalue = temporal_test$p.value) %>%
                                             add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      
    }
    
    if (id_max == nrow(d_mod)) break
  }
  
  return(list(
    Effect_size = all_results %>%
      add_column(., Stability_var = stability_var),
    Coeff_variance = all_variance %>%
      add_column(., Stability_var = stability_var),
    R2_model = all_R2 %>%
      add_column(., Stability_var = stability_var),
    Partial_residuals = all_partial_residuals %>%
      add_column(., Stability_var = stability_var),
    Spatial_autorrelation = all_spatial_autocorrelation %>%
      add_column(., Stability_var = stability_var),
    Temporal_autocorrelation = all_temporal_autocorrelation %>%
      add_column(., Stability_var = stability_var)
  ))
}


SPEI_timescale = c( "spei06")
stability_vars = c("Resistance_Isbell")
phylo_options = c(T)
datasets = c("HairPhae")
diversity_types = c("FD_mixed","FD_below","FD_above")
climate_params = c("SPEI_value", "longitude", "latitude","Temperature",
                   "Precipitation","Soil_moisture")
name_functions = c("FD") 


all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 2,
            variable_gradient = "LUI_global",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}

saveRDS(all_results,"./Results/Robustness_trait_axes_results.rds")

SPEI_timescale = c( "spei06")
stability_vars = c("Resistance_Isbell")
phylo_options = c(T)
datasets = c("HairPhae")
diversity_types = c("FD_mixed","FD_below","FD_above")
climate_params = c("LUI_global", "longitude", "latitude","Temperature",
                   "Precipitation","Soil_moisture")
name_functions = c("FD") 


all_results = list()
for (timescale_k in 1:length(SPEI_timescale)) {
  for (dataset in datasets) {
    for (diversity_type in diversity_types) {
      for (phylo in phylo_options) {
        
        
        results_list = list()
        
        for (f in 1:length(stability_vars)) {
          
          stability_var = stability_vars[f]
          
          # Calculate n_min_sites based on number of predictors
          result = Moving_window_gradient_diversity(
            dataset = dataset,
            stability_var = stability_var,
            n_min_sites = NULL,  
            n_site_increase = 4,
            variable_gradient = "SPEI_value",
            N_bootstrap = 100,
            SPEI = SPEI_timescale[timescale_k],
            climate_params = climate_params,
            random = "both",
            correct_clim_residuals = T,
            phylo = phylo,
            drought_category = c("Moderate drought", "Extreme drought"),
            diversity_type = diversity_type
          )
          results_list[[f]]=result
        }
        
        combo_id = paste(dataset, diversity_type, ifelse(phylo, "phylo", "no_phylo"), sep = "_")
        
        # Store in nested structure
        all_results[[SPEI_timescale[timescale_k]]][[combo_id]] = list(
          Partial_residuals = bind_rows(lapply(results_list, function(x) x$Partial_residuals)),
          Coeff_variance = bind_rows(lapply(results_list, function(x) x$Coeff_variance)),
          R2_model = bind_rows(lapply(results_list, function(x) x$R2_model)),
          Effect_size = bind_rows(lapply(results_list, function(x) x$Effect_size)),
          Temporal_autocorrelation = bind_rows(lapply(results_list, function(x) x$Temporal_autocorrelation)),
          Spatial_autocorrelation = bind_rows(lapply(results_list, function(x) x$Spatial_autorrelation))
        )
      }
    }
  }
}

saveRDS(all_results,"./Results/Robustness_trait_axes_results_SPEI.rds")


all_results=readRDS("./Results/Robustness_trait_axes_results.rds")
Plot_MV_figure=function(d,stability_var="Resistance_Isbell",CI=90,
                        grepl_character="CWM",
                        add_signif=T,phylo=F,alpha_=.08){
  
  d_plot = (d) %>%
    dplyr::filter(., Stability_var == stability_var) %>%
    dplyr::select(-Stability_var) %>%
    reshape2::melt(., id.vars = "Variable_gradient") %>%
    dplyr::rename(var = variable, value = value) %>%
    dplyr::filter(grepl(grepl_character, var)) %>%
    dplyr::filter(!is.na(value))
  
  p=ggplot(d_plot, aes(x = Variable_gradient, y = value, group = var)) +
    geom_smooth(aes(color = var),se = F) +
    stat_summary(aes(fill = var), 
                 fun.data = function(x) {
                   data.frame(ymin = quantile(x, (1-(CI/100))/2),
                              ymax = quantile(x, ((CI/100))+(1-(CI/100))/2))},
                 geom = "ribbon", alpha = alpha_) +
    facet_wrap(~var, scales = "free_y", nrow = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 1) +
    labs(x = "Land-use intensity", 
         y = "Effect size on drought resistance",
         title = "") +
    the_theme2 +
    theme(strip.background = element_rect(fill = "black"),
          strip.text = element_text(colour = "white"),
          legend.position = "none") +
    scale_color_manual(values = color_trait_axes) +
    scale_fill_manual(values = fill_trait_axes)
  
  if (phylo){
    p=p+  scale_color_manual(values=rep("black",3))+
      scale_fill_manual(values=rep("grey",3))
  }
  
  if (add_signif){
    
    
    # compute 90% CI per var x Variable_gradient, flag significance
    d_sig = d_plot %>%
      dplyr::group_by(var, Variable_gradient) %>%
      dplyr::summarise(lwr = quantile(value, (1-(CI/100))/2),
                       upr = quantile(value, ((CI/100))+(1-(CI/100))/2),
                       .groups = "drop") %>%
      dplyr::mutate(significant = (lwr > 0) | (upr < 0))
    
    # y-position for the significance marks: just above each facet's max value
    d_sig = d_sig %>%
      dplyr::group_by(var) %>%
      dplyr::mutate(y_mark = max(upr, na.rm = TRUE) * 1.08) %>%
      dplyr::ungroup()
    
    p=p+ geom_point(data = dplyr::filter(d_sig, significant),
                    aes(x = Variable_gradient, y = y_mark, group = var),
                    shape = 108, size = 3, color = "grey")
    
  }
  
  return(p)
}


p1=Plot_MV_figure((all_results$spei06$HairPhae_FD_below_phylo$Partial_residuals),grepl_character = "CWM",alpha_ = .6)
p2=Plot_MV_figure((all_results$spei06$HairPhae_FD_below_phylo$Partial_residuals),grepl_character = "FD",alpha_ = .6)

ggsave("./Figures/Trait_effects_resistance_LUI_BELOWGROUND.pdf",
       ggarrange(p1,p2,nrow=2,labels = letters[1:2]),width = 6,height = 7)


p1=Plot_MV_figure((all_results$spei06$HairPhae_FD_above_phylo$Partial_residuals),grepl_character = "CWM",alpha_ = .6)
p2=Plot_MV_figure((all_results$spei06$HairPhae_FD_above_phylo$Partial_residuals),grepl_character = "FD",alpha_ = .6)
ggsave("./Figures/Trait_effects_resistance_LUI_ABOVEGROUND.pdf",
       ggarrange(p1,p2,nrow=2,labels = letters[1:2]),width = 6,height = 7)

p1=Plot_variance_partitioning((all_results$spei06$HairPhae_FD_below_phylo$Coeff_variance))

## 5) Paired covariation----

paired_Rs = Moving_window_gradient_diversity_paired(
  dataset                = "HairPhae",
  stability_var          = "Resistance_Isbell",
  n_site_increase        = 2,
  variable_gradient      = "LUI_global",
  SPEI                   = "spei06",
  climate_params         = c("longitude", "latitude", "SPEI_value",
                             "Temperature", "Precipitation", "Soil_moisture"),
  correct_clim_residuals = TRUE,
  random                 = "both",
  phylo                  = TRUE,
  diversity_type         = "FD"
)

paired_Rl = Moving_window_gradient_diversity_paired(
  dataset                = "HairPhae",
  stability_var          = "Resilience_Isbell_abs",
  n_site_increase        = 2,
  variable_gradient      = "LUI_global",
  SPEI                   = "spei06",
  climate_params         = c("longitude", "latitude", "SPEI_value",
                             "Temperature", "Precipitation", "Soil_moisture"),
  correct_clim_residuals = TRUE,
  random                 = "both",
  phylo                  = TRUE,
  diversity_type         = "FD"
)

paired_Rs2 = Moving_window_gradient_diversity_paired(
  dataset                = "HairPhae",
  stability_var          = "Resistance_Isbell",
  n_site_increase        = 2,
  variable_gradient      = "SPEI_value",
  SPEI                   = "spei06",
  climate_params         = c("longitude", "latitude", "LUI_global",
                             "Temperature", "Precipitation", "Soil_moisture"),
  correct_clim_residuals = TRUE,
  random                 = "both",
  phylo                  = TRUE,
  diversity_type         = "FD"
)

paired_Rl2 = Moving_window_gradient_diversity_paired(
  dataset                = "HairPhae",
  stability_var          = "Resilience_Isbell_abs",
  n_site_increase        = 2,
  variable_gradient      = "SPEI_value",
  SPEI                   = "spei06",
  climate_params         = c("longitude", "latitude", "LUI_global",
                             "Temperature", "Precipitation", "Soil_moisture"),
  correct_clim_residuals = TRUE,
  random                 = "both",
  phylo                  = TRUE,
  diversity_type         = "FD"
)

saveRDS(list(Res_LUI=paired_Rs,
             Rel_LUI=paired_Rl,
             Res_SPEI=paired_Rs2,
             Rel_SPEI=paired_Rl2),
        "./Results/Paired_results.rds")



#null approach to control for correlation 

n_perm = 50

d_results=tibble()
for (k in seq_len(n_perm)) {
  
  message("Permutation ", k, " / ", n_perm)
  
  paired_Rs = Moving_window_gradient_diversity_paired(
    dataset                = "HairPhae",
    stability_var          = "Resistance_Isbell",
    n_site_increase        = 10,
    variable_gradient      = "LUI_global",
    SPEI                   = "spei06",
    climate_params         = c("longitude", "latitude", "SPEI_value",
                               "Temperature", "Precipitation", "Soil_moisture"),
    correct_clim_residuals = TRUE,
    random                 = "both",
    shuffle_traits         = TRUE,
    phylo                  = TRUE,
    diversity_type         = "FD"
  )
  
  paired_Rl = Moving_window_gradient_diversity_paired(
    dataset                = "HairPhae",
    stability_var          = "Resilience_Isbell_abs",
    n_site_increase        = 10,
    variable_gradient      = "LUI_global",
    SPEI                   = "spei06",
    climate_params         = c("longitude", "latitude", "SPEI_value",
                               "Temperature", "Precipitation", "Soil_moisture"),
    correct_clim_residuals = TRUE,
    random                 = "both",
    shuffle_traits         = TRUE,
    phylo                  = TRUE,
    diversity_type         = "FD"
  )
  
  d_slope_Rs = paired_Rs$Partial_residuals %>%
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
  
  d_slope_Rl = paired_Rl$Partial_residuals %>%
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
  
  d_slope = dplyr::inner_join(d_slope_Rs, d_slope_Rl,
                              by = c("Predictor", "Variable_gradient")) %>%
    dplyr::filter(grepl("CWM|FD|Species|PSV|MNTD", Predictor)) %>%
    dplyr::group_by(Predictor) %>%
    dplyr::summarise(
      slope     = coef(lm(Effect_size_Rl ~ Effect_size_Rs))[2],
      intercept = coef(lm(Effect_size_Rl ~ Effect_size_Rs))[1],
      r         = cor(Effect_size_Rs, Effect_size_Rl,
                      method = "spearman", use = "complete.obs"),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(Perm = k)
  
  d_results=rbind(d_results,d_slope)
}

saveRDS(d_results, "./Results/Null_paired_LUI.rds")




