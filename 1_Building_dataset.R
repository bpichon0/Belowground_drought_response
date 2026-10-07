rm(list=ls())
source("0_Functions.R")
# ---------------------- Step 1: Getting SPEI & resistance/resilience  ----

## >> 1) Thresholds for droughts (Moderate & extreme) ----
Sites_info=read.table("./Data/Geography_names_sites.csv",sep=";")

Get_thresholds=readRDS("./Data/SPEI_data.rds")%>%
  lapply(.,function(x){
    
    colnames(x)[1:150]=Sites_info$Name #plot names
    x$Time=as.Date(x$Time)
    x$Month=as.numeric(format(x$Time, "%m"))
    
    moderate_drought_threshold=apply(x[,1:150],2,quantile,.25,na.rm=T) #moderate drought, quantile 25%
    extreme_drought_threshold=apply(x[,1:150],2,quantile,.1,na.rm=T) #extreme drought, quantile 10%
    moderate_wet_threshold=apply(x[,1:150],2,quantile,.75,na.rm=T) #same for wet
    extreme_wet_threshold=apply(x[,1:150],2,quantile,.9,na.rm=T) #same for wet
    
    return(data.frame(Moderate_threshold_drought=moderate_drought_threshold,
                      Extreme_threshold_drought=extreme_drought_threshold,
                      Extreme_threshold_wet=extreme_wet_threshold,
                      Moderate_threshold_wet=moderate_wet_threshold,
                      Site=Sites_info$Name)
    )})


Info_plots=list(Info_site=Sites_info,
                Threshold_droughts=Get_thresholds)

saveRDS(Info_plots,"./Data/Info_plots_thresholds.rds")

## >> 2) Merging SPEI and vegetation data ----
Sites_info=read.table("./Data/Geography_names_sites.csv",sep=";")


Info_plots=readRDS("./Data/Info_plots_thresholds.rds")
d_spei=readRDS("./Data/SPEI_data.rds")

veg_data=read.table("./Data/Vegetation/vegetation_species.csv",sep=",",header = T)%>%
  dplyr::group_by(., Useful_EP_PlotID,Year)%>%
  dplyr::summarise(., .groups = "keep",Cover=sum(Cover,na.rm = T))%>%
  merge(., readxl::read_xlsx("./Data/Vegetation/vegetation_community.xlsx")%>%
          dplyr::mutate(., Year=as.numeric(format(Year, "%Y")))%>%
          dplyr::select(., Year,day_of_year,biomass,Useful_EP_PlotID)%>%
          dplyr::rename(., Total_biomass=biomass,Sampling_day=day_of_year),by=c("Useful_EP_PlotID","Year"),all.x=T)%>%
  dplyr::rename(., Site=Useful_EP_PlotID)%>%
  dplyr::mutate(., Month=sapply(1:nrow(.),function(x){
    return(get_month_from_day(as.numeric(.$Sampling_day[x]),as.numeric(.$Year[x])))}))

#Adding drought categories into SPEI list
spei_cat=lapply(1:length(d_spei), FUN = function(x) {
  
  df_x=d_spei[[x]]
  
  colnames(df_x)[1:150]=Sites_info$Name #plot names
  df_x$Time=as.Date(df_x$Time)
  df_x$Month=as.numeric(format(df_x$Time, "%m"))
  
  df_x=df_x%>%
    dplyr::filter(., Year >= range(veg_data$Year)[1] & Year <= range(veg_data$Year)[2])%>%
    dplyr::filter(., Month %in% c(4,5,6,7))  #april, may, june, july
  
  d_melt=reshape2::melt(df_x,id.vars = c("Month","Year","Time"),variable.name = "Site")
  
  d_melt = merge(d_melt, Info_plots$Threshold_droughts[[x]], by = "Site", all.x = TRUE)
  
  d_melt$drought_category = with(d_melt, 
                                 ifelse(value <= Extreme_threshold_drought, "Extreme drought",
                                        ifelse(value <= Moderate_threshold_drought, "Moderate drought",
                                               ifelse(value <= Moderate_threshold_wet, "Non-drought",
                                                      ifelse(value <= Extreme_threshold_wet, "Moderate wet","Extreme wet")))))
  return(d_melt)
})
names(spei_cat)=names(d_spei)

standardize_site_names = function(x) {
  prefix = gsub("[0-9]", "", x)  
  number = as.numeric(gsub("[A-Z]", "", x))
  paste0(prefix, sprintf("%02d", number))
}

Merged_data=lapply(1:length(spei_cat), function(spei_timescale_k){
  
  spei_cat_k=spei_cat[[spei_timescale_k]]
  spei_cat_k$Site = standardize_site_names(spei_cat_k$Site)
  
  colnames(spei_cat_k)[colnames(spei_cat_k) == "value"] = "SPEI_value"
  colnames(spei_cat_k)[colnames(spei_cat_k) == "drought_category"] = "Drought_category"
  
  veg_data_k=merge(veg_data,spei_cat_k%>%
                     dplyr::select(., Site,Month,Year,SPEI_value,Drought_category),
                   by=c("Month","Year","Site"))
  
  return(veg_data_k)  
})
names(Merged_data)=names(spei_cat)
saveRDS(Merged_data,"./Data/Vegetation_SPEI.rds")


#add temperature
Merged_data=readRDS("./Data/Vegetation_SPEI.rds")
Info_plots=readRDS("./Data/Info_plots_thresholds_MAT.rds")
All_clim_temperature=read.table("./Data/plots_SM.csv",sep=",",header = T)

Merged_data=lapply(Merged_data, function(x){
  
  x=x%>%merge(.,All_clim_temperature[,c("Temperature","plotID","Year","Month")],
                   by.x=c("Month","Year","Site"),by.y=c("Month","Year","plotID"))
  
  Info_plots$Threshold_MAT$Site=str_c(
    str_extract(Info_plots$Threshold_MAT$Site, "^[A-Z]+"),  
    str_pad(str_extract(Info_plots$Threshold_MAT$Site, "\\d+"), width = 2, pad = "0"))
  x=merge(x,Info_plots$Threshold_MAT,by="Site")
  
  x$Temperature_category = with(x, 
                                 ifelse(Temperature <= Extreme_threshold_cold, "Extreme cold",
                                        ifelse(Temperature <= Moderate_threshold_cold, "Moderate cold",
                                               ifelse(Temperature <= Moderate_threshold_hot, "Non-temperature",
                                                      ifelse(Temperature <= Extreme_threshold_hot, "Moderate hot","Extreme got")))))
  return(x%>%dplyr::select(., -Moderate_threshold_cold,-Extreme_threshold_cold,
                           -Moderate_threshold_hot,-Extreme_threshold_hot))
})
saveRDS(Merged_data,"./Data/Vegetation_SPEI.rds")


## >> 3) Compute resilience and resistance and recovery and constancy (SPEI) ----

#stability indices using community biomass
Merged_data=readRDS("./Data/Vegetation_SPEI.rds")%>%
  lapply(.,function(spei_k){
    
    data_k=spei_k%>%
      dplyr::mutate(., Total_biomass=as.numeric(Total_biomass))%>%
      dplyr::group_by(Site) %>%  
      dplyr::arrange(Year, .by_group = TRUE) %>%
      dplyr::mutate(Total_biomass = zoo::na.approx(Total_biomass, method = "linear", na.rm = FALSE)) %>% #in case of missing values
      dplyr::ungroup(.)
    
    Temporal_stab=data_k%>%  #getting temporal stability
      dplyr::group_by(., Site)%>%
      dplyr::summarise(., .groups = "keep",Temporal_stability=sd(Total_biomass,na.rm = T)/mean(Total_biomass,na.rm=T))
    
    Temporal_stab_detrended =data_k %>% #same but detrending the biomass to avoid bias on the SD
      dplyr::group_by(Site) %>%
      dplyr::arrange(Year, .by_group = TRUE) %>%
      dplyr::mutate(time = dplyr::row_number(),  
                    biomass_detrended_loess = Total_biomass - predict(loess(Total_biomass ~ time, na.action = na.exclude)),
                    biomass_detrended_lm = residuals(lm(Total_biomass ~ time, na.action = na.exclude))) %>%
      dplyr::summarise(
        .groups = "keep",
        Temporal_stability_loess = sd(biomass_detrended_loess, na.rm = TRUE)/mean(Total_biomass, na.rm = TRUE),
        Temporal_stability_lm = sd(biomass_detrended_lm, na.rm = TRUE)/mean(Total_biomass, na.rm = TRUE)
      )
    
    #as in JOE paper 2024 on explo 10.1111/1365-2745.14288, though its hard to understand what does two log ratio correspond to ecologically
    LogRatios_t = data_k%>% 
      dplyr::arrange(Site, Year)%>%
      dplyr::mutate(.,Total_biomass_lag1 = lag(Total_biomass))%>%
      dplyr::mutate(.,LogR = log(Total_biomass / Total_biomass_lag1))
    
    LogRatios_ref = data_k%>% #kind of a resistance definition
      dplyr::arrange(Site, Year)%>%
      dplyr::group_by(., Site)%>%
      dplyr::mutate(., Median_biomass = median(Total_biomass, na.rm = TRUE))%>%
      dplyr::mutate(.,LogR_ref_plot = log(Total_biomass / Median_biomass))
    
    #Average under non-drought conditions
    Normal_means = data_k %>%
      dplyr::filter(Drought_category == "Non-drought") %>%
      dplyr::group_by(Site) %>%
      dplyr::summarise(Y_normal_mean = mean(Total_biomass, na.rm = TRUE))
    
    #RESISTANCE DROUGHT AND WET EVENTS
    Resistance_Isbell = data_k %>%
      dplyr::left_join(Normal_means, by = "Site") %>%
      dplyr::mutate(.,
                    Resistance_Isbell=case_when(Drought_category %in% c("Moderate drought","Extreme drought") ~ Y_normal_mean / (Y_normal_mean-Total_biomass),
                                                Drought_category %in% c("Moderate wet","Extreme wet") ~ Y_normal_mean / (Total_biomass-Y_normal_mean),
                                                Drought_category =="Non-drought" ~ NA))%>% #we don't take absolute value as we wanna focus on drought that have impact on productivity (negative ones)
      dplyr::mutate(.,Resistance_Isbell_abs=Y_normal_mean / abs(Y_normal_mean-Total_biomass))
    
    
    #RESILIENCE DROUGHT AND WET EVENTS
    Recovery = data_k %>%
      dplyr::left_join(Normal_means, by = "Site")%>%
      dplyr::group_by(Site) %>%
      dplyr::arrange(Site) %>%
      dplyr::mutate(
        is_drought = Drought_category %in% c("Moderate drought","Extreme drought"),
        Y_next = lead(Total_biomass),
        drought_next = lead(as.character(Drought_category)),
        drought_previous = lag(as.character(Drought_category)))%>%
      dplyr::mutate(Recovery = (Y_normal_mean / abs((Y_next - Y_normal_mean))))
    
    #RESILIENCE DROUGHT AND WET EVENTS
    Resilience_Isbell = data_k %>%
      dplyr::left_join(Normal_means, by = "Site")%>%
      dplyr::group_by(Site) %>%
      dplyr::arrange(Site) %>%
      dplyr::mutate(
        is_drought = Drought_category %in% c("Moderate drought","Extreme drought"),
        Y_next = lead(Total_biomass),
        drought_next = lead(as.character(Drought_category)),
        drought_previous = lag(as.character(Drought_category)))
    
    #to make sure
    Resilience_Isbell$Y_next[which(is.na(Resilience_Isbell$Y_next))]=Resilience_Isbell$Total_biomass[which(is.na(Resilience_Isbell$Y_next))+1]
    Resilience_Isbell$drought_next[which(is.na(Resilience_Isbell$drought_next))]=Resilience_Isbell$Drought_category[which(is.na(Resilience_Isbell$drought_next))+1]
    
    Resilience_Isbell=Resilience_Isbell%>%
      dplyr::mutate(Resilience_Isbell = ((Total_biomass - Y_normal_mean) / (Y_next - Y_normal_mean)))%>%
      dplyr::mutate(Resilience_Isbell_abs = abs((Total_biomass - Y_normal_mean) / (Y_next - Y_normal_mean)))%>%
      dplyr::mutate(.,is_drought_next=drought_next %in% c("Moderate drought","Extreme drought"))%>%
      dplyr::mutate(.,Recovery_context = case_when(
        is_drought_next ~ "Consecutive_drought",
        isFALSE(is_drought_next) ~ "Non_consecutive_drought"))
    
    #merging everything
    
    d_all=merge(LogRatios_t%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Total_biomass,Total_biomass_lag1,LogR,Drought_category,SPEI_value),
                LogRatios_ref%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,LogR_ref_plot),by=c("Site","Year","Month"))%>%
      merge(Resilience_Isbell%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Recovery_context,Resilience_Isbell,Resilience_Isbell_abs,drought_next,drought_previous),by=c("Site","Year","Month"))%>%
      merge(Resistance_Isbell%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Resistance_Isbell,Resistance_Isbell_abs,Sampling_day),by=c("Site","Year","Month"))%>%
      merge(Recovery%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Recovery),by=c("Site","Year","Month"))%>%
      merge(Temporal_stab%>%ungroup(.),by=c("Site"))%>%
      merge(Temporal_stab_detrended%>%ungroup(.),by=c("Site"))%>%
      dplyr::rename(.,Drought_category_next=drought_next,Drought_category_previous=drought_previous,Logbiomass_t_minus_1=LogR,Logbiomass_t_median=LogR_ref_plot)%>%
      dplyr::select(., 
                    Month,Year,Site,Sampling_day,#site and years
                    Total_biomass,Total_biomass_lag1, #biomass
                    SPEI_value,Drought_category,Drought_category_next,Drought_category_previous, #SPEI related
                    Temporal_stability,Temporal_stability_loess,Temporal_stability_lm, #temporalstability
                    Logbiomass_t_median,Logbiomass_t_minus_1, #other metrics
                    Resistance_Isbell, Resistance_Isbell_abs, #resistance
                    Resilience_Isbell,Resilience_Isbell_abs,Recovery,Recovery_context#resilience and is there a drought after (to avoid legacy effects)
      )
    return(d_all)
  })


#stability indices using community cover
Merged_data_cover=readRDS("./Data/Vegetation_SPEI.rds")%>%
  lapply(.,function(spei_k){
    
    data_k=spei_k%>%
      dplyr::mutate(., Total_biomass=as.numeric(Cover))%>% #CHANGING BIOMASS TO COVER
      dplyr::mutate(.,    Total_biomass = ifelse(Total_biomass == 0, NA, Total_biomass))%>%
      dplyr::group_by(Site) %>%  
      dplyr::arrange(Year, .by_group = TRUE) %>%
      dplyr::mutate(Total_biomass = zoo::na.approx(Total_biomass, method = "linear", na.rm = FALSE)) %>%
      dplyr::ungroup(.)
    
    Temporal_stab=data_k%>%  #getting temporal stability
      dplyr::group_by(., Site)%>%
      dplyr::summarise(., .groups = "keep",Temporal_stability=mean(Total_biomass,na.rm=T)/sd(Total_biomass,na.rm = T))
    
    Temporal_stab_detrended =data_k %>% #same but detrending the biomass to avoid bias on the SD
      dplyr::group_by(Site) %>%
      dplyr::arrange(Year, .by_group = TRUE) %>%
      dplyr::mutate(time = dplyr::row_number(),  
                    biomass_detrended_loess = Total_biomass - predict(loess(Total_biomass ~ time, na.action = na.exclude)),
                    biomass_detrended_lm = residuals(lm(Total_biomass ~ time, na.action = na.exclude))) %>%
      dplyr::summarise(
        .groups = "keep",
        Temporal_stability_loess = mean(Total_biomass, na.rm = TRUE) / sd(biomass_detrended_loess, na.rm = TRUE),
        Temporal_stability_lm = mean(Total_biomass, na.rm = TRUE) / sd(biomass_detrended_lm, na.rm = TRUE)
      )
    
    
    #as in JOE paper 2024 on explo 10.1111/1365-2745.14288, though its hard to understand what does two log ratio correspond to ecologically
    LogRatios_t = data_k%>% 
      dplyr::arrange(Site, Year)%>%
      dplyr::mutate(.,Total_biomass_lag1 = lag(Total_biomass))%>%
      dplyr::mutate(.,LogR = log(Total_biomass / Total_biomass_lag1))
    
    LogRatios_ref = data_k%>% #kind of a resistance definition
      dplyr::arrange(Site, Year)%>%
      dplyr::group_by(., Site)%>%
      dplyr::mutate(., Median_biomass = median(Total_biomass, na.rm = TRUE))%>%
      dplyr::mutate(.,LogR_ref_plot = log(Total_biomass / Median_biomass))
    
    
    #Average under non-drought conditions
    Normal_means = data_k %>%
      dplyr::filter(Drought_category == "Non-drought") %>%
      dplyr::group_by(Site) %>%
      dplyr::summarise(Y_normal_mean = mean(Total_biomass, na.rm = TRUE))
    
    
    
    #RESISTANCE DROUGHT AND WET EVENTS
    Resistance_Isbell = data_k %>%
      dplyr::left_join(Normal_means, by = "Site") %>%
      dplyr::mutate(.,
                    Resistance_Isbell=case_when(Drought_category %in% c("Moderate drought","Extreme drought") ~ Y_normal_mean / (Y_normal_mean-Total_biomass),
                                                Drought_category %in% c("Moderate wet","Extreme wet") ~ Y_normal_mean / (Total_biomass-Y_normal_mean),
                                                Drought_category =="Non-drought" ~ NA))%>% #we don't take absolute value as we wanna focus on drought that have impact on productivity (negative ones)
      dplyr::mutate(.,Resistance_Isbell_abs=Y_normal_mean / abs(Y_normal_mean-Total_biomass))
    
    
    #RESILIENCE DROUGHT AND WET EVENTS
    Recovery = data_k %>%
      dplyr::left_join(Normal_means, by = "Site")%>%
      dplyr::group_by(Site) %>%
      dplyr::arrange(Site) %>%
      dplyr::mutate(
        is_drought = Drought_category %in% c("Moderate drought","Extreme drought"),
        Y_next = lead(Total_biomass),
        drought_next = lead(as.character(Drought_category)),
        drought_previous = lag(as.character(Drought_category)))%>%
      dplyr::mutate(Recovery = (Y_normal_mean / abs((Y_next - Y_normal_mean))))
    
    #RESILIENCE DROUGHT AND WET EVENTS
    Resilience_Isbell = data_k %>%
      dplyr::left_join(Normal_means, by = "Site")%>%
      dplyr::group_by(Site) %>%
      dplyr::arrange(Site) %>%
      dplyr::mutate(
        is_drought = Drought_category %in% c("Moderate drought","Extreme drought"),
        Y_next = lead(Total_biomass),
        drought_next = lead(as.character(Drought_category)),
        drought_previous = lag(as.character(Drought_category)))
    
    #to make sure
    Resilience_Isbell$Y_next[which(is.na(Resilience_Isbell$Y_next))]=Resilience_Isbell$Total_biomass[which(is.na(Resilience_Isbell$Y_next))+1]
    Resilience_Isbell$drought_next[which(is.na(Resilience_Isbell$drought_next))]=Resilience_Isbell$Drought_category[which(is.na(Resilience_Isbell$drought_next))+1]
    
    Resilience_Isbell=Resilience_Isbell%>%
      dplyr::mutate(Resilience_Isbell = ((Total_biomass - Y_normal_mean) / (Y_next - Y_normal_mean)))%>%
      dplyr::mutate(Resilience_Isbell_abs = abs((Total_biomass - Y_normal_mean) / (Y_next - Y_normal_mean)))%>%
      dplyr::mutate(.,is_drought_next=drought_next %in% c("Moderate drought","Extreme drought"))%>%
      dplyr::mutate(.,Recovery_context = case_when(
        is_drought_next ~ "Consecutive_drought",
        isFALSE(is_drought_next) ~ "Non_consecutive_drought"))
    
    #merging everything
    
    d_all=merge(LogRatios_t%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Total_biomass,Total_biomass_lag1,LogR,Drought_category,SPEI_value),
                LogRatios_ref%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,LogR_ref_plot),by=c("Site","Year","Month"))%>%
      merge(Resilience_Isbell%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Recovery_context,Resilience_Isbell,Resilience_Isbell_abs,drought_next,drought_previous),by=c("Site","Year","Month"))%>%
      merge(Resistance_Isbell%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Resistance_Isbell,Resistance_Isbell_abs,Sampling_day),by=c("Site","Year","Month"))%>%
      merge(Recovery%>%ungroup(.)%>%dplyr::select(.,Month,Year,Site,Recovery),by=c("Site","Year","Month"))%>%
      merge(Temporal_stab%>%ungroup(.),by=c("Site"))%>%
      merge(Temporal_stab_detrended%>%ungroup(.),by=c("Site"))%>%
      dplyr::rename(.,Drought_category_next=drought_next,Drought_category_previous=drought_previous,Logbiomass_t_minus_1=LogR,Logbiomass_t_median=LogR_ref_plot)%>%
      dplyr::select(., 
                    Month,Year,Site,Sampling_day,#site and years
                    Total_biomass,Total_biomass_lag1, #biomass
                    SPEI_value,Drought_category,Drought_category_next,Drought_category_previous, #SPEI related
                    Temporal_stability,Temporal_stability_loess,Temporal_stability_lm, #temporalstability
                    Logbiomass_t_median,Logbiomass_t_minus_1, #other metrics
                    Resistance_Isbell, Resistance_Isbell_abs, #resistance
                    Resilience_Isbell,Resilience_Isbell_abs,Recovery,Recovery_context#resilience and is there a drought after (to avoid legacy effects)
      )
    return(d_all)
  })

saveRDS(list(Cover=Merged_data_cover,Biomass=Merged_data),"./Data/Stability_indices.rds")


# ---------------------- Step 2: Phylogenetic diversity ----

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
  dplyr::mutate(., Species=gsub("_aggr.","",Species))

Phylo_explo = readRDS("./Data/Phylo_explo.rds")
Phylo_explo=ape::as.phylo(Phylo_explo$scenario.3)

matrix_data = veg_data %>%
  dplyr::group_by(., Year,Site,Species)%>%
  dplyr::summarise(., Cover=sum(Cover),.groups = "drop")%>%
  dplyr::mutate(., Com=paste0(Site,"_",Year))%>%
  dplyr::select(.,Cover,Species,Com)%>%
  dplyr::filter(.,Species %in% Phylo_explo$tip.label)%>%
  pivot_wider(names_from = Species, values_from = Cover,values_fill = list(Cover=0))%>%
  as.data.frame(.)

rownames(matrix_data)=matrix_data$Com
matrix_data=as.matrix(matrix_data[,-1])
matrix_data[is.na(matrix_data)]=0

species_matrix = colnames(matrix_data)
species_phylo = Phylo_explo$tip.label
setdiff(species_matrix, species_phylo)  # Species in matrix not in phylogeny: OKKK all good
setdiff(species_phylo, species_matrix)  # Species in phylogeny not in matrix : OKKK all good
tree_pruned=ape::drop.tip(Phylo_explo, setdiff(species_phylo, species_matrix))

psd_explo=picante::psd(matrix_data, vcv.phylo(tree_pruned, corr = TRUE))
mntd_explo_weighted=mntd(matrix_data, cophenetic(tree_pruned), abundance.weighted=TRUE)
mntd_explo_non_weighted=mntd(matrix_data, cophenetic(tree_pruned), abundance.weighted=F)
mpd_explo_weighted=mpd(matrix_data, cophenetic(tree_pruned), abundance.weighted=T)
mpd_explo_non_weighted=mpd(matrix_data, cophenetic(tree_pruned), abundance.weighted=F)
faith_pd=pd(matrix_data,tree_pruned,include.root = F)

if (is.ultrametric(tree_pruned)) Rao_pd=raoD(matrix_data,tree_pruned) #Getting Rao PD

Phylo_div_metrics=list(PSD=psd_explo%>%add_column(., Complete_name=rownames(psd_explo)),
                       MNTD_weighted=data.frame(MNTD_weighted=mntd_explo_weighted,Complete_name=rownames(matrix_data)),
                       MNTD_non_weighted=data.frame(MNTD_nonweighted=mntd_explo_non_weighted,Complete_name=rownames(matrix_data)),
                       MPD_weighted=data.frame(MPD_weighted=mpd_explo_weighted,Complete_name=rownames(matrix_data)),
                       MPD_no_weighted=data.frame(MPD_nonweighted=mpd_explo_non_weighted,Complete_name=rownames(matrix_data)),
                       Faith_FD=data.frame(Faith_FD=faith_pd$PD,Complete_name=rownames(matrix_data)))

saveRDS(Phylo_div_metrics,"./Data/Phylogenetic_diversity.rds")


# ---------------------- Step 3: Traits ----
## >> 1) TRY TRAITS dataset ----
# ------ representativity of traits 

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
  dplyr::summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop")

trait_sp=read.table("./Data/Trait/Trait_exploratories.csv",sep=",",header = T)%>% #TRY
  dplyr::mutate(
    scientificName = case_when(
      scientificName == "Capsella_bursa_pastoris" ~ "Capsella_bursa-pastoris",
      scientificName == "Capsella_bursa.pastoris" ~ "Capsella_bursa-pastoris",
      scientificName == "Galeopsis_bifida" ~ "Galeopsis_cf_bifida",
      scientificName == "Geum_rivale" ~ "Geum_rivale/urbanum_aggr.",
      scientificName == "Persicaria_amphibia" ~ "Persicaria_amphibium",
      scientificName == "Primula_elatior_veris_aggr." ~ "Primula_elatior/veris_aggr.",
      scientificName == "Persicaria_lapathifolia" ~ "Ploygonum_lapathifolium",
      scientificName == "Tripleurospermum_perforatum" ~ "Tripleurospermum_inodorum",
      scientificName == "Trifolium_campestre" ~ "Trifolium_campestre/dubium_aggr.",
      scientificName == "Viola_alba" ~ "Viola_cf_alba",
      scientificName == "Viola_tricolor" ~ "Viola_tricolor_aggr.",
      scientificName == "Lythrum_salicaria" ~ "cf_Lythrum_salcaria",
      scientificName == "Ononis_repens_spinosa_aggr." ~ "Ononis_repens/spinosa_aggr.",
      TRUE ~ scientificName))%>%
  dplyr::filter(., scientificName %in% veg_data$Species)

unique(veg_data$Species)[which(unique(veg_data$Species) %!in% trait_sp$scientificName)]

#looking at species that have no recorded traits (0.7% of total cover)
100*sum(veg_data$Cover[which(veg_data$Species %in% unique(veg_data$Species)
                             [which(!unique(veg_data$Species) %in% trait_sp$scientificName)])],na.rm = T)/
  sum(veg_data$Cover,na.rm = T)

all_traits=lapply(unique(trait_sp$traitName),function(trait_k){
  cover_by_site = veg_data %>%
    dplyr::mutate(has_trait = Species %in% trait_sp$scientificName[which(trait_sp$verbatimTraitName==trait_k)],
                  no_trait = !has_trait) %>%
    dplyr::group_by(Site,Year) %>%
    dplyr::summarise(
      total_cover = sum(Cover, na.rm = TRUE),
      no_trait_cover = sum(Cover[no_trait], na.rm = TRUE),
      percent_no_trait = 100 * (no_trait_cover / total_cover),.groups = 'keep'
    )
})
names(all_traits)=unique(trait_sp$traitName)

pdf("./Figures/Trait_representativity_aboveground_TRY.pdf",height = 12,width = 10)
par(mfrow=c(4,3))
for (k in 1:10){
  hist(all_traits[[k]]$percent_no_trait,title="",xlab="",ylab="",main="")
  abline(v = 30,lwd=4,col="blue") 
  mtext(paste0("Trait = ",unique(trait_sp$traitName)[k]))
}
dev.off()

trait_sp=trait_sp%>%
  dplyr::filter(., scientificName%in%veg_data$Species)
saveRDS(trait_sp,"./Data/Trait/TRAIT_try_cleaned.rds")


## >> 2) Hairphae ----

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
  dplyr::summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop")

d=read.table("./Data/Trait/hairphae_species_traits_inclrootN_notcleaned.csv",sep = ",",header = T)
unique(d$species_bexis)[which(unique(d$species_bexis) %!in% veg_data$Species)] #one species is not found

d=d%>%reshape2::melt(.,id.vars=c("species_bexis"))

all_traits=lapply(unique(d$variable)[2:8],function(trait_k){
  cover_by_site = veg_data %>%
    dplyr::mutate(has_trait = Species %in% d$species_bexis[which(d$variable==trait_k)],
                  no_trait = !has_trait) %>%
    dplyr::group_by(Site,Year) %>%
    dplyr::summarise(
      total_cover = sum(Cover, na.rm = TRUE),
      no_trait_cover = sum(Cover[no_trait], na.rm = TRUE),
      percent_no_trait = 100 * (no_trait_cover / total_cover),.groups = 'keep'
    )
})
names(all_traits)=unique(d$variable)[2:8]

pdf("./Figures/Trait_representativity_belowground_Hairphae.pdf",height = 6,width = 14)
par(mfrow=c(2,4))
for (k in 1:7) {
  hist(all_traits[[k]]$percent_no_trait,title="",xlab="",ylab="",main="")
  abline(v = 30,lwd=4,col="blue") 
  mtext(paste0("Trait = ",unique(d$variable)[2:8][k]))
}
dev.off()

unlist(lapply(all_traits, function(x){
  sum(x$percent_no_trait<=30,na.rm = T)/nrow(x)
}))



## >> 3) Merging trait datasets ----

#Merging trait datasets into a global trait dataset

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

Hairphae = read.table("./Data/Trait/hairphae_species_traits_inclrootN_notcleaned.csv", sep=",", header=T) %>%
  dplyr::select(-species_matched, -X, -X.1, -X.2, -X.3) %>%
  dplyr::rename(Species = species_bexis) %>%
  dplyr::mutate(Species = gsub("_aggr.", "", Species))%>%
  dplyr::filter(., !is.na(Species))

TRY_explo = readRDS("./Data/Trait/TRAIT_try_cleaned.rds") %>%
  dplyr::select(scientificName, traitValue, traitName) %>%
  dplyr::filter(traitName %!in% c("SLA_with_petiole", "Myco_type", "SSD")) %>%
  dplyr::rename(Species = scientificName, Value = traitValue, Name = traitName) %>%
  dplyr::mutate(Value = as.numeric(Value)) %>%
  pivot_wider(id_cols = Species, names_from = Name, values_from = Value, values_fn = mean) %>%
  dplyr::mutate(Species = gsub("_aggr.", "", Species))

#all good
species_Hairphae = unique(Hairphae$Species)
species_TRY = unique(TRY_explo$Species)
all_species = unique(c(species_Hairphae, species_TRY))

all_species[which(all_species %!in% veg_data$Species)] #all species are in the cover daata, good
unique(veg_data$Species)[which(unique(veg_data$Species) %!in% all_species)] #all species are in the cover daata, good

colnames(Hairphae)[2:8]=paste0(colnames(Hairphae)[2:8],"_HP")
traits_Hairphae = setdiff(colnames(Hairphae), "Species")
traits_TRY = setdiff(colnames(TRY_explo), "Species")

all_traits = unique(c(traits_Hairphae, traits_TRY))
global_traits = data.frame(Species = all_species,stringsAsFactors = FALSE)

for(trait in all_traits) global_traits[[trait]] = NA_real_

merge_trait_source = function(global_df, 
                              source_df, 
                              source_name) {
  
  #rownmaes as species
  source_df = as.data.frame(source_df)
  rownames(source_df) = source_df$Species
  source_df=source_df%>%dplyr::select(., -Species)
  
  traits = colnames(source_df)
  
  # We fill species values for each trait
  for(i in 1:nrow(global_df)) { #i = species
    sp = global_df$Species[i]
    
    if(sp %in% rownames(source_df)) {
      for(trait in traits) {
        if(!is.na(source_df[sp, trait])) {
          global_df[i, trait] = source_df[sp, trait]
        }
      }
    }
  }
  
  return(global_df)
}

global_traits = merge_trait_source(global_traits, Hairphae, "Hairphae")
global_traits = merge_trait_source(global_traits, TRY_explo, "TRY")

global_traits=as.data.frame(global_traits)
rownames(global_traits)=global_traits$Species
global_traits=global_traits%>%dplyr::select(., -Species)
write.table(global_traits,"./Data/Global_trait_dataset.csv",sep=";")


#Checking representativity

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

global_traits = read.table("./Data/Global_trait_dataset.csv", sep=";")

all_traits = lapply(colnames(global_traits), function(trait_k) {
  species_with_trait = rownames(global_traits)[which(!is.na(global_traits[[trait_k]]))]
  
  veg_data %>%
    dplyr::mutate(
      has_trait  = Species %in% species_with_trait,
      no_trait   = !has_trait) %>%
    dplyr::group_by(Site, Year) %>%
    dplyr::summarise(
      total_cover             = sum(Cover, na.rm = TRUE),
      trait_cover             = sum(Cover[has_trait],  na.rm = TRUE),
      no_trait_cover          = sum(Cover[no_trait],   na.rm = TRUE),
      percent_cover_with_trait = 100 * (trait_cover / total_cover),
      .groups = "drop") %>%
    dplyr::mutate(trait = trait_k)
})
names(all_traits) = colnames(global_traits)
all_traits_df = dplyr::bind_rows(all_traits)%>%
  dplyr::group_by(., trait) %>%  
  dplyr::summarise(
    count_above = sum(percent_cover_with_trait >= 70, na.rm = TRUE),
    fraction_above = sum(percent_cover_with_trait >= 70, na.rm = TRUE)/n(),  )

## >> 4) PCA, and community CWM and diversity ----
## RDepth, RDia, RTD, RN, SRL (underplot)
# ------ PCA traits 

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

global_traits=read.table("./Data/Global_trait_dataset.csv",sep=";")%>%
  Closer_to_normal_traits(.)%>%
  dplyr::select(., all_of(c(colnames(.)[grep("_HP",colnames(.))],"Seed_mass",
                            # "Myco_intensity",
                            "SLA_all","Height","LDMC","LeafN","LeafP")))%>%
  drop_na(.)

all_traits =  veg_data %>%
  dplyr::mutate(
    has_trait  = Species %in% rownames(global_traits),
    no_trait   = !has_trait) %>%
  dplyr::group_by(Site, Year) %>%
  dplyr::summarise(
    total_cover             = sum(Cover, na.rm = TRUE),
    trait_cover             = sum(Cover[has_trait],  na.rm = TRUE),
    no_trait_cover          = sum(Cover[no_trait],   na.rm = TRUE),
    percent_cover_with_trait = 100 * (trait_cover / total_cover),
    .groups = "drop") 

dplyr::bind_rows(all_traits)%>%
  dplyr::summarise(
    count_above = sum(percent_cover_with_trait >= 70, na.rm = TRUE),
    fraction_above = sum(percent_cover_with_trait >= 70, na.rm = TRUE)/n())

colnames(global_traits)=c("SRL","RTD","AD","RN","RHL",'AMF','RHI',"SM","SLA","Height","LDMC","LNC","LPC")

paranResult = paran::paran(global_traits,
                           iterations = 5000, centile = 95, quietly = TRUE, status = FALSE)

paranResult$Retained
fitAll = psych::principal(r = global_traits,nfactors = 3,rotate = "varimax")
p1=plotPCA(fitAll$scores,
           fitAll$loadings,
           fitAll$Vaccounted,annotateFactor = 3,
           xIndex = 1,yIndex = 2,
           xLim = c(-3,3),yLim = c(-3,3))

p2=plotPCA(fitAll$scores,
           fitAll$loadings,
           fitAll$Vaccounted,annotateFactor = 3,
           xIndex = 1,yIndex = 3,
           xLim = c(-3,3),yLim = c(-3,3))

ggsave("./Figures/PCA_traits_above_below_Hairphae.pdf",
       ggarrange(p1,p2),width = 10,height = 4)


global_traits$PC1=fitAll$scores[,"RC1"]
global_traits$PC2=fitAll$scores[,"RC2"]
global_traits$PC3=fitAll$scores[,"RC3"]
write.table(global_traits,"./Data/Global_trait_dataset_Hairphae.csv",sep=";")


# ------ Trait diversity and CWM (above + belowground dAta)  


Trait_communities = read.table("./Data/Vegetation/vegetation_species.csv", sep=",", header = T) %>%
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
  dplyr::mutate(Species = gsub("_aggr.", "", Species))%>%
  dplyr::filter(., !is.na(Cover))%>%
  dplyr::mutate(.,dataset="Hairphae")%>%
  dplyr::group_by(., Year,Site)%>%
  dplyr::do(., Compute_functional_structure(.$Cover,.$Species,.$dataset,.$ID_year_plot))

write.table(Trait_communities,"./Data/Trait/Community_diversity_HairPhae.csv",sep=";")



# ---------------------- Step 4: Merging final data set (with LUI) ----
## >> 1) PCA traits ----

timescale=readRDS("./Data/Timescale_SPEI.rds")
stability=readRDS("./Data/Stability_indices.rds")
Traits=read.table("./Data/Trait/Community_diversity_HairPhae.csv",sep=";")
Phylo_div=readRDS("./Data/Phylogenetic_diversity.rds")
Geo_clim=read.table("./Data/Geographical_climatic.csv",sep=";")
Explo_clim=read.table("./Data/plots_SM.csv",sep=",",header = T)%>%
  dplyr::group_by(plotID,Month) %>%
  dplyr::mutate(
    Soil_moisture = ifelse(is.na(Soil_moisture), mean(Soil_moisture, na.rm = TRUE), Soil_moisture),
    Temperature = ifelse(is.na(Temperature), mean(Temperature, na.rm = TRUE), Temperature),
    Precipitation = ifelse(is.na(Precipitation), mean(Precipitation, na.rm = TRUE), Precipitation)) %>%
  ungroup(.)

LUI_regional=read.table("./Data/LUI/LUI_default components set_regional_separately_2026-06-03.txt",sep=",",header = T)%>%
  dplyr::arrange(., PLOTID)%>%
  dplyr::filter(., YEAR != "separately(2025)")%>%
  dplyr::mutate(., YEAR=as.numeric(gsub("separately","",gsub(paste(c("[(]", "[)]"), collapse = "|"), "", .$YEAR))))
LUI_global=read.table("./Data/LUI/LUI_default components set_global_separately_2026-03-15.txt",sep=",",header = T)%>%
  dplyr::arrange(., PLOTID)%>%
  dplyr::mutate(., YEAR=as.numeric(gsub("separately","",gsub(paste(c("[(]", "[)]"), collapse = "|"), "", .$YEAR))))
LUI_detail=read.table("./Data/LUI/LUI_detail.csv",sep=",",header = T)%>%
  dplyr::arrange(., EP_PlotID)%>%
  dplyr::mutate(.,
                TotalFertilization=(TotalFertilization/mean(TotalFertilization,na.rm=T)),
                TotalMowing=(TotalMowing/mean(TotalMowing,na.rm=T)),
                TotalGrazing=(TotalGrazing/mean(TotalGrazing,na.rm=T)))


LUI_global$PLOTID=str_c(
  str_extract(LUI_global$PLOTID, "^[A-Z]+"),  
  str_pad(str_extract(LUI_global$PLOTID, "\\d+"), width = 2, pad = "0"))
LUI_regional$PLOTID=str_c(
  str_extract(LUI_regional$PLOTID, "^[A-Z]+"),  
  str_pad(str_extract(LUI_regional$PLOTID, "\\d+"), width = 2, pad = "0"))

LUI_regional=read.table("./Data/LUI/LUI_default components set_regional_separately_2026-06-03.txt",sep=",",header = T)%>%
  dplyr::arrange(., PLOTID)%>%
  dplyr::filter(., YEAR != "separately(2025)")%>%
  dplyr::mutate(., YEAR=as.numeric(gsub("separately","",gsub(paste(c("[(]", "[)]"), collapse = "|"), "", .$YEAR))))
LUI_global=read.table("./Data/LUI/LUI_default components set_global_separately_2026-03-15.txt",sep=",",header = T)%>%
  dplyr::arrange(., PLOTID)%>%
  dplyr::mutate(., YEAR=as.numeric(gsub("separately","",gsub(paste(c("[(]", "[)]"), collapse = "|"), "", .$YEAR))))

LUI_global$PLOTID=str_c(
  str_extract(LUI_global$PLOTID, "^[A-Z]+"),  
  str_pad(str_extract(LUI_global$PLOTID, "\\d+"), width = 2, pad = "0"))
LUI_regional$PLOTID=str_c(
  str_extract(LUI_regional$PLOTID, "^[A-Z]+"),  
  str_pad(str_extract(LUI_regional$PLOTID, "\\d+"), width = 2, pad = "0"))

Merged_dataset=lapply(1:24, function(timescale_k){
  
  final_df=stability$Biomass[[timescale_k]]%>%
    merge(., Traits,by=c("Year","Site"),all.x=T)%>%
    merge(., Geo_clim,by.x=c("Site"),by.y=c("site_name"))%>%
    merge(.,Explo_clim%>%dplyr::select(., -datetime),by.x=c("Site","Year","Month"),by.y=c("plotID","Year","Month"))%>%
    merge(., timescale$Best_timescale_biomass%>%
            dplyr::select(., Site,Correlation_biomass)%>%
            dplyr::rename(.,Sensitivity_water_biomass=Correlation_biomass),by="Site",all.x=T)%>%
    merge(., timescale$Best_timescale_cover%>%
            dplyr::select(., Site,Correlation_biomass)%>%
            dplyr::rename(., Sensitivity_water_cover=Correlation_biomass),by="Site",all.x=T)%>%
    merge(., LUI_global%>%
            dplyr::select(., PLOTID,YEAR,LUI)%>%
            dplyr::rename(., LUI_global=LUI),by.x=c("Site","Year"),by.y=c("PLOTID","YEAR"))%>%
    merge(., LUI_regional%>%
            dplyr::select(., PLOTID,YEAR,LUI)%>%
            dplyr::rename(., LUI_regional=LUI),by.x=c("Site","Year"),by.y=c("PLOTID","YEAR"))%>%
    merge(., LUI_detail%>%
            dplyr::select(., ,Year,EP_PlotID,TotalFertilization,TotalMowing,TotalGrazing),
          by.x=c("Site","Year"),by.y=c("EP_PlotID","Year"),all.x=T)%>%
    dplyr::mutate(., Complete_name=paste0(Site,"_",Year))%>%
    merge(., Phylo_div$PSD%>%
            dplyr::select(., PSV,PSE,PSR,PSC,Complete_name),by.x=c("Complete_name"),by.y=c("Complete_name"))%>%
    merge(., Phylo_div$MNTD_weighted,by.x=c("Complete_name"),by.y=c("Complete_name"))%>%
    merge(., Phylo_div$MNTD_non_weighted,by.x=c("Complete_name"),by.y=c("Complete_name"))%>%
    merge(., Phylo_div$MPD_weighted,by.x=c("Complete_name"),by.y=c("Complete_name"))%>%
    merge(., Phylo_div$MPD_no_weighted,by.x=c("Complete_name"),by.y=c("Complete_name"))%>%
    merge(., Phylo_div$Faith_FD,by.x=c("Complete_name"),by.y=c("Complete_name"))%>%
    dplyr::mutate(., 
                  RaoQ_PC1=sqrt(RaoQ_PC1),
                  RaoQ_PC2=sqrt(RaoQ_PC2),
                  RaoQ_PC3=sqrt(RaoQ_PC3))%>%
    dplyr::mutate(., 
                  MNTD_weighted=log(MNTD_weighted),
                  FD_PC1=log(FD_PC1),
                  FD_PC2=log(FD_PC2),
                  FD_PC3=log(FD_PC3))
  
  
  return(final_df)
})
names(Merged_dataset)=names(stability$Cover)
saveRDS(Merged_dataset,"./Data/Merged_datasets_Hairphae.rds")

