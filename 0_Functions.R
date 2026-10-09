x=c("sf","raster","elevatr","ggplot2","ggpubr","rnaturalearth","V.PhyloMaker",
    "rnaturalearthhires","ggspatial","reshape2","tidyverse","FD","car","funrar",
    "missForest","ape","lme4","LMERConvenienceFunctions","lmerTest","lmeresampler",
    "MuMIn","boot","nlme","boot","visreg","arm","purrr","DHARMa","spdep","nlme","mgcv","zoo",
    "ggraph","psych","paran")
lapply(x, require, character.only = TRUE)
# lapply(x, install.packages, character.only = TRUE)
get_month_from_day = function(day_of_year, year) {
  date = as.Date(day_of_year - 1, origin = paste0(year, "-01-01"))
  return(as.integer(format(date, "%m")))
}
"%!in%" = Negate("%in%")

color_trait_axes=c( "#8F9DD8","#E2AA84","#82B7A4")
fill_trait_axes=c( "#C7CEEA","#FFDAC1","#B5EAD7")


## UTILS FUNCTIONS ----

Plot_correlation_variables=function(d){
  d=d[,colnames(dplyr::select_if (d, is.numeric))]
  corr_pred=corr.test(d,use = "pairwise.complete.obs",adjust = "none")
  
  corr_pred$r=round(corr_pred$r,2)
  corr_pred$r[lower.tri(corr_pred$r)]=NA
  diag(corr_pred$r)=NA
  
  p=ggplot(corr_pred$r%>%
             reshape2::melt(.)%>%
             add_column(., pval=sapply(1:nrow(.), function(x){
               if (is.na(.$value[x])){
                 return(NA)
               }else{return(reshape2::melt(corr_pred$p)$value[x])}
             })))+
    geom_tile(aes(x=Var1,Var2,fill=ifelse(pval<.05,value,0)))+
    theme_classic()+
    geom_text(aes(x=Var1,Var2,label=ifelse(pval<.05,round(value,2),"X")),size=3)+
    scale_fill_gradient2(low="red",mid="white",high = "blue",midpoint = 0,na.value = "white")+
    theme(axis.text.x = element_text(angle=60,hjust=1))+
    labs(x="",y="",fill="")
  return(p)
}

the_theme2 = theme_classic() + theme(
  legend.position = "bottom",
  strip.background = element_rect(fill = "transparent",color="transparent"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 10),title = element_text(size=8),
  axis.title.y=element_text(size = 10),
  axis.title.x=element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10), 
  text = element_text(family = "NewCenturySchoolbook")
)
Rename_PCA_axes = function(df,
                           axis_names = c(PC1 = "Economics (slow - fast)",
                                          PC2 = "Absorption (AMF - hairs)",
                                          PC3 = "Exploration (SLR - AD)"),
                           melt = FALSE,
                           predictor_col = "Predictor") {
  
  rename_one = function(x) {
    for (pc in names(axis_names)) {
      if (grepl(paste0("_", pc, "$"), x)) {
        prefix = sub(paste0("_", pc, "$"), "", x)
        return(paste0(prefix, " ", axis_names[[pc]]))
      }
    }
    if (x == "Species_richness") return("Sp. richness")
    if (x == "PSV")             return("PSV (phylo)")
    if (x == "MNTD_weighted")   return("MNTD (phylo)")
    return(x)
  }
  
  if (melt) {
    df[[predictor_col]] = vapply(df[[predictor_col]], rename_one, character(1))
  } else {
    colnames(df) = vapply(colnames(df), rename_one, character(1))
  }
  
  return(df)
}


Numeric_scaling_df=function(d){
  d[,colnames(dplyr::select_if (d, is.numeric))]=
    apply(d[,colnames(dplyr::select_if (d, is.numeric))],2,
          function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm = T))})
  
  return(d)
}

Plot_distributions_variables=function(d){
  
  p=ggplot(d[,colnames(dplyr::select_if (d, is.numeric))]%>%
             reshape2::melt(.))+
    geom_histogram(aes(x=value),fill="grey")+
    the_theme2+
    facet_wrap(.~variable,scales = "free")
  
  return(p)
}

## TRAITS FUNCTIONS ----

RaoQ_traits = function(abund, trait, Hill = TRUE, scale = FALSE, method = "default") {
  
  abund = as.matrix(abund)
  trait = as.matrix(trait)
  abund = abund / rowSums(abund)
  
  if (method == "default") {
    Q = apply(abund, 1, function(x) crossprod(x, trait %*% x))
  }
  
  if (method == "divc") {
    Q = apply(abund, 1, function(x) x %*% trait^2 %*% (x / 2 / sum(x)^2))
  }
  
  if (Hill == TRUE) Q = 1 / (1 - Q)
  
  if (scale == TRUE) Q = Q / max(Q)
  
  names(Q) = rownames(abund)
  
  return(Q)
}


plotPCA = function(fitScores, fitLoadings, fitVaccounted, xIndex, yIndex, xLim, yLim, annotateFactor = 2.6,
                   colorLow = "ivory", colorHigh = "purple", colorMid = "pink", midPoint = 0.09,
                   pointSize = 1, pointAlpha = 1, pointColor = "grey85", xTimeSegment = 2.5,
                   yTimeSegment = 2.5, sizeAnnotation = 5) {
  p=ggplot2::ggplot(data = NULL, ggplot2::aes(x = fitScores[, xIndex], y = fitScores[, yIndex])) +
    ggplot2::theme_classic() +
    ggplot2::stat_density_2d(ggplot2::aes(fill = ..level..), geom = "polygon") +
    ggplot2::geom_point(size = pointSize, alpha = pointAlpha, color = "grey85") +
    ggplot2::geom_segment(data = NULL, ggplot2::aes(x = 0, y = 0, xend = (fitLoadings[, xIndex] * xTimeSegment),
                                                    yend = (fitLoadings[, yIndex] * yTimeSegment)),
                          arrow = ggplot2::arrow(length = ggplot2::unit(1 / 2, units = "picas")),
                          color = "black") +
    ggplot2::annotate("text", label = rownames(fitLoadings), size = sizeAnnotation,
                      x = (fitLoadings[, xIndex] * annotateFactor), y = (fitLoadings[, yIndex] * annotateFactor)) +
    ggplot2::scale_fill_gradient2(low = colorLow, high = colorHigh, mid = colorMid,  midpoint = midPoint) +
    ggplot2::xlim(xLim) +
    ggplot2::ylim(yLim) +
    ggplot2::xlab(paste0("PC ", xIndex, " (",
                         round(x = fitVaccounted["Proportion Var", xIndex] * 100, digits = 1), "%)")) +
    ggplot2::ylab(paste0("PC ", yIndex, " (",
                         round(x = fitVaccounted["Proportion Var", yIndex] * 100, digits = 1), "%)")) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank()) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 14),
                   axis.title.y = ggplot2::element_text(size = 14))
  return(p)
}
Closer_to_normal_traits=function(d){
  return(d%>%
           dplyr::mutate(., 
                         SRL_HP=log(SRL_HP),
                         RTD_HP=log(RTD_HP),
                         AD_HP=log(AD_HP),
                         RN_HP=bcPower(RN_HP+.1,.9),
                         RHL_HP=sqrt(RHL_HP),
                         RHI_HP=sqrt(RHI_HP),
                         LeafP=sqrt(LeafP),
                         LeafN=log(LeafN),
                         Height=bcPower(Height,-.3),
                         Seed_mass=bcPower(Seed_mass,0),
                         LDMC=sqrt(LDMC)))
}


run_pca_scores = function(df) {
  
  pr  = paran::paran(df, iterations = 5000, centile = 95,
                     quietly = TRUE, status = FALSE)
  nfac = max(pr$Retained, 1)
  fit = psych::principal(r = df, nfactors = nfac, rotate = "varimax")
  
  scores = as.data.frame(fit$scores)
  colnames(scores) = paste0("PC", seq_len(ncol(scores)))
  scores$scientificName = rownames(df)
  scores
}

Compute_functional_structure=function(Cover_k,
                                      Species_k,
                                      dataset,
                                      ID_year_plot){
  print(unique(ID_year_plot))
  dataset=unique(dataset)
  
  trait_sp=read.table(paste0("./Data/Global_trait_dataset_",ifelse(dataset=="Underplot","Underplot","Hairphae"),".csv"),sep=";")
  
  trait_sp$scientificName=rownames(trait_sp)
  
  if (!any(Cover_k) | length(Cover_k)==0){
    return(tibble(
      wFEve=NA,
      FEve=NA,
      wFDis=NA,
      FDis=NA,
      wFDiv=NA,
      FDiv=NA,
      RaoQ=NA,
      CWM_PC1=NA,
      CWvar_PC1=NA,
      CWskew_PC1=NA,
      CWkurt_PC1=NA,
      CWM_PC2=NA,
      CWvar_PC2=NA,
      CWskew_PC2=NA,
      CWkurt_PC2=NA,
      CWM_PC3=NA,
      CWvar_PC3=NA,
      CWskew_PC3=NA,
      CWkurt_PC3=NA,
      FD_PC1=NA,
      FD_PC2=NA,
      FD_PC3=NA,
      Trait_representativity=NA))
    
  }else{
    if (any(Cover_k==0)){ #remove species with cover = 0
      which_sp=which(Cover_k==0)
      Cover_k=Cover_k[-which_sp]
      Species_k=Species_k[-which_sp]
    }
    
    #Functional diversity
    Trait_representativity=100
    if (any(Species_k %!in% trait_sp$scientificName)){ #if species have no traits, remove them
      which_sp=which(Species_k %!in% trait_sp$scientificName)
      Trait_representativity=100*(1-sum(Cover_k[which_sp])/sum(Cover_k))
      Cover_k=Cover_k[-which_sp]
      Species_k=Species_k[-which_sp]
    }
    
    trait_sp_k=dplyr::filter(trait_sp,scientificName %in% Species_k)
    
    Cover_k=as.matrix(t(Cover_k))
    colnames(Cover_k)=Species_k
    
    Trait_matrix=as.matrix(trait_sp_k[,colnames(trait_sp_k)[grep(pattern = "PC",colnames(trait_sp_k))]]) #removing cover and name
    rownames(Trait_matrix)=colnames(Cover_k)
    
    #Diversity indices
    Dist_matrix=compute_dist_matrix(Trait_matrix,metric = "euclidean")
    Dist_matrix=Dist_matrix/max(Dist_matrix,na.rm = T)
    
    FD_trait_weighted=FD::dbFD(Dist_matrix,Cover_k)
    for (k in 1:ncol(Trait_matrix)) {
      tryCatch({
        
        Dist_matrix=compute_dist_matrix(Trait_matrix[, paste0("PC", k)],metric = "euclidean")
        Dist_matrix=Dist_matrix/max(Dist_matrix,na.rm = T)
        
        assign(paste0("FD_trait_weighted_", k),FD::dbFD(Dist_matrix,Cover_k),envir = .GlobalEnv)
      }, error = function(e) {
        assign(paste0("FD_trait_weighted_", k), list(FEve=NA,
                                                     FDiv=NA,
                                                     FDis=NA,
                                                     RaoQ=NA),envir = .GlobalEnv)
      })
    }
    
    # FD_trait_nonweighted=FD::dbFD(compute_dist_matrix(Trait_matrix,metric = "euclidean"),Cover_k,w.abun = F)
    
    Dist_matrix=compute_dist_matrix(Trait_matrix,metric = "euclidean")
    Dist_matrix=Dist_matrix/max(Dist_matrix,na.rm = T)
    
    RaoQ_traits_all=RaoQ_traits(Cover_k,Dist_matrix)
    RaoQ_traits_all[RaoQ_traits_all==1]=NA
    
    for (k in 1:ncol(Trait_matrix)){
      Dist_matrix=compute_dist_matrix(Trait_matrix[, paste0("PC", k)],metric = "euclidean")
      Dist_matrix=Dist_matrix/max(Dist_matrix,na.rm = T)
      
      RaoQ_traits_x=RaoQ_traits(Cover_k,Dist_matrix)
      RaoQ_traits_x[RaoQ_traits_x==1]=NA
      assign(paste0("RaoQ_traits_x_",k),RaoQ_traits_x)
    }    
    
    d_CWT_FD=trait_sp_k%>%
      dplyr::mutate(., Cover=as.numeric(Cover_k))%>%
      dplyr::mutate(., Cover=Cover/sum(Cover))%>%
      reshape2::melt(., measure.vars = colnames(trait_sp_k)[grep(pattern = "PC",colnames(trait_sp_k))])%>%
      dplyr::group_by(., variable)%>%
      dplyr::summarize(., .groups = "keep",
                       CWM = sum(Cover*value,na.rm = T), # community weighted mean
                       CW_Var = sum(Cover*((value-sum(Cover*value,na.rm = T))**2),na.rm = T), # community weighted variance
                       CW_skew = sum(Cover*((value-sum(Cover*value,na.rm = T))**3)/(sum(Cover*((value-sum(Cover*value,na.rm = T))**2),na.rm = T)**(3/2)),na.rm = T), # community weighted skewness
                       CW_kurt = sum(Cover*((value-sum(Cover*value,na.rm = T))**4)/(sum(Cover*((value-sum(Cover*value,na.rm = T))**2),na.rm = T)**2),na.rm = T), # community weighted kurtosis
                       FD = sum(Cover * (abs(value-CWM)/sum(abs(value-CWM))))
                       
      )
    
    if (dataset %in% c("Underplot","")){
      return(tibble(
        Species_richness=length(Species_k),
        FEve_all=FD_trait_weighted$FEve,
        FDis_all=FD_trait_weighted$FDis,
        FDiv_all=FD_trait_weighted$FDiv,
        RaoQ_all=RaoQ_traits_all,
        FEve_PC1=FD_trait_weighted_1$FEve,
        FDis_PC1=FD_trait_weighted_1$FDis,
        FDiv_PC1=FD_trait_weighted_1$FDiv,
        RaoQ_PC1=RaoQ_traits_x_1,
        FEve_PC2=FD_trait_weighted_2$FEve,
        FDis_PC2=FD_trait_weighted_2$FDis,
        FDiv_PC2=FD_trait_weighted_2$FDiv,
        RaoQ_PC2=RaoQ_traits_x_2,
        RaoQ_all_FD=FD_trait_weighted$RaoQ,
        RaoQ_PC1_FD=FD_trait_weighted_1$RaoQ,
        RaoQ_PC2_FD=FD_trait_weighted_2$RaoQ,
        CWM_PC1=d_CWT_FD$CWM[which(d_CWT_FD$variable=="PC1")],
        CWvar_PC1=d_CWT_FD$CW_Var[which(d_CWT_FD$variable=="PC1")],
        CWskew_PC1=d_CWT_FD$CW_skew[which(d_CWT_FD$variable=="PC1")],
        CWkurt_PC1=d_CWT_FD$CW_kurt[which(d_CWT_FD$variable=="PC1")],
        FD_PC1=d_CWT_FD$FD[which(d_CWT_FD$variable=="PC1")],
        CWM_PC2=d_CWT_FD$CWM[which(d_CWT_FD$variable=="PC2")],
        CWvar_PC2=d_CWT_FD$CW_Var[which(d_CWT_FD$variable=="PC2")],
        CWskew_PC2=d_CWT_FD$CW_skew[which(d_CWT_FD$variable=="PC2")],
        CWkurt_PC2=d_CWT_FD$CW_kurt[which(d_CWT_FD$variable=="PC2")],
        FD_PC2=d_CWT_FD$FD[which(d_CWT_FD$variable=="PC2")],
        Trait_representativity=Trait_representativity))
      
    }else{
      return(tibble(
        Species_richness=length(Species_k),
        FEve_all=FD_trait_weighted$FEve,
        FDis_all=FD_trait_weighted$FDis,
        FDiv_all=FD_trait_weighted$FDiv,
        RaoQ_all=RaoQ_traits_all,
        FEve_PC1=FD_trait_weighted_1$FEve,
        FDis_PC1=FD_trait_weighted_1$FDis,
        FDiv_PC1=FD_trait_weighted_1$FDiv,
        RaoQ_PC1=RaoQ_traits_x_1,
        FEve_PC2=FD_trait_weighted_2$FEve,
        FDis_PC2=FD_trait_weighted_2$FDis,
        FDiv_PC2=FD_trait_weighted_2$FDiv,
        RaoQ_PC2=RaoQ_traits_x_2,
        FEve_PC3=FD_trait_weighted_3$FEve,
        FDis_PC3=FD_trait_weighted_3$FDis,
        FDiv_PC3=FD_trait_weighted_3$FDiv,
        RaoQ_PC3=RaoQ_traits_x_3,
        RaoQ_all_FD=FD_trait_weighted$RaoQ,
        RaoQ_PC1_FD=FD_trait_weighted_1$RaoQ,
        RaoQ_PC2_FD=FD_trait_weighted_2$RaoQ,
        RaoQ_PC3_FD=FD_trait_weighted_3$RaoQ,
        CWM_PC1=d_CWT_FD$CWM[which(d_CWT_FD$variable=="PC1")],
        CWvar_PC1=d_CWT_FD$CW_Var[which(d_CWT_FD$variable=="PC1")],
        CWskew_PC1=d_CWT_FD$CW_skew[which(d_CWT_FD$variable=="PC1")],
        CWkurt_PC1=d_CWT_FD$CW_kurt[which(d_CWT_FD$variable=="PC1")],
        FD_PC1=d_CWT_FD$FD[which(d_CWT_FD$variable=="PC1")],
        CWM_PC2=d_CWT_FD$CWM[which(d_CWT_FD$variable=="PC2")],
        CWvar_PC2=d_CWT_FD$CW_Var[which(d_CWT_FD$variable=="PC2")],
        CWskew_PC2=d_CWT_FD$CW_skew[which(d_CWT_FD$variable=="PC2")],
        CWkurt_PC2=d_CWT_FD$CW_kurt[which(d_CWT_FD$variable=="PC2")],
        FD_PC2=d_CWT_FD$FD[which(d_CWT_FD$variable=="PC2")],
        CWM_PC3=d_CWT_FD$CWM[which(d_CWT_FD$variable=="PC3")],
        CWvar_PC3=d_CWT_FD$CW_Var[which(d_CWT_FD$variable=="PC3")],
        CWskew_PC3=d_CWT_FD$CW_skew[which(d_CWT_FD$variable=="PC3")],
        CWkurt_PC3=d_CWT_FD$CW_kurt[which(d_CWT_FD$variable=="PC3")],
        FD_PC3=d_CWT_FD$FD[which(d_CWT_FD$variable=="PC3")],
        Trait_representativity=Trait_representativity))
    }
  }
}



## MOVING WINDOWS -----
Plot_MV_figure = function(d, stability_var = "Resistance_Isbell",
                          CI_inner = 90,
                          CI_outer = 95,
                          grepl_character = "CWM",
                          add_signif = TRUE, phylo = FALSE, alpha_ = .08,
                          binary = FALSE, gradientname = NULL,
                          same_panel = FALSE) {
  
  d_plot = Rename_PCA_axes(d) %>%
    dplyr::filter(., Stability_var == stability_var) %>%
    dplyr::select(-Stability_var) %>%
    reshape2::melt(., id.vars = "Variable_gradient") %>%
    dplyr::rename(var = variable, value = value) %>%
    dplyr::filter(grepl(grepl_character, var)) %>%
    dplyr::filter(!is.na(value))
  
  if (binary) {
    
    d_plot = d_plot %>%
      dplyr::mutate(., Variable_gradient = recode_factor(Variable_gradient,
                                                         "NO"  = paste0("No ",   gradientname),
                                                         "YES" = paste0("With ", gradientname))) %>%
      dplyr::mutate(Variable_gradient = factor(Variable_gradient,
                                               levels = c(paste0("No ",   gradientname),
                                                          paste0("With ", gradientname))))
    
    d_sum = d_plot %>%
      dplyr::group_by(var, Variable_gradient) %>%
      dplyr::summarise(
        Stat    = median(value, na.rm = TRUE),
        q1_in   = quantile(value, (1 - CI_inner/100)/2, na.rm = TRUE),
        q3_in   = quantile(value, (1 - (1 - CI_inner/100)/2), na.rm = TRUE),
        q1_out  = quantile(value, (1 - CI_outer/100)/2, na.rm = TRUE),
        q3_out  = quantile(value, (1 - (1 - CI_outer/100)/2), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(significant = (q1_out > 0) | (q3_out < 0))
    
    d_sum = d_sum %>%
      dplyr::group_by(var) %>%
      dplyr::mutate(y_mark = max(q3_out, na.rm = TRUE) * 1.08) %>%
      dplyr::ungroup()
    
    p = ggplot(d_sum, aes(x = Variable_gradient, y = Stat, color = var)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_linerange(aes(ymin = q1_out, ymax = q3_out),
                     lwd = .5, position = position_jitterdodge(seed = 123)) +
      geom_linerange(aes(ymin = q1_in, ymax = q3_in),
                     lwd = 1.5, position = position_jitterdodge(seed = 123)) +
      geom_pointrange(aes(ymin = q1_in, ymax = q3_in),
                      color = "transparent",
                      size = .15, position = position_jitterdodge(seed = 123), shape = 21) +
      geom_point(color = "black", fill = "white",
                 size = 4, position = position_jitterdodge(seed = 123), shape = 21) +
      labs(x = "", y = "Effect size on drought resistance", title = "") +
      the_theme2 +
      theme(legend.position = "none") +
      scale_color_manual(values = color_trait_axes)
    
    if (!same_panel) {
      p = p + facet_wrap(~var, scales = "free_y", nrow = 1) +
        theme(strip.background = element_rect(fill = "black"),
              strip.text       = element_text(colour = "white"))
    } else {
      p = p + theme(legend.position = "right")
    }
    
    if (phylo) {
      p = p + scale_color_manual(values = rep("black", 3))
    }
    
    return(p)
    
  } else {
    
    p = ggplot(d_plot, aes(x = Variable_gradient, y = value, group = var)) +
      geom_smooth(aes(color = var), se = FALSE) +
      stat_summary(aes(fill = var),
                   fun.data = function(x) {
                     data.frame(ymin = quantile(x, (1-(CI_outer/100))/2),
                                ymax = quantile(x, ((CI_outer/100))+(1-(CI_outer/100))/2))
                   },
                   geom = "ribbon", alpha = alpha_) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 1) +
      labs(x = "Land-use intensity",
           y = "Effect size on drought resistance",
           title = "") +
      the_theme2 +
      scale_color_manual(values = color_trait_axes) +
      scale_fill_manual(values  = fill_trait_axes)
    
    if (!same_panel) {
      p = p + facet_wrap(~var, scales = "free_y", nrow = 1) +
        theme(strip.background = element_rect(fill = "black"),
              strip.text       = element_text(colour = "white"),
              legend.position  = "none")
    } else {
      p = p + theme(legend.position = "right")
    }
    
    if (phylo) {
      p = p + scale_color_manual(values = rep("black", 3)) +
        scale_fill_manual(values = rep("grey", 3))
    }
    
    if (add_signif) {
      d_sig = d_plot %>%
        dplyr::group_by(var, Variable_gradient) %>%
        dplyr::summarise(lwr = quantile(value, (1-(CI_outer/100))/2),
                         upr = quantile(value, ((CI_outer/100))+(1-(CI_outer/100))/2),
                         .groups = "drop") %>%
        dplyr::mutate(significant = (lwr > 0) | (upr < 0))
      
      d_sig = d_sig %>%
        dplyr::group_by(var) %>%
        dplyr::mutate(y_mark = max(upr, na.rm = TRUE) * 1.08) %>%
        dplyr::ungroup()
      
      p = p + geom_point(data = dplyr::filter(d_sig, significant),
                         aes(x = Variable_gradient, y = y_mark, group = var),
                         shape = 108, size = 3, color = "grey")
    }
    
    return(p)
  }
}

Plot_MV_paired_by_axis_join = function(d_Rs, d_Rl,
                                       grepl_character = "CWM",
                                       CI_inner = 90,
                                       CI_outer = 95,
                                       alpha_ = .15,
                                       same_y = TRUE,
                                       negative_x=F,
                                       color_Rs = "#3366cc",
                                       color_Rl = "#FF6699") {
  
  prep = function(df, resp_label, stab) {
    df %>%
      dplyr::filter(Stability_var == stab) %>%
      dplyr::select(-Stability_var) %>%
      tidyr::pivot_longer(-Variable_gradient,
                          names_to = "Predictor", values_to = "Effect_size") %>%
      dplyr::filter(grepl(grepl_character, Predictor),
                    !is.na(Effect_size)) %>%
      dplyr::mutate(
        Axis = dplyr::case_when(
          grepl("PC1", Predictor) ~ "PC1: Economics (slow - fast)",
          grepl("PC2", Predictor) ~ "PC2: Absorption (AMF - hairs)",
          grepl("PC3", Predictor) ~ "PC3: Exploration (AD - SRL)",
          TRUE ~ "Other"
        ),
        Response = resp_label
      )
  }
  
  d_Rs_sub = prep(d_Rs, "Drought resistance", "Resistance_Isbell")
  d_Rl_sub = prep(d_Rl, "Drought resilience", "Resilience_Isbell_abs")
  
  d_plot = dplyr::bind_rows(d_Rs_sub, d_Rl_sub) %>%
    dplyr::filter(Axis != "Other") %>%
    dplyr::mutate(
      Response = factor(Response, levels = c("Drought resistance", "Drought resilience")),
      Variable_gradient = as.numeric(Variable_gradient)
    )
  
  if (negative_x) d_plot=d_plot%>%dplyr::mutate(., Variable_gradient=-Variable_gradient)
  
  if (nrow(d_plot) == 0) stop("No rows after filtering; check Predictor names.")
  
  d_sum = d_plot %>%
    dplyr::group_by(Axis, Response, Variable_gradient) %>%
    dplyr::summarise(
      q1_in  = quantile(Effect_size, (1 - CI_inner/100)/2, na.rm = TRUE),
      q3_in  = quantile(Effect_size, (1 - (1 - CI_inner/100)/2), na.rm = TRUE),
      q1_out = quantile(Effect_size, (1 - CI_outer/100)/2, na.rm = TRUE),
      q3_out = quantile(Effect_size, (1 - (1 - CI_outer/100)/2), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(significant = (q1_out > 0) | (q3_out < 0))
  
  d_sig = d_sum %>%
    dplyr::filter(significant) %>%
    dplyr::group_by(Axis) %>%
    dplyr::mutate(
      y_top  = max(q3_out, na.rm = TRUE),
      y_mark = dplyr::case_when(
        Response == "Drought resistance" ~ y_top * 1.06,
        Response == "Drought resilience" ~ y_top * 1.18
      )
    ) %>%
    dplyr::ungroup()
  
  p = ggplot(d_plot, aes(x = Variable_gradient, y = Effect_size,
                         group = Response, color = Response, fill = Response)) +
    geom_smooth(se = FALSE, linewidth = 1) +
    stat_summary(
      fun.data = function(x) {
        data.frame(ymin = quantile(x, (1 - CI_outer/100)/2, na.rm = TRUE),
                   ymax = quantile(x, (1 - (1 - CI_outer/100)/2), na.rm = TRUE))
      },
      geom = "ribbon", alpha = alpha_, color = NA
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_linerange(data = d_sig,
                   aes(x = Variable_gradient,
                       ymin = y_mark - 0.02 * abs(y_mark),
                       ymax = y_mark + 0.02 * abs(y_mark),
                       color = Response),
                   linewidth = 1.5, show.legend = FALSE, inherit.aes = FALSE) +
    geom_linerange(data = d_sig,
                   aes(x = Variable_gradient,
                       ymin = y_mark - 0.008 * abs(y_mark),
                       ymax = y_mark + 0.008 * abs(y_mark),
                       color = Response),
                   linewidth = 0.5, show.legend = FALSE, inherit.aes = FALSE) +
    scale_color_manual(values = c("Drought resistance" = color_Rs,
                                  "Drought resilience" = color_Rl)) +
    scale_fill_manual(values  = c("Drought resistance" = color_Rs,
                                  "Drought resilience" = color_Rl)) +
    facet_wrap(~Axis, nrow = 1, scales = if (same_y) "fixed" else "free_y") +
    labs(x = "Land-use intensity", y = "Effect size",
         color = "", fill = "") +
    the_theme2 +
    theme(strip.background = element_rect(fill = "black"),
          strip.text       = element_text(colour = "white"),
          legend.position  = "bottom")
  
  return(p)
}


smooth_rolling_median = function(raw, group_col, window = 3, n_grid = 300) {
  
  if (window < 3) window <- 3
  if (window %% 2 == 0) window <- window + 1
  
  raw %>%
    dplyr::group_split(.data[[group_col]]) %>%
    purrr::map_dfr(function(df) {
      df = df %>% dplyr::arrange(Variable_gradient)
      if (nrow(df) < window) return(df)
      
      xg   = seq(min(df$Variable_gradient), max(df$Variable_gradient),
                 length.out = n_grid)
      y_sm = stats::runmed(df$fraction, k = window, endrule = "median")
      
      tibble::tibble(
        Variable_gradient = xg,
        Stability_var     = unique(df$Stability_var),
        !!group_col       := unique(df[[group_col]]),
        fraction          = approx(df$Variable_gradient, y_sm,
                                   xout = xg, rule = 2)$y
      )
    }) %>%
    dplyr::group_by(Variable_gradient, Stability_var) %>%
    dplyr::mutate(fraction = pmax(fraction, 0),
                  fraction = fraction / sum(fraction, na.rm = TRUE)) %>%
    dplyr::ungroup()
}

Plot_variance_partitioning = function(d, alpha1 = .85, alpha2 = .5,
                                      negative_x     = F,
                                      stability_var  = "Resistance_Isbell",
                                      binary         = FALSE,
                                      smooth         = TRUE,
                                      smooth_window  = 3,
                                      n_grid         = 300) {
  
  frac_var_explained = d %>%
    dplyr::rename(var = Predictor) %>%
    dplyr::mutate(
      PC_group = dplyr::case_when(
        grepl("PC1|PC2|PC3", var) ~ "PC_axes",
        TRUE ~ "Other"),
      var_type = dplyr::case_when(
        grepl("CWM", var) ~ "Identity dominant sp.",
        grepl("RaoQ|FDis|FEve|FDiv|FD", var) ~ "Functional div.",
        grepl("PSV|MNTD", var) ~ "Phylogenetic div.",
        grepl("Species_richness", var) ~ "Taxonomic div.",
        TRUE ~ "Geography & climate"),
      PC_axis = dplyr::case_when(
        grepl("PC1", var) ~ "PC1: Economics (slow - fast)",
        grepl("PC2", var) ~ "PC2: Absorption (AMF - hairs)",
        grepl("PC3", var) ~ "PC3: Exploration (AD - SRL)",
        TRUE ~ "Other"),
      PC_axis2 = dplyr::case_when(
        grepl("PC1", var) & !grepl("CWM_PC1", var) ~ "Diversity (Economics)",
        grepl("PC2", var) & !grepl("CWM_PC2", var) ~ "Diversity (Absorption)",
        grepl("PC3", var) & !grepl("CWM_PC3", var) ~ "Diversity (Exploration)",
        grepl("CWM_PC1", var) ~ "Identity (Economics)",
        grepl("CWM_PC2", var) ~ "Identity (Absorption)",
        grepl("CWM_PC3", var) ~ "Identity (Exploration)",
        grepl("Species_richness|PSV|MNTD", var) ~ "Phylogenetic & taxo div.",
        TRUE ~ "Other"),
      PC_axis3 = dplyr::case_when(
        grepl("PC1", var) ~ "PC1: Economics (slow - fast)",
        grepl("PC2", var) ~ "PC2: Absorption (AMF - hairs)",
        grepl("PC3", var) ~ "PC3: Exploration (AD - SRL)",
        grepl("Species_richness|PSV|MNTD", var) ~ "Phylogenetic & taxo div.",
        TRUE ~ "Other")
    )
  
  base_data = frac_var_explained %>%
    dplyr::mutate(Effect_size = Effect_size^2) %>%
    dplyr::filter(Stability_var == stability_var, !is.na(var_type))
  
  pastel_vartype = c(
    "Identity dominant sp." = "#AEC6CF",
    "Functional div."       = "#FFB7B2",
    "Phylogenetic div."     = "#B5EAD7",
    "Taxonomic div."        = "#FFDAC1",
    "Geography & climate"   = "grey20"
  )
  pastel_pcaxis = c(
    "PC1: Economics (slow - fast)"  = "#FFDAC1",
    "PC2: Absorption (AMF - hairs)" = "#FF6699",
    "PC3: Exploration (AD - SRL)"   = "#3366cc"
  )
  pastel_pcaxis3 = c(
    "PC1: Economics (slow - fast)"  = "#C7CEEA",
    "PC2: Absorption (AMF - hairs)" = "#FFDAC1",
    "PC3: Exploration (AD - SRL)"   = "#B5EAD7",
    "Phylogenetic & taxo div."      = "grey"
  )
  pastel_pcaxis2 = c(
    "Diversity (Exploration)"  = "#82B7A4",
    "Diversity (Economics)"    = "#8F9DD8",
    "Diversity (Absorption)"   = "#E2AA84",
    "Identity (Economics)"     = "#C7CEEA",
    "Identity (Absorption)"    = "#FFDAC1",
    "Phylogenetic & taxo div." = "grey",
    "Identity (Exploration)"   = "#B5EAD7"
  )
  
  if (negative_x) {
    base_data = base_data %>% dplyr::mutate(Variable_gradient = -Variable_gradient)
  }
  
  if (binary) {
    base_data = base_data %>%
      dplyr::mutate(Variable_gradient = factor(Variable_gradient,
                                               levels = c("NO", "YES")))
    
    frac_var_type = base_data %>%
      dplyr::group_by(Variable_gradient, Stability_var, var_type) %>%
      dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
      dplyr::group_by(Variable_gradient, Stability_var) %>%
      dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    p1 = ggplot(frac_var_type,
                aes(x = Variable_gradient, y = fraction, fill = var_type)) +
      geom_col(alpha = alpha1, width = 0.6) +
      scale_fill_manual(values = pastel_vartype) +
      labs(x = "Fertilization treatment",
           y = "Fraction of explained variance", fill = "") +
      the_theme2 + guides(fill = guide_legend(nrow = 2))
    
    frac_var_pc = base_data %>%
      dplyr::filter(PC_group == "PC_axes") %>%
      dplyr::group_by(Variable_gradient, Stability_var, PC_axis) %>%
      dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
      dplyr::group_by(Variable_gradient, Stability_var) %>%
      dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    p2 = ggplot(frac_var_pc,
                aes(x = Variable_gradient, y = fraction, fill = PC_axis)) +
      geom_col(alpha = alpha1, width = 0.6) +
      scale_fill_manual(values = pastel_pcaxis) +
      labs(x = "Fertilization treatment",
           y = "Fraction of explained variance", fill = "") +
      the_theme2 + guides(fill = guide_legend(nrow = 2))
    
    frac_var_pc2 = base_data %>%
      dplyr::filter(PC_axis2 != "Other") %>%
      dplyr::group_by(Variable_gradient, Stability_var, PC_axis2) %>%
      dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
      dplyr::group_by(Variable_gradient, Stability_var) %>%
      dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    p3 = ggplot(frac_var_pc2,
                aes(x = Variable_gradient, y = fraction, fill = PC_axis2)) +
      geom_col(alpha = alpha2, width = 0.6) +
      scale_fill_manual(values = pastel_pcaxis2) +
      labs(x = "Fertilization treatment",
           y = "Fraction of explained variance", fill = "") +
      the_theme2 + guides(fill = guide_legend(nrow = 3))
    
    frac_var_pc3 = base_data %>%
      dplyr::filter(PC_axis3 != "Other") %>%
      dplyr::group_by(Variable_gradient, Stability_var, PC_axis3) %>%
      dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
      dplyr::group_by(Variable_gradient, Stability_var) %>%
      dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
      dplyr::ungroup()
    
    p4 = ggplot(frac_var_pc3,
                aes(x = Variable_gradient, y = fraction, fill = PC_axis3)) +
      geom_col(alpha = alpha1, width = 0.6) +
      scale_fill_manual(values = pastel_pcaxis3) +
      labs(x = "Fertilization treatment",
           y = "Fraction of explained variance", fill = "") +
      the_theme2 + guides(fill = guide_legend(nrow = 2))
    
    return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
  }
  
  frac_var_explained_vartype = base_data %>%
    dplyr::group_by(Variable_gradient, Stability_var, var_type) %>%
    dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
    dplyr::group_by(Variable_gradient, Stability_var) %>%
    dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  if (smooth) {
    frac_var_explained_vartype = smooth_rolling_median(
      frac_var_explained_vartype, group_col = "var_type",
      window = smooth_window, n_grid = n_grid)
  }
  
  p1 = ggplot(frac_var_explained_vartype,
              aes(x = Variable_gradient, y = fraction, fill = var_type)) +
    geom_area(alpha = alpha1, position = "stack") +
    scale_fill_manual(values = pastel_vartype) +
    labs(x = "LUI", y = "Fraction of explained variance", fill = "") +
    the_theme2 + guides(fill = guide_legend(nrow = 2))
  
  frac_pcaxis = base_data %>%
    dplyr::filter(PC_group == "PC_axes") %>%
    dplyr::group_by(Variable_gradient, Stability_var, PC_axis) %>%
    dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
    dplyr::group_by(Variable_gradient, Stability_var) %>%
    dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  if (smooth) {
    frac_pcaxis = smooth_rolling_median(
      frac_pcaxis, group_col = "PC_axis",
      window = smooth_window, n_grid = n_grid)
  }
  
  p2 = ggplot(frac_pcaxis,
              aes(x = Variable_gradient, y = fraction, fill = PC_axis)) +
    geom_area(alpha = alpha1, position = "stack") +
    scale_fill_manual(values = pastel_pcaxis) +
    labs(x = "LUI", y = "Fraction of explained variance", fill = "") +
    the_theme2 + guides(fill = guide_legend(nrow = 2))
  
  frac_pcaxis2 = base_data %>%
    dplyr::filter(PC_axis2 != "Other") %>%
    dplyr::group_by(Variable_gradient, Stability_var, PC_axis2) %>%
    dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
    dplyr::group_by(Variable_gradient, Stability_var) %>%
    dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  if (smooth) {
    frac_pcaxis2 = smooth_rolling_median(
      frac_pcaxis2, group_col = "PC_axis2",
      window = smooth_window, n_grid = n_grid)
  }
  
  p3 = ggplot(frac_pcaxis2,
              aes(x = Variable_gradient, y = fraction, fill = PC_axis2)) +
    geom_area(alpha = alpha2, position = "stack") +
    scale_fill_manual(values = pastel_pcaxis2) +
    labs(x = "LUI", y = "Fraction of explained variance", fill = "") +
    the_theme2 + guides(fill = guide_legend(nrow = 3))
  
  frac_pcaxis3 = base_data %>%
    dplyr::filter(PC_axis3 != "Other") %>%
    dplyr::group_by(Variable_gradient, Stability_var, PC_axis3) %>%
    dplyr::summarise(total_var = sum(Effect_size, na.rm = TRUE), .groups = "drop_last") %>%
    dplyr::group_by(Variable_gradient, Stability_var) %>%
    dplyr::mutate(fraction = total_var / sum(total_var, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  if (smooth) {
    frac_pcaxis3 = smooth_rolling_median(
      frac_pcaxis3, group_col = "PC_axis3",
      window = smooth_window, n_grid = n_grid)
  }
  
  p4 = ggplot(frac_pcaxis3,
              aes(x = Variable_gradient, y = fraction, fill = PC_axis3)) +
    geom_area(alpha = alpha1, position = "stack") +
    scale_fill_manual(values = pastel_pcaxis3) +
    labs(x = "LUI", y = "Fraction of explained variance", fill = "") +
    the_theme2 + guides(fill = guide_legend(nrow = 2))
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}

Indiv_trait="AD"
dataset = "HairPhae"
stability_var="Resistance_Isbell"
n_min_sites = NULL
n_site_increase = 10
variable_gradient = "LUI_global"
N_bootstrap = 100
SPEI = "spei06"
climate_params = c("longitude","latitude","SPEI_value")
random = "both"
phylo=T
drought_category = c("Moderate drought","Extreme drought")
diversity_type = c("FD")
correct_clim_residuals=T

Moving_window_gradient_diversity = function(dataset = "Underplot",
                                            stability_var,
                                            n_min_sites = NULL,
                                            n_site_increase = 10,
                                            variable_gradient = "LUI_global",
                                            N_bootstrap = 100,
                                            shuffle_gradient=F,
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
  } else if (diversity_type=="FD") {  
    diversity_predictors = c("Species_richness",
                             "FD_PC1", "FD_PC2", "FD_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
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
  print(n_min_sites)
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
                  all_of(diversity_predictors)) 
  
  if (shuffle_gradient) d_mod$Variable_gradient=sample(d_mod$Variable_gradient,size = nrow(d_mod),replace = F)
  
  
  d_mod=d_mod%>%
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
      
      # Store results
      all_variance = rbind(all_variance, 
                           data.frame(Effect_size = colMeans(boot_mod$t),
                                      Predictor = colnames(boot_mod$t)) %>%
                             add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
      
      
      if (isFALSE(shuffle_gradient)){
        
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



Moving_window_gradient_diversity_binary = function(dataset = "Underplot",
                                                   stability_var,
                                                   n_min_sites = NULL,
                                                   n_site_increase = 10,
                                                   variable_gradient = "LUI_global",
                                                   variable_gradient_type = c("continuous", "binned"),
                                                   N_bootstrap = 100,
                                                   SPEI = "spei06",
                                                   climate_params = c(""),
                                                   correct_clim_residuals = TRUE,
                                                   random = "both",
                                                   phylo = TRUE,
                                                   drought_category = c("Moderate drought", "Extreme drought"),
                                                   diversity_type = c("FD")) {
  
  
  
  # ---- diversity predictors  
  if (diversity_type == "RaoQ") {
    diversity_predictors = c("Species_richness",
                             "RaoQ_PC1", "RaoQ_PC2", "RaoQ_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type == "FD") {
    diversity_predictors = c("Species_richness",
                             "FD_PC1", "FD_PC2", "FD_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type == "FDis") {
    diversity_predictors = c("Species_richness",
                             "FDis_PC1", "FDis_PC2", "FDis_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type == "all") {
    diversity_predictors = c("Species_richness",
                             "FDis_PC1", "FDis_PC2", "FDis_PC3",
                             "FEve_PC1", "FEve_PC2", "FEve_PC3",
                             "FDiv_PC1", "FDiv_PC2", "FDiv_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  }
  
  if (phylo) {
    diversity_predictors = c(diversity_predictors, "PSV", "MNTD_weighted")
  }
  
  if (is.null(n_min_sites)) {
    n_min_sites = (length(diversity_predictors) +
                     ifelse(length(climate_params) == 1 | correct_clim_residuals,
                            0, length(climate_params)) + 3) * 10
  }
  
  data_path = paste0("./Data/Merged_datasets_", dataset, ".rds")
  
  if (stability_var == "Temporal_stability") {
    d_file = readRDS(data_path)[[SPEI]] %>%
      dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
      dplyr::filter(., Trait_representativity > 70) %>%
      dplyr::mutate(., Species_richness = sqrt(Species_richness)) %>%
      dplyr::distinct(., .keep_all = TRUE, Site)
    
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
                    Recovery = bcPower(Recovery, -.2),
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
    
  } else if (stability_var == "Resilience_wet") {
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
  
  d_mod = d_file %>%
    dplyr::mutate(.,
                  Variable_gradient = as.numeric(dplyr::pull(., dplyr::all_of(variable_gradient))),
                  Response_var = as.numeric(dplyr::pull(., dplyr::all_of(stability_var))))
  
  if (variable_gradient_type == "binned") {
    d_mod$Variable_gradient = d_mod$Variable_gradient == 0   # TRUE = unfertilized, ungrazed, unmown
  }
  
  d_mod = d_mod %>%
    dplyr::select(.,
                  Response_var,
                  Variable_gradient,
                  Site, Region, any_of(climate_params),
                  Year,
                  all_of(diversity_predictors)) %>%
    dplyr::arrange(., Variable_gradient) %>%
    dplyr::mutate(., Variable_gradient = as.character(Variable_gradient)) %>%
    drop_na(.)
  
  if (climate_params[1] != "" & correct_clim_residuals) {
    d_mod$Response_var = residuals(lm(as.formula(paste0("Response_var ~ ",
                                                        paste0(climate_params, collapse = "+"))),
                                      data = d_mod))
    d_mod = d_mod %>%
      dplyr::select(., -any_of(climate_params))
  }
  
  all_results = all_variance = all_R2 =
    all_partial_residuals = all_spatial_autocorrelation =
    all_temporal_autocorrelation = tibble()
  
  print(dim(d_mod))
  
  if (variable_gradient_type == "binned") {
    
    group_levels = c("TRUE", "FALSE")
    group_labels = c("NO", "YES")   
    
    for (g in seq_along(group_levels)) {
      
      save_d = d_mod[d_mod$Variable_gradient == group_levels[g], ]%>%
        dplyr::select(.,-Variable_gradient)
      
      # Scale
      data_model = save_d %>%
        Numeric_scaling_df(.)
      
      # Fit LMM 
      if (random == "both") {
        mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Site) + (1|Year) + ",
                                             paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                           data = data_model, REML = TRUE,
                           control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
      } else if (random %in% c("Site", "site")) {
        mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Site) + ",
                                             paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                           data = data_model, REML = TRUE,
                           control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
      } else if (random %in% c("Year", "year")) {
        mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Year) + ",
                                             paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                           data = data_model, REML = TRUE,
                           control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
      } else if (random %in% c("Region", "region")) {
        mod_drought = lmer(as.formula(paste0("Response_var ~ (1|Region) + ",
                                             paste0(colnames(data_model %>% dplyr::select(., -Response_var, -Site, -Year, -Region)), collapse = "+"))),
                           data = data_model, REML = TRUE,
                           control = lmerControl(boundary.tol = 1e-6), na.action = na.fail)
      } else if (random == "lme") {
        mod_drought = gls(as.formula(paste0("Response_var ~ ",
                                            paste0(colnames(data_model %>%
                                                              dplyr::select(., -Response_var, -Site, -Year, -Region)),
                                                   collapse = "+"))),
                          data = data_model, method = "ML", na.action = na.fail,
                          control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
        mod_drought = update(mod_drought, correlation = corCAR1(form = as.formula(paste0("~ Year|Site"))))
      } else if (random == "GAM") {
        predictors = colnames(data_model %>% dplyr::select(-Response_var, -Site, -Region))
        spatial_predictors = c("longitude", "latitude")
        linear_predictors = predictors[!predictors %in% spatial_predictors]
        
        mod_drought = gam(
          as.formula(paste0("Response_var ~ ",
                            paste0(linear_predictors, collapse = " + "),
                            " + s(longitude, latitude)")),
          data = data_model, method = "REML", na.action = na.omit)
        mod_drought = gamm(
          as.formula(paste0("Response_var ~ ",
                            paste0(linear_predictors, collapse = " + "),
                            " + s(longitude, latitude) + s(Site, bs='re')")),
          data = data_model, method = "REML", na.action = na.omit)
      }
      
      if (random != "lme") {
        # R2
        all_R2 = rbind(all_R2,
                       tibble(R2  = r.squaredGLMM(mod_drought)[1],
                              R2c = r.squaredGLMM(mod_drought)[2]) %>%
                         add_column(., Variable_gradient = group_labels[g]))
        
        # Bootstrap coefficients
        if (random != "GAM") {
          boot_strapmod = fixef(arm::sim(mod_drought))[, -1]
          boot_mod = list()
          boot_mod$t = boot_strapmod
        } else {
          boot_mod = list()
          boot_mod$t = tibble()
        }
        
        # Partial residuals
        df = tibble()
        for (pred in diversity_predictors) {
          partial_res = visreg::visreg(mod_drought, pred, plot = FALSE)
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
        
        all_partial_residuals = rbind(
          all_partial_residuals,
          df %>% add_column(., Variable_gradient = group_labels[g])
        )
        
        all_variance = rbind(
          all_variance,
          data.frame(Effect_size = colMeans(boot_mod$t),
                     Predictor   = colnames(boot_mod$t)) %>%
            add_column(., Variable_gradient = group_labels[g])
        )
        
        all_results = rbind(
          all_results,
          data.frame(boot_mod$t) %>%
            add_column(., Variable_gradient = group_labels[g])
        )
        
        # Spatial autocorrelation
        data_model$location_id = paste(d_file$longitude[match(data_model$Site, d_file$Site)],
                                       d_file$latitude[match(data_model$Site, d_file$Site)],
                                       sep = "_")
        sim_res_agg = recalculateResiduals(simulateResiduals(mod_drought),
                                           group = data_model$location_id)
        unique_coords = data_model %>%
          dplyr::mutate(.,
                        longitude = d_file$longitude[match(data_model$Site, d_file$Site)],
                        latitude  = d_file$latitude[match(data_model$Site, d_file$Site)]) %>%
          dplyr::distinct(longitude, latitude) %>%
          dplyr::arrange(longitude, latitude)
        spatial_test = testSpatialAutocorrelation(sim_res_agg,
                                                  x = unique_coords$longitude,
                                                  y = unique_coords$latitude, plot = FALSE)
        all_spatial_autocorrelation = rbind(
          all_spatial_autocorrelation,
          tibble(Statistic = spatial_test$statistic[1],
                 Pvalue    = spatial_test$p.value) %>%
            add_column(., Variable_gradient = group_labels[g])
        )
        
        # Temporal autocorrelation
        sim_res_agg_time = recalculateResiduals(simulateResiduals(mod_drought),
                                                group = data_model$Year)
        unique_years = data_model %>%
          dplyr::distinct(Year) %>%
          dplyr::arrange(Year)
        temporal_test = testTemporalAutocorrelation(sim_res_agg_time,
                                                    time = unique_years$Year, plot = FALSE)
        all_temporal_autocorrelation = rbind(
          all_temporal_autocorrelation,
          tibble(Statistic = temporal_test$statistic[1],
                 Pvalue    = temporal_test$p.value) %>%
            add_column(., Variable_gradient = group_labels[g])
        )
      }
    }
    
    return(list(
      Effect_size             = all_results             %>% add_column(., Stability_var = stability_var),
      Coeff_variance          = all_variance            %>% add_column(., Stability_var = stability_var),
      R2_model                = all_R2                  %>% add_column(., Stability_var = stability_var),
      Partial_residuals       = all_partial_residuals   %>% add_column(., Stability_var = stability_var),
      Spatial_autorrelation   = all_spatial_autocorrelation  %>% add_column(., Stability_var = stability_var),
      Temporal_autocorrelation= all_temporal_autocorrelation %>% add_column(., Stability_var = stability_var)
    ))
  }
  return(list(
    Effect_size             = all_results             %>% add_column(., Stability_var = stability_var),
    Coeff_variance          = all_variance            %>% add_column(., Stability_var = stability_var),
    R2_model                = all_R2                  %>% add_column(., Stability_var = stability_var),
    Partial_residuals       = all_partial_residuals   %>% add_column(., Stability_var = stability_var),
    Spatial_autorrelation   = all_spatial_autocorrelation  %>% add_column(., Stability_var = stability_var),
    Temporal_autocorrelation= all_temporal_autocorrelation %>% add_column(., Stability_var = stability_var)
  ))
}




Moving_window_gradient_diversity_paired = function(dataset = "HairPhae",
                                                   stability_var = "Resistance_Isbell",
                                                   n_min_sites = NULL,
                                                   n_site_increase = 2,
                                                   variable_gradient = "LUI_global",
                                                   N_bootstrap = 100,
                                                   shuffle_gradient = FALSE,
                                                   SPEI = "spei06",
                                                   climate_params = c(""),
                                                   correct_clim_residuals = TRUE,
                                                   random = "both",
                                                   phylo = TRUE,
                                                   shuffle_traits=F,
                                                   drought_category = c("Moderate drought", "Extreme drought"),
                                                   diversity_type = c("FD")) {
  
  
  
  
  
  if (diversity_type == "RaoQ") {
    diversity_predictors = c("Species_richness",
                             "RaoQ_PC1", "RaoQ_PC2", "RaoQ_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type == "FD") {
    diversity_predictors = c("Species_richness",
                             "FD_PC1", "FD_PC2", "FD_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type == "FDis") {
    diversity_predictors = c("Species_richness",
                             "FDis_PC1", "FDis_PC2", "FDis_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  } else if (diversity_type == "all") {
    diversity_predictors = c("Species_richness",
                             "FDis_PC1", "FDis_PC2", "FDis_PC3",
                             "FEve_PC1", "FEve_PC2", "FEve_PC3",
                             "FDiv_PC1", "FDiv_PC2", "FDiv_PC3",
                             "CWM_PC1", "CWM_PC2", "CWM_PC3")
  }
  
  if (phylo) {
    diversity_predictors = c(diversity_predictors, "PSV", "MNTD_weighted")
  }
  
  if (is.null(n_min_sites)) {
    n_min_sites = (length(diversity_predictors) +
                     ifelse(length(climate_params) == 1 | correct_clim_residuals,
                            0, length(climate_params)) + 3) * 10
  }
  
  # which response is the focal one
  if (stability_var == "Resistance_Isbell") {
    focal_col  = "Response_Rs"
    other_col  = "Response_Rl"
    response_label = "Rs"
  } else if (stability_var == "Resilience_Isbell_abs") {
    focal_col  = "Response_Rl"
    other_col  = "Response_Rs"
    response_label = "Rl"
  } else {
    stop("Only 'Resistance_Isbell' and 'Resilience_Isbell_abs' are supported.")
  }
  
  data_path = paste0("./Data/Merged_datasets_", dataset, ".rds")
  
  # common event table: keep events where BOTH responses are computable
  d_file = readRDS(data_path)[[SPEI]] %>%
    dplyr::mutate(., Region = sub("\\d+$", "", .$Site)) %>%
    dplyr::filter(.,
                  Trait_representativity > 70,
                  Drought_category %in% drought_category,
                  Resistance_Isbell > 0,
                  !is.na(Resistance_Isbell),
                  !is.na(Resilience_Isbell_abs),
                  is.na(Recovery_context)) %>%
    dplyr::mutate(.,
                  Resistance_Isbell     = bcPower(Resistance_Isbell, -0.7),
                  Resilience_Isbell_abs = log(Resilience_Isbell_abs),
                  Species_richness      = sqrt(Species_richness))
  
  d_mod = d_file %>%
    dplyr::mutate(.,
                  Variable_gradient = as.numeric(dplyr::pull(., dplyr::all_of(variable_gradient))),
                  Response_Rs       = as.numeric(Resistance_Isbell),
                  Response_Rl       = as.numeric(Resilience_Isbell_abs)) %>%
    dplyr::select(.,
                  Response_Rs, Response_Rl,
                  Variable_gradient,
                  Site, Region, any_of(climate_params),
                  Year,
                  all_of(diversity_predictors))
  
  if (shuffle_traits){
    trait_cols = grep("CWM_|FD_|Species_richness|PSV|MNTD", colnames(d_mod), value = TRUE)
    d_mod[trait_cols] = d_mod[sample.int(nrow(d_mod)), trait_cols]
  }
  
  
  if (shuffle_gradient) {
    d_mod$Variable_gradient = sample(d_mod$Variable_gradient,
                                     size = nrow(d_mod), replace = FALSE)
  }
  
  
  d_mod=d_mod%>%
    dplyr::arrange(., Variable_gradient) %>%
    dplyr::mutate(., Variable_gradient = as.character(Variable_gradient)) %>%
    drop_na(.)
  
  
  if (climate_params[1] != "" & correct_clim_residuals) {
    d_mod$Response_Rs = residuals(lm(
      as.formula(paste0("Response_Rs ~ ", paste0(climate_params, collapse = "+"))),
      data = d_mod))
    d_mod$Response_Rl = residuals(lm(
      as.formula(paste0("Response_Rl ~ ", paste0(climate_params, collapse = "+"))),
      data = d_mod))
    d_mod = d_mod %>% dplyr::select(., -any_of(climate_params))
  }
  
  all_results = all_variance = all_R2 =
    all_partial_residuals = all_spatial_autocorrelation =
    all_temporal_autocorrelation = tibble()
  
  print(dim(d_mod))
  
  for (index in 1:326) {
    
    if (index %% 20 == 0) print(index)
    
    id_max = ifelse((n_min_sites + n_site_increase * (index - 1)) > nrow(d_mod),
                    nrow(d_mod),
                    (n_min_sites + n_site_increase * (index - 1)))
    id_min = ifelse((id_max - n_min_sites) < 1, 1, id_max - n_min_sites)
    save_d = d_mod[id_min:id_max, ]
    
    data_model = d_mod[id_min:id_max, ] %>%
      dplyr::mutate(., Variable_gradient = as.numeric(Variable_gradient)) %>%
      Numeric_scaling_df(.)
    
    predictors_formula = paste0(
      setdiff(colnames(data_model),
              c("Response_Rs", "Response_Rl",
                "Site", "Year", "Region")),
      collapse = "+")
    
    fit_one = function(resp) {
      f  = paste0("Response ~ (1|Site) + (1|Year) + ", predictors_formula)
      dd = data_model
      dd$Response = dd[[resp]]
      tryCatch(
        lmer(as.formula(f), data = dd, REML = TRUE,
             control = lmerControl(boundary.tol = 1e-6), na.action = na.fail),
        error = function(e) NULL)
    }
    
    m_focal = fit_one(focal_col)
    m_other = fit_one(other_col)
    
    if (is.null(m_focal)) {
      if (id_max == nrow(d_mod)) break
      next
    }
    
    boot_focal = fixef(arm::sim(m_focal))[, -1]
    boot_mod = list()
    boot_mod$t = boot_focal
    
    all_R2 = rbind(all_R2,
                   tibble(R2  = r.squaredGLMM(m_focal)[1],
                          R2c = r.squaredGLMM(m_focal)[2]) %>%
                     add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
    
    all_variance = rbind(all_variance,
                         data.frame(Effect_size = colMeans(boot_mod$t),
                                    Predictor   = colnames(boot_mod$t)) %>%
                           add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient))))
    if (isFALSE(shuffle_gradient)) {
      df = tibble()
      for (pred in diversity_predictors) {
        partial_res = visreg::visreg(m_focal, pred, plot = FALSE)
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
      
      all_partial_residuals = rbind(
        all_partial_residuals,
        df %>% add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient)))
      )
      
      all_results = rbind(
        all_results,
        data.frame(boot_mod$t) %>%
          add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient)))
      )
      
      data_model$location_id = paste(d_file$longitude[match(data_model$Site, d_file$Site)],
                                     d_file$latitude[match(data_model$Site, d_file$Site)],
                                     sep = "_")
      
      sim_res_agg = recalculateResiduals(simulateResiduals(m_focal),
                                         group = data_model$location_id)
      
      unique_coords = data_model %>%
        dplyr::mutate(.,
                      longitude = d_file$longitude[match(data_model$Site, d_file$Site)],
                      latitude  = d_file$latitude[match(data_model$Site, d_file$Site)]) %>%
        dplyr::distinct(longitude, latitude) %>%
        dplyr::arrange(longitude, latitude)
      
      spatial_test = testSpatialAutocorrelation(sim_res_agg,
                                                x = unique_coords$longitude,
                                                y = unique_coords$latitude, plot = FALSE)
      
      all_spatial_autocorrelation = rbind(
        all_spatial_autocorrelation,
        tibble(Statistic = spatial_test$statistic[1],
               Pvalue    = spatial_test$p.value) %>%
          add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient)))
      )
      
      sim_res_agg_time = recalculateResiduals(simulateResiduals(m_focal),
                                              group = data_model$Year)
      
      unique_years = data_model %>%
        dplyr::distinct(Year) %>%
        dplyr::arrange(Year)
      
      temporal_test = testTemporalAutocorrelation(sim_res_agg_time,
                                                  time = unique_years$Year, plot = FALSE)
      
      all_temporal_autocorrelation = rbind(
        all_temporal_autocorrelation,
        tibble(Statistic = temporal_test$statistic[1],
               Pvalue    = temporal_test$p.value) %>%
          add_column(., Variable_gradient = mean(as.numeric(save_d$Variable_gradient)))
      )
    }    
    if (id_max == nrow(d_mod)) break
  }
  
  return(list(
    Effect_size              = all_results %>%
      add_column(., Stability_var = stability_var),
    Coeff_variance           = all_variance %>%
      add_column(., Stability_var = stability_var),
    R2_model                 = all_R2 %>%
      add_column(., Stability_var = stability_var),
    Partial_residuals        = all_partial_residuals %>%
      add_column(., Stability_var = stability_var),
    Spatial_autorrelation    = all_spatial_autocorrelation %>%
      add_column(., Stability_var = stability_var),
    Temporal_autocorrelation = all_temporal_autocorrelation %>%
      add_column(., Stability_var = stability_var)
  ))
}


boot_function_lm_partialres = function(formula, data, indices) {
  d = data[indices,] 
  fit = lm(formula, data=d) 
  return(summary(fit)$coefficient[2,1])
}




# TRADE OFF ----


plot_tradeoff = function(d_pair, gradient_label = "LUI") {
  ggplot(d_pair,
         aes(x = Effect_size_Rs, y = Effect_size_Rl, color = Variable_gradient)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.8, size = 2) +
    facet_wrap(~Predictor, scales = "fixed", nrow = 2) +
    scale_color_viridis_c(name = gradient_label) +
    labs(x = "Effect size on drought resistance",
         y = "Effect size on drought resilience") +
    the_theme2 +
    theme(strip.background = element_rect(fill = "black"),
          strip.text       = element_text(colour = "white"),
          legend.position  = "right")
}
join_responses = function(res_obj, rel_obj) {
  
  res = res_obj$Partial_residuals %>%
    dplyr::filter(Stability_var == "Resistance_Isbell") %>%
    dplyr::select(-Stability_var) %>%
    tidyr::pivot_longer(-Variable_gradient,
                        names_to = "Predictor",
                        values_to = "Effect_size_Rs") %>%
    dplyr::filter(grepl("CWM|FD", Predictor), !is.na(Effect_size_Rs)) %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::summarise(Effect_size_Rs = mean(Effect_size_Rs, na.rm = TRUE),
                     .groups = "drop")
  
  rel = rel_obj$Partial_residuals %>%
    dplyr::filter(Stability_var == "Resilience_Isbell_abs") %>%
    dplyr::select(-Stability_var) %>%
    tidyr::pivot_longer(-Variable_gradient,
                        names_to = "Predictor",
                        values_to = "Effect_size_Rl") %>%
    dplyr::filter(grepl("CWM|FD", Predictor), !is.na(Effect_size_Rl)) %>%
    dplyr::group_by(Predictor, Variable_gradient) %>%
    dplyr::summarise(Effect_size_Rl = mean(Effect_size_Rl, na.rm = TRUE),
                     .groups = "drop")
  
  dplyr::inner_join(res, rel, by = c("Predictor", "Variable_gradient"))
}
relabel_predictors = function(d_pair) {
  d_pair %>%
    dplyr::mutate(
      Predictor = dplyr::recode(Predictor,
                                "CWM_PC1" = "CWM PC1: Economics (slow - fast)",
                                "CWM_PC2" = "CWM PC2: Absorption (AMF - hairs)",
                                "CWM_PC3" = "CWM PC3: Exploration (AD - SRL)",
                                "FD_PC1"  = "FD PC1: Economics (slow - fast)",
                                "FD_PC2"  = "FD PC2: Absorption (AMF - hairs)",
                                "FD_PC3"  = "FD PC3: Exploration (AD - SRL)"
      )
    )
}

quadrant_tab = function(d_pair) {
  d_pair %>%
    dplyr::mutate(quadrant = dplyr::case_when(
      Effect_size_Rs < 0 & Effect_size_Rl > 0 ~ "Trade-off (+ on resilience, - on resistance)",
      Effect_size_Rs > 0 & Effect_size_Rl < 0 ~ "Trade-off (- on resilience, + on resistance)",
      Effect_size_Rs > 0 & Effect_size_Rl > 0 ~ "Synergy (+ on both)",
      Effect_size_Rs < 0 & Effect_size_Rl < 0 ~ "Synergy (- on both)",
      TRUE ~ "Neutral"
    )) %>%
    dplyr::count(Predictor, quadrant) %>%
    dplyr::group_by(Predictor) %>%
    dplyr::mutate(fraction = n / sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(Predictor,
                    quadrant = c("Trade-off (+ on resilience, - on resistance)",
                                 "Trade-off (- on resilience, + on resistance)",
                                 "Synergy (+ on both)",
                                 "Synergy (- on both)"),
                    fill = list(n = 0, fraction = 0))
}


predictor_colors = c(
  "CWM PC1: Economics (slow - fast)"     = "#3D405B",
  "CWM PC2: Absorption (AMF - hairs)"    = "#E07A5F",
  "CWM PC3: Exploration (AD - SRL)"    = "#81B29A",
  "FD PC1: Economics (slow - fast)"      = "#F2CC8F",
  "FD PC2: Absorption (AMF - hairs)"     = "#9A8C98",
  "FD PC3: Exploration (AD - SRL)"     = "#22223B"
)
