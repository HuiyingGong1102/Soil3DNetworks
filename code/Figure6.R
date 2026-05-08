


library(dplyr)
library(purrr)
library(randomForest)
library(ggplot2)
library(tidytext) 

g1 <- paste0("R", c(20,16,31,4,15,24,18,23,33,9,26,17,25))
g2 <- paste0("R", c(32,1,22,6,11,27,28))
g3 <- paste0("R", c(3,30,13,21,14,19,29,7,10,2,12,5,8))

df_long <- bind_rows(df_all, .id = "Region_ID")

df_plot <- df_long %>%
  mutate(
    Regime_Class = case_when(
      Region_ID %in% g1 ~ "Regime 1",
      Region_ID %in% g2 ~ "Regime 2",
      Region_ID %in% g3 ~ "Regime 3",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Regime_Class != "Other") 

df_plot$Regime_Class <- factor(
  df_plot$Regime_Class,
  levels = c("Regime 1", 
             "Regime 2", 
             "Regime 3")
)

colnames(df_plot) <- c("Region_ID",colname,"EMF","Regime_Class")
df_plot <- df_plot[,c("Region_ID",group1_features, group2_features, group3_features, group4_features,"EMF","Regime_Class")]


library(dplyr)
library(purrr)
library(randomForest)
library(ggplot2)
library(tidytext)
library(viridis)

set.seed(2026)

n_repeat <- 30

rf_results <- df_plot %>%
  dplyr::select(-Region_ID) %>%
  na.omit() %>%
  group_split(Regime_Class) %>%
  set_names(map_chr(., ~ as.character(unique(.x$Regime_Class)))) %>%
  map(function(group_data) {
    
    model_data <- group_data %>% dplyr::select(-Regime_Class)
    
    imp_all <- map_dfr(1:n_repeat, function(i) {
      set.seed(2026 + i)
      
      rf_model <- randomForest(
        EMF ~ .,
        data = model_data,
        ntree = 1000,
        importance = TRUE
      )
      
      imp_df <- as.data.frame(importance(rf_model))
      imp_df$Feature <- rownames(imp_df)
      imp_df$Repeat <- i
      
      imp_df
    })
    
    imp_summary <- imp_all %>%
      group_by(Feature) %>%
      summarise(
        mean_IncMSE = mean(`%IncMSE`, na.rm = TRUE),
        sd_IncMSE = sd(`%IncMSE`, na.rm = TRUE),
        .groups = "drop"
      )
    
    top10_freq_df <- imp_all %>%
      group_by(Repeat) %>%
      arrange(desc(`%IncMSE`), .by_group = TRUE) %>%
      slice_head(n = 10) %>%
      ungroup() %>%
      count(Feature, name = "top10_count") %>%
      mutate(top10_freq = top10_count / n_repeat)
    
    final_summary <- imp_summary %>%
      left_join(top10_freq_df, by = "Feature") %>%
      mutate(
        top10_count = ifelse(is.na(top10_count), 0, top10_count),
        top10_freq = ifelse(is.na(top10_freq), 0, top10_freq)
      ) %>%
      arrange(desc(mean_IncMSE)) %>%
      slice(1:10) %>%
      mutate(Regime = unique(group_data$Regime_Class))
    
    return(final_summary)
  })

all_top_features <- bind_rows(rf_results)

p_rf_compare <- ggplot(
  all_top_features,
  aes(x = reorder_within(Feature, mean_IncMSE, Regime), y = mean_IncMSE)
) +
  geom_segment(
    aes(
      xend = reorder_within(Feature, mean_IncMSE, Regime),
      y = 0,
      yend = mean_IncMSE
    ),
    color = "grey60",
    linewidth = 1
  ) +
  geom_point(
    aes(color = mean_IncMSE, size = top10_freq),
    alpha = 0.95
  ) +
  facet_wrap(~ Regime, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(2.5, 5)) +
  theme_bw(base_size = 13) +
  labs(
    x = "Network topological features",
    y = "Mean importance (% increase in MSE)",
    size = "Top-10 frequency"
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.y = element_text(color = "black", face = "bold", size = 9)
  )

print(p_rf_compare)








library(dplyr)
library(purrr)
library(randomForest)
library(pdp)
library(ggplot2)

top_features_list <- list(
  "Regime 1" = c("withinL3_pos_beta0", "betweenLag0_pos_beta0"),
  "Regime 2" = c("withinL6_neg_beta0", "betweenLag0_pos_beta0"),
  "Regime 3" = c("betweenLag1_pos_beta0", "betweenLag0_pos_beta0")
)

set.seed(2026)

pdp_data_all <- map_dfr(names(top_features_list), function(regime_name) {
  
  group_data <- df_plot %>% 
    filter(Regime_Class == regime_name) %>%
    dplyr::select(EMF, withinL1_pos_beta0:betweenLag3_neg_beta0) %>% 
    na.omit()
  
  rf_model <- randomForest(
    EMF ~ ., 
    data = group_data, 
    ntree = 1000,
    importance = TRUE
  )
  
  feature_vec <- top_features_list[[regime_name]]
  
  map_dfr(feature_vec, function(feature_name) {
    
    pd <- partial(
      rf_model,
      pred.var = feature_name,
      train = group_data
    )
    
    pd_clean <- pd %>%
      rename(Feature_Value = all_of(feature_name)) %>%
      mutate(
        Regime = regime_name,
        Feature_Name = feature_name
      )
    
    return(pd_clean)
  })
})

pdp_data_all <- pdp_data_all %>%
  mutate(
    Facet_ID = paste(Regime, Feature_Name, sep = "__")
  )

facet_levels <- unlist(
  lapply(names(top_features_list), function(rg) {
    paste(rg, top_features_list[[rg]], sep = "__")
  })
)

pdp_data_all$Facet_ID <- factor(pdp_data_all$Facet_ID, levels = facet_levels)

p_pdp <- ggplot(pdp_data_all, aes(x = Feature_Value, y = yhat)) +
  
  geom_smooth(
    aes(color = Regime, fill = Regime),
    method = "loess",
    span = 0.6,
    alpha = 0.2,
    linewidth = 1.2,
    se = TRUE
  ) +
  
  geom_line(
    aes(color = Regime),
    linewidth = 0.7,
    linetype = "dashed",
    alpha = 0.6
  ) +
  
  facet_wrap(
    ~ Facet_ID,
    scales = "free",
    nrow = 2,
    ncol = 3,
    dir = "v",
    labeller = as_labeller(function(x) sub("^.*__", "", x))
  ) +
  
  scale_color_manual(values = c(
    "Regime 1" = "firebrick",
    "Regime 2" = "goldenrod",
    "Regime 3" = "forestgreen"
  )) +
  scale_fill_manual(values = c(
    "Regime 1" = "firebrick",
    "Regime 2" = "goldenrod",
    "Regime 3" = "forestgreen"
  )) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "",
    x = "Value of the keystone topological feature",
    y = "Predicted ecosystem multifunctionality (EMF)"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 0, hjust = 1)
  )

print(p_pdp)










library(dplyr)
library(purrr)
library(randomForest)
library(ggplot2)
library(tidytext) 

set.seed(2026)

n_repeat <- 30

rf_results <- df_plot %>%
  dplyr::select(-Region_ID) %>%
  na.omit() %>%
  group_split(Regime_Class) %>%
  set_names(map_chr(., ~ as.character(unique(.x$Regime_Class)))) %>%
  map(function(group_data) {
    
    model_data <- group_data %>% dplyr::select(-Regime_Class)
    
    imp_all <- map_dfr(1:n_repeat, function(i) {
      set.seed(2026 + i)
      
      rf_model <- randomForest(
        EMF ~ .,
        data = model_data,
        ntree = 1000,
        importance = TRUE
      )
      
      imp_df <- as.data.frame(importance(rf_model))
      imp_df$Feature <- rownames(imp_df)
      imp_df$Repeat <- i
      
      imp_df
    })
    
    imp_summary <- imp_all %>%
      group_by(Feature) %>%
      summarise(
        mean_IncMSE = mean(`%IncMSE`, na.rm = TRUE),
        sd_IncMSE = sd(`%IncMSE`, na.rm = TRUE),
        .groups = "drop"
      )
    
    top10_freq_df <- imp_all %>%
      group_by(Repeat) %>%
      arrange(desc(`%IncMSE`), .by_group = TRUE) %>%
      slice_head(n = 10) %>%
      ungroup() %>%
      count(Feature, name = "top10_count") %>%
      mutate(top10_freq = top10_count / n_repeat)
    
    final_summary <- imp_summary %>%
      left_join(top10_freq_df, by = "Feature") %>%
      mutate(
        top10_count = ifelse(is.na(top10_count), 0, top10_count),
        top10_freq = ifelse(is.na(top10_freq), 0, top10_freq)
      ) %>%
      arrange(desc(mean_IncMSE)) %>%
      slice(1:20) %>%
      mutate(Regime = unique(group_data$Regime_Class))
    
    return(final_summary)
  })


Reduce(intersect, list(rf_results[[1]]$Feature,rf_results[[2]]$Feature,rf_results[[3]]$Feature))


library(lavaan)
library(semPlot)

topo_sem_model <- '
  withinL1_pos_beta0 ~ withinL2_pos_beta0
  withinL2_pos_beta0 ~ withinL3_pos_beta0
  withinL3_pos_beta0 ~ withinL5_pos_beta1
  withinL5_neg_beta0 ~ withinL6_neg_beta0
  
  betweenLag1_pos_beta0 ~ betweenLag0_pos_beta0
  betweenLag2_pos_beta0 ~ betweenLag1_pos_beta0

  betweenLag0_pos_beta0 ~ withinL2_pos_beta0 + withinL3_pos_beta0
  betweenLag1_pos_beta0 ~ withinL3_pos_beta0 + withinL5_pos_beta1
  betweenLag2_pos_beta0 ~ withinL6_neg_beta0 + withinL5_neg_beta0 + withinL5_pos_beta1

  EMF ~ betweenLag2_pos_beta0 + withinL2_pos_beta0 + withinL1_pos_beta0 + 
  betweenLag1_pos_beta0+ betweenLag0_pos_beta0 + withinL3_pos_beta0+withinL5_pos_beta1+
  withinL5_neg_beta0+withinL6_neg_beta0

  withinL5_pos_beta1 ~~ withinL5_neg_beta0
'

df_sem <- df_plot %>%
  dplyr::select(
    EMF,
    withinL1_pos_beta0,
    withinL2_pos_beta0,
    withinL3_pos_beta0,
    withinL5_pos_beta1,
    withinL5_neg_beta0,
    withinL6_neg_beta0,
    betweenLag0_pos_beta0,
    betweenLag1_pos_beta0,
    betweenLag2_pos_beta0,
    Regime_Class
  ) %>%
  na.omit()

df_plot_scaled <- as.data.frame(
  scale(df_sem %>% dplyr::select(-Regime_Class))
)

df_plot_scaled$Regime_Class <- df_sem$Regime_Class

fit_strict <- sem(topo_sem_model, data = df_plot_scaled, estimator = "MLM",auto.cov.y = FALSE)
summary(fit_strict, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_strict, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "df"))

mi <- modindices(fit_strict)

top_missing_paths <- mi %>%
  filter(op == "~") %>%          
  arrange(desc(mi)) %>%          
  filter(mi > 10) %>%            
  dplyr::select(Target = lhs, Predictor = rhs, MI_Score = mi, Expected_Change = epc)

print(top_missing_paths)

revised_model <- '
  withinL1_pos_beta0 ~ withinL2_pos_beta0+withinL3_pos_beta0+withinL6_neg_beta0+withinL5_neg_beta0
  withinL2_pos_beta0 ~ withinL3_pos_beta0+withinL6_neg_beta0
  withinL3_pos_beta0 ~ withinL5_pos_beta1+withinL6_neg_beta0+withinL5_neg_beta0
  withinL5_neg_beta0 ~ withinL6_neg_beta0
  withinL1_pos_beta0 ~ betweenLag0_pos_beta0
  
  betweenLag1_pos_beta0 ~ betweenLag0_pos_beta0
  betweenLag2_pos_beta0 ~ betweenLag1_pos_beta0
  betweenLag0_pos_beta0 ~ withinL5_pos_beta1

  betweenLag0_pos_beta0 ~ withinL2_pos_beta0 + withinL3_pos_beta0
  betweenLag1_pos_beta0 ~ withinL3_pos_beta0 + withinL5_pos_beta1
  betweenLag2_pos_beta0 ~ withinL6_neg_beta0 + withinL5_neg_beta0 + withinL5_pos_beta1

  EMF ~ betweenLag2_pos_beta0 + withinL2_pos_beta0 + withinL1_pos_beta0 + 
  betweenLag1_pos_beta0+ betweenLag0_pos_beta0 + withinL3_pos_beta0+withinL5_pos_beta1+
  withinL5_neg_beta0+withinL6_neg_beta0

  withinL5_pos_beta1 ~~ withinL5_neg_beta0
'

fit_revised <- sem(revised_model, data = df_plot_scaled, estimator = "MLM",auto.cov.y = FALSE)
summary(fit_revised, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_revised, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "df"))

mi <- modindices(fit_revised)

top_missing_paths <- mi %>%
  filter(op == "~") %>%          
  arrange(desc(mi)) %>%          
  filter(mi > 10) %>%            
  dplyr::select(Target = lhs, Predictor = rhs, MI_Score = mi, Expected_Change = epc)

print(top_missing_paths)



parameterEstimates(fit_revised, standardized = TRUE) %>%
  dplyr::filter(op == "~") %>%
  dplyr::select(lhs, rhs, est, std.all, pvalue)


final_model <- '
  withinL1_pos_beta0 ~ withinL2_pos_beta0+withinL3_pos_beta0+withinL6_neg_beta0+withinL5_neg_beta0
  withinL2_pos_beta0 ~ withinL3_pos_beta0+withinL6_neg_beta0
  withinL3_pos_beta0 ~ withinL6_neg_beta0+withinL5_neg_beta0
  withinL5_neg_beta0 ~ withinL6_neg_beta0
  withinL1_pos_beta0 ~ betweenLag0_pos_beta0
  
  betweenLag1_pos_beta0 ~ betweenLag0_pos_beta0
  betweenLag2_pos_beta0 ~ betweenLag1_pos_beta0
  betweenLag0_pos_beta0 ~ withinL5_pos_beta1

  # 3. within -> between
  betweenLag0_pos_beta0 ~ withinL2_pos_beta0 + withinL3_pos_beta0

  EMF ~ withinL1_pos_beta0 + 
  betweenLag1_pos_beta0+ withinL3_pos_beta0+withinL5_neg_beta0+withinL6_neg_beta0

  withinL5_pos_beta1 ~~ withinL5_neg_beta0
'

fit_final <- sem(final_model, data = df_plot_scaled, estimator = "MLM",auto.cov.y = FALSE)
summary(fit_final, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_final, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr", "df"))

mi <- modindices(fit_final)

top_missing_paths <- mi %>%
  filter(op == "~") %>%         
  arrange(desc(mi)) %>%          
  filter(mi > 10) %>%            
  dplyr::select(Target = lhs, Predictor = rhs, MI_Score = mi, Expected_Change = epc)

print(top_missing_paths)


final_paths <- parameterEstimates(fit_final, standardized = TRUE) %>%
  dplyr::filter(op == "~") %>%
  dplyr::select(
    Target = lhs,
    Predictor = rhs,
    Estimate = est,
    Std_Coef = std.all,
    P_value = pvalue
  ) %>%
  mutate(
    Significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

final_paths


fit_mg_final <- sem(final_model, 
                    data = df_plot_scaled, 
                    group = "Regime_Class", 
                    estimator = "MLM",auto.cov.y = FALSE)
fitMeasures(fit_mg_final, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr"))


mg_paths <- parameterEstimates(fit_mg_final, standardized = TRUE) %>%
  dplyr::filter(op == "~") %>%
  dplyr::select(group, lhs, rhs, est, std.all, pvalue)

mg_paths
