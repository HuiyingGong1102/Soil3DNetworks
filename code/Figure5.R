
library(jsonlite)
filename <- "./json/k33-1.json"
json_data <- fromJSON(filename)
module<-json_data$max_omega_logi
table(module)
length(table(module))
range(table(module))
k <- 33

json_data <- fromJSON("./json/k434-1.json")
ki <- 434
modulei<-json_data$max_omega_logi

datt <- rbind(depth1,depth2,depth3,depth4,depth5,depth6)
select_datt<- datt[,c("BD","CEC","SOC","TK","TN","TP")]

scale_datt <- apply(select_datt, 2, scale)
scale_datt[,1] <- -1*scale_datt[,1]

emf1 <- rowMeans(scale_datt)
matrix_emf1 <- matrix(emf1, nrow = 279120)

emf2 <- rowMeans(matrix_emf1)


EMF_all <- list()
for (i in 1:k) {
  EMF <- c()
  index<-which(module==(((1:k))-1)[i])
  for (ii in index) {
    EMF <- c(EMF,mean(emf2[which(modulei==(((1:ki))-1)[ii])]))
  }
  EMF_all[[i]] <- EMF
}
names(EMF_all) <- paste0("R", 1:33)



df_all <- list()
for (i2 in 1:33) {
  i=i2
  EMF <- c()
  index<-which(module==(((1:k))-1)[i])
  for (ii in index) {
    feature_all <- result_mat[index,]
    EMF <- c(EMF,mean(emf2[which(modulei==(((1:ki))-1)[ii])]))
  }
  df <- feature_all
  df <- as.data.frame(df)
  df$EMF <- EMF
  df_all[[i2]] <- df
}


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

names(EMF_all) <- paste0("R", 1:33)

df_long <- tibble::enframe(EMF_all, name = "Region_ID", value = "EMF") %>%
  unnest(EMF)

head(df_long)

kw_result <- kruskal.test(EMF ~ Region_ID, data = df_long)

print(kw_result)


df_long$Region_ID <- reorder(df_long$Region_ID, df_long$EMF, FUN = median)

p_boxplot <- ggplot(df_long, aes(x = Region_ID, y = EMF, fill = Region_ID)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    legend.position = "none", 
    panel.grid.major.x = element_blank() 
  ) +
  labs(
    # title = "Variation of Ecosystem Multifunctionality Across 33 Regions",
    x = "Region (Ordered by Median EMF)",
    y = "Comprehensive EMF"
  ) +
  stat_compare_means(
    method = "kruskal.test", 
    # label.y.npc = "top",    
    label.y.npc = 0.96,
    label.x.npc = 0.02,     
    hjust = 0,            
    size = 4
  )




library(dplyr)
library(purrr) 



cor_matrix <- c()
for (i in 1:33) {
  df=df_all[[i]]
  
  cor_results <- df %>%
    dplyr::select(where(is.numeric)) %>% 
    dplyr::select(-EMF) %>% 
    purrr::map_dbl(~ cor(.x, df$EMF, use = "pairwise.complete.obs", method = "spearman"))
  
  cor_matrix <- cbind(cor_matrix,cor_results)
}

prefixes <- c(
  paste0("withinL", 1:6),paste0("betweenLag", 0:4)
)

types <- c("pos", "neg")
var_names <- as.vector(outer(prefixes, types, paste, sep = "_"))
colname <- c()
for (nm in var_names) {
  colname <- c(colname,paste0(nm,"_beta", 0:3))
}
rownames(cor_matrix) <- colname

cor_matrix[is.na(cor_matrix)] <- 0
cor_matrix <- cor_matrix[rowMeans(cor_matrix == 0, na.rm = TRUE) < 0.95, ]
dim(cor_matrix)
colnames(cor_matrix) <- paste0("R", 1:33) 
rownames(cor_matrix)


library(pheatmap)

group1_features <- c("withinL1_pos_beta0","withinL1_pos_beta1","withinL1_pos_beta3",
                     "withinL2_pos_beta0","withinL2_pos_beta1","withinL2_pos_beta3",
                     "withinL3_pos_beta0","withinL3_pos_beta1","withinL3_pos_beta2","withinL3_pos_beta3",
                     "withinL4_pos_beta0","withinL4_pos_beta1","withinL4_pos_beta2","withinL4_pos_beta3",
                     "withinL5_pos_beta0","withinL5_pos_beta1","withinL5_pos_beta2","withinL5_pos_beta3",
                     "withinL6_pos_beta0","withinL6_pos_beta1","withinL6_pos_beta2","withinL6_pos_beta3")

group2_features <- c("withinL1_neg_beta0","withinL1_neg_beta1","withinL1_neg_beta3",
                     "withinL2_neg_beta0","withinL2_neg_beta1",
                     "withinL3_neg_beta0","withinL3_neg_beta1",
                     "withinL4_neg_beta0","withinL4_neg_beta1",
                     "withinL5_neg_beta0",
                     "withinL6_neg_beta0","withinL6_neg_beta1")

group3_features <- c("betweenLag0_pos_beta0","betweenLag0_pos_beta1","betweenLag0_pos_beta2",
                     "betweenLag1_pos_beta0",
                     "betweenLag2_pos_beta0",
                     "betweenLag3_pos_beta0",
                     "betweenLag4_pos_beta0")

group4_features <- c("betweenLag0_neg_beta0","betweenLag0_neg_beta1",
                     "betweenLag1_neg_beta0",
                     "betweenLag2_neg_beta0",
                     "betweenLag3_neg_beta0")

ordered_features <- c(group1_features, group2_features, group3_features, group4_features)
cor_matrix_ordered <- cor_matrix[ordered_features, ]

col_annotation <- data.frame(
  Feature_Category = c(rep("Horizontal-Pos", length(group1_features)),
                       rep("Horizontal-Neg", length(group2_features)),
                       rep("Vertical-Pos", length(group3_features)),
                       rep("Vertical-Neg", length(group4_features)))
)
rownames(col_annotation) <- ordered_features 

ann_colors <- list(
  Feature_Category = c(
    "Horizontal-Pos" = "#ed3525",
    "Horizontal-Neg" = "#5e9bd6",
    "Vertical-Pos" = "#ed3525",
    "Vertical-Neg" = "#5e9bd6"
  )
)





pheatmap(cor_matrix_ordered,
         cluster_cols = TRUE,               
         cluster_rows = FALSE,                
         annotation_row = col_annotation,    
         annotation_names_row = FALSE,
         annotation_legend = FALSE,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         breaks = seq(-1, 1, length.out = 51), 
         clustering_distance_cols = "correlation",
         cutree_cols = 3,
         clustering_method = "ward.D2",
         fontsize_col = 10.5,
         fontsize_row = 12,
         angle_col = 45,
         border_color = NA) 




library(dplyr)
library(ggplot2)

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
  )

df_plot$Regime_Class <- factor(
  df_plot$Regime_Class,
  levels = c("Regime 1", 
             "Regime 2", 
             "Regime 3")
)

colnames(df_plot) <- c("Region_ID",colname,"EMF","Regime_Class")
df_plot <- df_plot[,c("Region_ID",group1_features, group2_features, group3_features, group4_features,"EMF","Regime_Class")]


for (i in 2:47) {
  point_colors <- c(
    "Regime 1" = "darkgreen",    
    "Regime 2" = "#ff7000",
    "Regime 3" = "darkred"
  )
  line_colors <- c(
    "Regime 1" = "darkgreen",    
    "Regime 2" = "#ff7000",
    "Regime 3" = "darkred"
  )
  
  p_reversal <- ggplot(df_plot, aes(x = EMF, y = df_plot[, i])) +
    geom_point(aes(color = Regime_Class), alpha = 0.2, size = 1.8) +
    
    geom_smooth(aes(color = Regime_Class, fill = Regime_Class), 
                method = "lm", se = TRUE, linewidth = 1.2,alpha = 0.3) +
    
    facet_wrap(~ Regime_Class, scales = "free_x") +
    
    scale_color_manual(values = point_colors) +
    scale_fill_manual(values = line_colors) +
    
    theme_bw(base_size = 14) +
    labs(x = "EMF", y = colnames(df_plot)[i]) +
    theme(
      legend.position = "none",          
      strip.text = element_blank(),      
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1)
    )
  print(p_reversal)
}
