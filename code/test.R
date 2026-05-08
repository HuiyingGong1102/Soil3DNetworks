library(orthopolynom)
library(parallel)
library(pbapply)
library(grpreg)
library(glmnet)
library(polynom)
library(caret)
library(jsonlite)
library(reshape2)
library(dplyr)

data <- readRDS("D:/Soil3DNetworks/data/matrix_all.rds")
dim(data)

cleaned_data <- data[!duplicated(data), ]
dim(cleaned_data)

############################################

depth1 <- cleaned_data[,1:15]
colnames(depth1) <- c("BD","Clay","Silt","Sand","CEC",
                      "CF","pH","SOC","SOCD","TK",
                      "TKD","TN","TND","TP","TPD")

Normalized_depth1 <- t(sweep(depth1, 2, apply(depth1, 2, max), FUN = "/"))

total <- apply(Normalized_depth1,1,sum)
index1 <- total[order(total)]
data_order1 <- Normalized_depth1[order(total),]

fit1 <- power_equation_fit2(
  x = index1,
  y = data_order1,
  X_smooth = seq(index1[1], index1[length(index1)], length.out = 15),
  thread = 24)

Normalized_depth1 <- Normalized_depth1+1
total <- apply(Normalized_depth1,1,sum)
log_data_order1 <- log(Normalized_depth1[order(total),])
log_index1 <- log(total[order(total)])


###########################################

depth2 <- cleaned_data[,16:30]
colnames(depth2) <- c("BD","Clay","Silt","Sand","CEC",
                      "CF","pH","SOC","SOCD","TK",
                      "TKD","TN","TND","TP","TPD")

Normalized_depth2 <- t(sweep(depth2, 2, apply(depth2, 2, max), FUN = "/"))

total <- apply(Normalized_depth2,1,sum)
index2 <- total[order(total)]
data_order2 <- Normalized_depth2[order(total),]

fit2 <- power_equation_fit2(
  x = index2,
  y = data_order2,
  X_smooth = seq(index2[1], index2[length(index2)], length.out = 15),
  thread = 8)


Normalized_depth2 <- Normalized_depth2+1
total <- apply(Normalized_depth2,1,sum)
log_data_order2 <- log(Normalized_depth2[order(total),])
log_index2 <- log(total[order(total)])


###########################################


depth3 <- cleaned_data[,31:45]
colnames(depth3) <- c("BD","Clay","Silt","Sand","CEC",
                      "CF","pH","SOC","SOCD","TK",
                      "TKD","TN","TND","TP","TPD")

Normalized_depth3 <- t(sweep(depth3, 2, apply(depth3, 2, max), FUN = "/"))

total <- apply(Normalized_depth3,1,sum)
index3 <- total[order(total)]
data_order3 <- Normalized_depth3[order(total),]

fit3 <- power_equation_fit2(
  x = index3,
  y = data_order3,
  X_smooth = seq(index3[1], index3[length(index3)], length.out = 15),
  thread = 24)

Normalized_depth3 <- Normalized_depth3+1
total <- apply(Normalized_depth3,1,sum)
log_data_order3 <- log(Normalized_depth3[order(total),])
log_index3 <- log(total[order(total)])


###########################################

depth4 <- cleaned_data[,46:60]
colnames(depth4) <- c("BD","Clay","Silt","Sand","CEC",
                      "CF","pH","SOC","SOCD","TK",
                      "TKD","TN","TND","TP","TPD")

Normalized_depth4 <- t(sweep(depth4, 2, apply(depth4, 2, max), FUN = "/"))

total <- apply(Normalized_depth4,1,sum)
index4 <- total[order(total)]
data_order4 <- Normalized_depth4[order(total),]

fit4 <- power_equation_fit2(
  x = index4,
  y = data_order4,
  X_smooth = seq(index4[1], index4[length(index4)], length.out = 15),
  thread = 24)

Normalized_depth4 <- Normalized_depth4+1
total <- apply(Normalized_depth4,1,sum)
log_data_order4 <- log(Normalized_depth4[order(total),])
log_index4 <- log(total[order(total)])


###########################################

depth5 <- cleaned_data[,61:75]
colnames(depth5) <- c("BD","Clay","Silt","Sand","CEC",
                      "CF","pH","SOC","SOCD","TK",
                      "TKD","TN","TND","TP","TPD")

Normalized_depth5 <- t(sweep(depth5, 2, apply(depth5, 2, max), FUN = "/"))

total <- apply(Normalized_depth5,1,sum)
index5 <- total[order(total)]
data_order5 <- Normalized_depth5[order(total),]

fit5 <- power_equation_fit2(
  x = index5,
  y = data_order5,
  X_smooth = seq(index5[1], index5[length(index5)], length.out = 15),
  thread = 24)

Normalized_depth5 <- Normalized_depth5+1
total <- apply(Normalized_depth5,1,sum)
log_data_order5 <- log(Normalized_depth5[order(total),])
log_index5 <- log(total[order(total)])


###########################################

depth6 <- cleaned_data[,76:90]
colnames(depth6) <- c("BD","Clay","Silt","Sand","CEC",
                      "CF","pH","SOC","SOCD","TK",
                      "TKD","TN","TND","TP","TPD")

Normalized_depth6 <- t(sweep(depth6, 2, apply(depth6, 2, max), FUN = "/"))

total <- apply(Normalized_depth6,1,sum)
index6 <- total[order(total)]
data_order6 <- Normalized_depth6[order(total),]

fit6 <- power_equation_fit2(
  x = index6,
  y = data_order6,
  X_smooth = seq(index6[1], index6[length(index6)], length.out = 15),
  thread = 24)

Normalized_depth6 <- Normalized_depth6+1
total <- apply(Normalized_depth6,1,sum)
log_data_order6 <- log(Normalized_depth6[order(total),])
log_index6 <- log(total[order(total)])



###########################################

datt <- rbind(depth1,depth2,depth3,depth4,depth5,depth6)
Normalized_depth <- t(sweep(datt, 2, apply(datt, 2, max), FUN = "/"))

data_all <- rbind(Normalized_depth[,1:279120],Normalized_depth[,(1+279120):(2*279120)],
                  Normalized_depth[,(1+2*279120):(3*279120)],Normalized_depth[,(1+3*279120):(4*279120)],
                  Normalized_depth[,(1+4*279120):(5*279120)],Normalized_depth[,(1+5*279120):(6*279120)])
dim(data_all)
total <- apply(data_all,1,sum)
index <- total[order(total)]
data_order <- data_all[order(total),]

Normalized_depth <- data_all+1
total <- apply(Normalized_depth,1,sum)
log_data_order <- log(Normalized_depth[order(total),])
log_index <- log(total[order(total)])



