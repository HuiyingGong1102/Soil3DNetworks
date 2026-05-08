library(orthopolynom)
library(parallel)
library(pbapply)
library(grpreg)
library(glmnet)
library(polynom)
library(caret)
library(gglasso)
library(ggplot2)
library(reshape2)
library(readxl)
library(dplyr)

distance <- c(-(0- 5)/2, 5-(5- 15)/2, 15-(15- 30)/2,
              30-(30- 60)/2, 60-(60- 100)/2, 100-(100- 200)/2)

datt <- rbind(depth1,depth2,depth3,depth4,depth5,depth6)
Normalized_depth <- t(sweep(datt, 2, apply(datt, 2, max), FUN = "/"))

data_all <- rbind(Normalized_depth[,1:279120],Normalized_depth[,(1+279120):(2*279120)],
                  Normalized_depth[,(1+2*279120):(3*279120)],Normalized_depth[,(1+3*279120):(4*279120)],
                  Normalized_depth[,(1+4*279120):(5*279120)],Normalized_depth[,(1+5*279120):(6*279120)])
library(jsonlite)
json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer_all/json/k434-1.json")
k <- 434
module<-json_data$max_omega_logi
length(table(module))
table(module)

for (ii in 1:k) {
  group <- which(module==(ii-1))
  
  net_data <- data_all[,group]
  
  net_data <- t(net_data)
  total <- apply(net_data,1,sum)
  index <- total[order(total)]
  data_order <- net_data[order(total),]
  
  fit <- power_equation_fit2(
    x = index,
    y = data_order,
    X_smooth = seq(index[1], index[length(index)], length.out = 200),
    thread = 24)
  
  power_par <- fit$power_par
  # time <- fit$Time
  time <- fit$smooth_Time
  smooth_time <- fit$smooth_Time
  effecti=fit$smooth_fit
  colnames(effecti) <- as.vector(t(outer((1:6), 1:15, 
                                         FUN = function(x, y) paste(x, y, sep = "-"))))
  
  
  relationship_all <- c()
  for (i in 1:6) {
    relationship <- pblapply(1:15,function(c) get_interaction(effecti, col=c, depth=i))
    relationship_all <- c(relationship_all,relationship)
  }
  
  
  order_ind <- 0
  order_dep <- 0
  cl <- makeCluster(getOption("cl.cores", 24))
  clusterEvalQ(cl, {require(orthopolynom)})
  clusterExport(cl, c("power_par","order_ind","order_dep","relationship_all",
                      "get_value","ode_optimize","get_effect","effecti",
                      "time","smooth_time","power.equation"),
                envir=environment())
  lop_par <- pblapply(1:90,function(c)
    get_value(effecti,power_par,time,smooth_time,relationship_all[[c]],order_ind,order_dep),cl=cl)
  stopCluster(cl)
  
  
  links_table <- c()
  ij <- 1
  for (i in 1:6) {
    for (j in 1:15) {
      
      d2 <- get_ode_output(relationship=relationship_all[[ij]],
                           par=lop_par[[ij]],
                           effecti=effecti,
                           power_par=power_par,
                           time=smooth_time,
                           order_ind=0,
                           order_dep=0)
      
      for (depi in 1:length(d2[[2]])) {
        links_table <- rbind(links_table,c(d2[[2]][depi],d2[[1]],
                                           as.numeric(d2[[5]][1+depi])))
      }
      
      ij <- ij+1
    }
  }
  
  colnames(links_table)<-c("source","target","weight")
  links_table<-as.data.frame(links_table)
  links_table$edge_type<-ifelse(links_table$weight>0,1,2)
  links_table$weight <- abs(as.numeric(links_table$weight))
  rownames(links_table)<-rep(1:nrow(links_table))
  links_table$weight <- as.numeric(links_table$weight)^0.1*5
  
  library(openxlsx)
  
  linkname <- paste("./network/links_",ii,".xlsx",sep = "")
  write.xlsx(links_table, linkname)
  
  
  
  
  
  x <- c(900,372,900,768,504,636,240,504,636,768,636,372,504,1032,768)
  y <- c(-490,-280,-385,-280,-280,-280,-280,-385,-385,-385,-490,-385,-490,-490,-490)
  
  
  x_all <- rep(x,times=6)
  y_all <- c(y,y+500,y+500*2,y+500*3,y+500*4,y+500*5)
  
  name <- c(15,14,13,12,11,10,9,8,7,6,3,5,2,1,4)
  s_name <- c("TPD","TP","TND","TN","TKD","TK","SOCD","SOC","pH","CF","Silt","CEC",
              "Clay","BD","Sand")
  
  node_tableall <- data.frame(c(paste(1,name, sep = "-"),
                                paste(2,name, sep = "-"),
                                paste(3,name, sep = "-"),
                                paste(4,name, sep = "-"),
                                paste(5,name, sep = "-"),
                                paste(6,name, sep = "-")),rep(s_name,6))
  colnames(node_tableall) <- c("shared name","pname")
  node_tableall$x_location <- x_all
  node_tableall$y_location <- y_all
  
  
  size_table <- c()
  ij <- 1
  for (i in 1:6) {
    for (j in 1:15) {
      d2 <- get_ode_output(relationship=relationship_all[[ij]],
                           par=lop_par[[ij]],
                           effecti=effecti,
                           power_par=power_par,
                           time=smooth_time,
                           order_ind=0,
                           order_dep=0)
      
      size_table <- rbind(size_table,c(d2[[1]],as.numeric(d2[[5]][1])))
      
      ij <- ij+1
    }
  }
  
  colnames(size_table)<-c("shared name","size")
  size_table<-as.data.frame(size_table)
  size_table$size <- as.numeric(size_table$size)^0.06*65
  
  merged_table <- merge(size_table, node_tableall, by = "shared name")
  head(merged_table)
  nodename <- paste("./network/nodes_",ii,".xlsx",sep = "")
  write.xlsx(merged_table, nodename)
  
  
  
  
  
  
  
  
  links <- links_table
  sepname <- as.vector(t(outer((1:6), 1:15, 
                               FUN = function(x, y) paste(x, y, sep = "-"))))
  
  matrix1 <- matrix(0.00001,nrow=2,ncol = 90)
  colnames(matrix1) <- 1:90
  for (i in 1:ncol(matrix1)) {
    for (j in 1:2) {
      nn <- which(links$source==sepname[i]&links$edge_type==j)
      if(length(nn)>0){
        matrix1[j,i] <- length(nn)
      }
    }
  }
  # matrix1 <- cbind(type=1:2,matrix1)
  df_j1<-melt(matrix1,id.vars='type', na.rm = FALSE)
  df_j1 <- data.frame(type=df_j1$Var1,sep=df_j1$Var2,val=df_j1$value)
  df_j1$direction=rep("outgoing")
  df_j1$depth <- as.numeric(rep(substr(sepname,1,1),each=2))
  df_j1$property <- rep(c("BD","Clay","Silt","Sand","CEC",
                          "CF","pH","SOC","SOCD","TK",
                          "TKD","TN","TND","TP","TPD")[as.numeric(substr(sepname, 3, nchar(sepname)))],each=2)
  # df_j1$x <- rep(as.numeric(substr(sepname, 3, nchar(sepname))),each=2)
  
  
  matrix2 <- matrix(0.00001,nrow=2,ncol = 90)
  colnames(matrix2) <- 1:90
  for (i in 1:ncol(matrix2)) {
    for (j in 1:2) {
      nn <- which(links$target==sepname[i]&links$edge_type==j)
      if(length(nn)>0){
        matrix2[j,i] <- -length(nn)
      }
    }
  }
  df_j2<-melt(matrix2,id.vars='type')
  df_j2 <- data.frame(type=df_j2$Var1,sep=df_j2$Var2,val=df_j2$value)
  df_j2$direction=rep("incoimg")
  df_j2$depth <- as.numeric(rep(substr(sepname,1,1),each=2))
  df_j2$property <- rep(c("BD","Clay","Silt","Sand","CEC",
                          "CF","pH","SOC","SOCD","TK",
                          "TKD","TN","TND","TP","TPD")[as.numeric(substr(sepname, 3, nchar(sepname)))],each=2)
  # df_j2$x <- rep(as.numeric(substr(sepname, 3, nchar(sepname))),each=2)
  
  merged_data <- merge(df_j1, df_j2, by = c("sep", "type","depth","property"), suffixes = c("_1", "_2"))
  merged_data$val_sum <- abs(merged_data$val_1) + abs(merged_data$val_2)
  
  data_order <- names(sort(tapply(abs(merged_data$val_sum), merged_data$property, sum),decreasing = TRUE))
  
  df <- rbind(df_j1,df_j2)
  
  df_summary <- df %>%
    group_by(direction, property,type,sep) %>%
    summarise(total_val = sum(val), .groups = "drop")
  df <- df_summary
  df <- as.data.frame(df)
  
  df$direction <- factor(df$direction,levels=c("outgoing","incoimg"))
  df$property <- factor(df$property,levels=data_order)
  
  
  
  
  
  
  ggplot() +
    geom_bar(data = df, 
             mapping = aes(x = property, y = total_val, fill = as.factor(type)), 
             stat="identity", 
             position='stack') + 
    geom_hline(yintercept = 0, linewidth = 0.5,col="black")+
    xlab("Nodes") + ylab("Frequency")+
    guides(fill="none",alpha="none")+
    scale_y_continuous(limits = c(-60,180),
                       expand = c(0, 0),
                       breaks = seq(-60,180,30),
                       labels = as.character(abs(seq(-60,180,30))))+
    scale_fill_manual(values=c(rgb(255,51,51,230,maxColorValue = 255),
                               rgb(0,153,255,230,maxColorValue = 255)))+
    theme_minimal()+
    theme(panel.background = element_rect(linewidth=1, color = 'black', fill = 'transparent'),
          axis.title.x = element_text(size=15,vjust = 1),
          axis.title.y = element_text(size=15),
          axis.text.x=element_text(angle = 45,hjust=0.5,size=12,color='black',vjust = 0.6),
          axis.text.y=element_text(vjust=0.5,size=12,color='black',hjust=0.5),
          axis.ticks = element_line(color = "black", size = 1),
          axis.ticks.length = unit(0.15, "cm"))
  dev.off()
  
}

