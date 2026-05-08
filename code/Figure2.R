library(minpack.lm)
# S(d) = S0 * exp(-k * d) + S_b

distance <- c(-(0- 5)/2, 5-(5- 15)/2, 15-(15- 30)/2,
              30-(30- 60)/2, 60-(60- 100)/2, 100-(100- 200)/2)

x <- distance
y <- c(210,197,178,172,164,164)

y_last <- tail(y, 1)

fit <- nlsLM(
  y ~ c + a * exp(-b * x),
  start = list(a = max(y) - min(y), b = 0.1, c = min(y)),
  lower = c(a = 0, b = 0, c = -Inf),
  upper = c(a = Inf, b = Inf, c = y_last),
  control = nls.lm.control(maxiter = 1000)
)

summary(fit)

par <- coef(fit)
l <- S_function(par,seq(0,300,1))

par(mar=c(4,4,2,1),las=1)
plot(distance,y,type = "n",xlab="",#xlab="Soil Classification",
     ylab = "",#ylab = "Soil Depth (m)",#yaxt = "n",
     lty=1,mgp=c(2.5,1,0),cex.lab=1.5,yaxs = "i",xaxs = "i",cex.axis=1.3,
     xlim = c(0,250),ylim = c(150,220),axes=FALSE)
rect(par("usr")[1],  par("usr")[3],5,  par("usr")[4], col = rgb(255,227,246,150,maxColorValue = 255), border = NA)
rect(5, par("usr")[3],15, par("usr")[4], col = rgb(234,193,255,150,maxColorValue = 255), border = NA)
rect(15, par("usr")[3],30, par("usr")[4], col = rgb(187,169,255,150,maxColorValue = 255), border = NA)
rect(30, par("usr")[3],60, par("usr")[4], col = rgb(124,149,255,150,maxColorValue = 255), border = NA)
rect(60, par("usr")[3],100, par("usr")[4], col = rgb(72,140,255,150,maxColorValue = 255), border = NA)
rect(100, par("usr")[3],200, par("usr")[4], col = rgb(0,125,189,150,maxColorValue = 255), border = NA)
box(lwd=1)
points(distance,y,cex=1.2,pch=19)
lines(distance,y,lwd=2,lty=1)
# lines(c(distance,200),c(140,160,190,190,190,190,190),lwd=2)
lines(c(0,2.5),c(y[1],y[1]),lwd=2,lty=3)
lines(c(0,10),c(y[2],y[2]),lwd=2,lty=3)
lines(c(0,22.5),c(y[3],y[3]),lwd=2,lty=3)
lines(c(0,45),c(y[4],y[4]),lwd=2,lty=3)
lines(c(0,150),c(y[5],y[5]),lwd=2,lty=3)
lines(seq(0,300,1),l,col="red",lwd=2)
axis(2,at = y,labels =y,cex.axis=1.5)
axis(1,at = c(0,50,100,150,200,250),labels =c(0,50,100,150,200,250),cex.axis=1.5)
colors <- c("black", "red")
labels <- c("Number of Modules", "Exponential Decay Fit")

legend("topright",
       legend = labels,
       col = colors,
       lty = 1,
       lwd = 2,
       bg = "#fdc38d",
       bty = "n")
dev.off()




json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer1/json/k210-1.json")
k1 <- 210
module1<-json_data$max_omega_logi
length(table(module1))

json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer2/json/k197-1.json")
k2 <- 197
module2<-json_data$max_omega_logi
length(table(module2))

json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer3/json/k178-1.json")
k3 <- 178
module3<-json_data$max_omega_logi
length(table(module3))

json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer4/json/k172-1.json")
k4 <- 172
module4<-json_data$max_omega_logi
length(table(module4))

json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer5/json/k164-1.json")
k5 <- 164
module5<-json_data$max_omega_logi
length(table(module5))

json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer6/json/k164-1.json")
k6 <- 164
module6<-json_data$max_omega_logi
length(table(module6))



jaccard_matrix12 <- outer(
  (1:k1) - 1,  
  (1:k2) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module1 == i)  
    cj <- which(module2 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)


###########

jaccard_matrix23 <- outer(
  (1:k2) - 1,  
  (1:k3) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module2 == i)  
    cj <- which(module3 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix34 <- outer(
  (1:k3) - 1,  
  (1:k4) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module3 == i)  
    cj <- which(module4 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix45 <- outer(
  (1:k4) - 1,  
  (1:k5) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module4 == i)  
    cj <- which(module5 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)


###########

jaccard_matrix56 <- outer(
  (1:k5) - 1,  
  (1:k6) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module5 == i)  
    cj <- which(module6 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)


###########

jaccard_matrix13 <- outer(
  (1:k1) - 1,  
  (1:k3) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module1 == i)  
    cj <- which(module3 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix14 <- outer(
  (1:k1) - 1,  
  (1:k4) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module1 == i)  
    cj <- which(module4 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))      
    insert_num / union_num                   
  })
)

###########

jaccard_matrix15 <- outer(
  (1:k1) - 1,  
  (1:k5) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module1 == i)  
    cj <- which(module5 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix16 <- outer(
  (1:k1) - 1,  
  (1:k6) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module1 == i)  
    cj <- which(module6 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix24 <- outer(
  (1:k2) - 1,  
  (1:k4) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module2 == i)  
    cj <- which(module4 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                  
  })
)

###########

jaccard_matrix25 <- outer(
  (1:k2) - 1,  
  (1:k5) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module2 == i)  
    cj <- which(module5 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix26 <- outer(
  (1:k2) - 1,  
  (1:k6) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module2 == i)  
    cj <- which(module6 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix35 <- outer(
  (1:k3) - 1,  
  (1:k5) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module3 == i)  
    cj <- which(module5 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix36 <- outer(
  (1:k3) - 1,  
  (1:k6) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module3 == i)  
    cj <- which(module6 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)

###########

jaccard_matrix46 <- outer(
  (1:k4) - 1,  
  (1:k6) - 1,       
  Vectorize(function(i, j) {
    ci <- which(module4 == i)  
    cj <- which(module6 == j)  
    insert_num <- length(intersect(ci, cj))  
    union_num <- length(union(ci, cj))       
    insert_num / union_num                   
  })
)


library(mclust)
library(aricode)

calculate_vertical_connectivity <- function(jaccard_mat, labels_1, labels_2) {
  
  max_j_cols <- apply(jaccard_mat, 2, max)
  colnames(jaccard_mat) <- (1:ncol(jaccard_mat))-1
  size_2 <- table(labels_2)[colnames(jaccard_mat)]
  
  score_2_to_1 <- sum(max_j_cols * size_2) / sum(size_2)
  
  return(score_2_to_1)
}


#####################
j12 <- calculate_vertical_connectivity(jaccard_matrix12, module1, module2)
ari_12 <- adjustedRandIndex(as.character(module1), as.character(module2))
nmi_12 <- NMI(as.character(module1), as.character(module2))

j23 <- calculate_vertical_connectivity(jaccard_matrix23, module2, module3)
ari_23 <- adjustedRandIndex(as.character(module2), as.character(module3))
nmi_23 <- NMI(as.character(module2), as.character(module3))

j34 <- calculate_vertical_connectivity(jaccard_matrix34, module3, module4)
ari_34 <- adjustedRandIndex(as.character(module3), as.character(module4))
nmi_34 <- NMI(as.character(module3), as.character(module4))

j45 <- calculate_vertical_connectivity(jaccard_matrix45, module4, module5)
ari_45 <- adjustedRandIndex(as.character(module4), as.character(module5))
nmi_45 <- NMI(as.character(module4), as.character(module5))

j56 <- calculate_vertical_connectivity(jaccard_matrix56, module5, module6)
ari_56 <- adjustedRandIndex(as.character(module5), as.character(module6))
nmi_56 <- NMI(as.character(module5), as.character(module6))

########################
j13 <- calculate_vertical_connectivity(jaccard_matrix13, module1, module3)
ari_13 <- adjustedRandIndex(as.character(module1), as.character(module3))
nmi_13 <- NMI(as.character(module1), as.character(module3))

j24 <- calculate_vertical_connectivity(jaccard_matrix24, module2, module4)
ari_24 <- adjustedRandIndex(as.character(module2), as.character(module4))
nmi_24 <- NMI(as.character(module2), as.character(module4))

j35 <- calculate_vertical_connectivity(jaccard_matrix35, module3, module5)
ari_35 <- adjustedRandIndex(as.character(module3), as.character(module5))
nmi_35 <- NMI(as.character(module3), as.character(module5))

j46 <- calculate_vertical_connectivity(jaccard_matrix46, module4, module6)
ari_46 <- adjustedRandIndex(as.character(module4), as.character(module6))
nmi_46 <- NMI(as.character(module4), as.character(module6))

########################
j14 <- calculate_vertical_connectivity(jaccard_matrix14, module1, module4)
ari_14 <- adjustedRandIndex(as.character(module1), as.character(module4))
nmi_14 <- NMI(as.character(module1), as.character(module4))

j25 <- calculate_vertical_connectivity(jaccard_matrix25, module2, module5)
ari_25 <- adjustedRandIndex(as.character(module2), as.character(module5))
nmi_25 <- NMI(as.character(module2), as.character(module5))

j36 <- calculate_vertical_connectivity(jaccard_matrix36, module3, module6)
ari_36 <- adjustedRandIndex(as.character(module3), as.character(module6))
nmi_36 <- NMI(as.character(module3), as.character(module6))

########################
j15 <- calculate_vertical_connectivity(jaccard_matrix15, module1, module5)
ari_15 <- adjustedRandIndex(as.character(module1), as.character(module5))
nmi_15 <- NMI(as.character(module1), as.character(module5))

j26 <- calculate_vertical_connectivity(jaccard_matrix26, module2, module6)
ari_26 <- adjustedRandIndex(as.character(module2), as.character(module6))
nmi_26 <- NMI(as.character(module2), as.character(module6))

########################
j16 <- calculate_vertical_connectivity(jaccard_matrix16, module1, module6)
ari_16 <- adjustedRandIndex(as.character(module1), as.character(module6))
nmi_16 <- NMI(as.character(module1), as.character(module6))




data <- data.frame(
  category = c(rep("IB2L",5),rep("IB3L",4),rep("IB4L",3),rep("IB5L",2),"IB6L"),
  name = c("1-2","2-3","3-4","4-5","5-6",
           "1-3","2-4","3-5","4-6",
           "1-4","2-5","3-6",
           "1-5","2-6",
           "1-6"),
  value = c(j12, j23,j34,j45,j56,j13,j24,j35,j46,j14,j25,j36,j15,j26,j16)
)
data$name <- factor(data$name, levels = c("1-2","2-3","3-4","4-5","5-6",
                                          "1-3","2-4","3-5","4-6",
                                          "1-4","2-5","3-6",
                                          "1-5","2-6",
                                          "1-6"))
data$category <- factor(data$category, levels = c("IB2L", "IB3L", "IB4L", "IB5L", "IB6L"))
data$xpos <- c(seq(1,5,1),seq(1,4,1)+5.5,seq(1,3,1)+10,seq(1,2,1)+13.5,17)
data$color <- c("#e18cb5","#e18cb5","#e18cb5","#e18cb5","#e18cb5",
                "#6bc2af","#6bc2af","#6bc2af","#6bc2af",
                "#e6ae7f","#e6ae7f","#e6ae7f",
                "#a9a3ca","#a9a3ca",
                "#89aad7")

ggplot(data, aes(x = xpos, y = value,fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = data$color) +
  scale_x_continuous(breaks = data$xpos, labels = data$name,expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0,0.31),expand = expansion(mult = c(0.01, 0.1)))+
  geom_text(aes(label = sprintf("%0.3f", value)), vjust = -0.3) + 
  labs(x = "", y = "", title = "") +
  theme_minimal()+
  theme(panel.background = element_rect(linewidth=1, color = 'black', fill = 'transparent'),
        panel.grid = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.title.x = element_text(size=15,vjust = 1),
        axis.title.y = element_text(size=15),
        axis.text.x=element_text(hjust=0.5,size=12,color='black',vjust = 0.6),
        axis.text.y=element_text(vjust=0.5,size=12,color='black',hjust=0.5),
        axis.ticks = element_line(color = "black", size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none")
dev.off()





data <- data.frame(
  category = c(rep("IB2L",5),rep("IB3L",4),rep("IB4L",3),rep("IB5L",2),"IB6L"),
  name = c("1-2","2-3","3-4","4-5","5-6",
           "1-3","2-4","3-5","4-6",
           "1-4","2-5","3-6",
           "1-5","2-6",
           "1-6"),
  value = c(ari_12, ari_23,ari_34,ari_45,ari_56,ari_13,ari_24,ari_35,ari_46,
            ari_14,ari_25,ari_36,ari_15,ari_26,ari_16)
)
data$name <- factor(data$name, levels = c("1-2","2-3","3-4","4-5","5-6",
                                          "1-3","2-4","3-5","4-6",
                                          "1-4","2-5","3-6",
                                          "1-5","2-6",
                                          "1-6"))
data$category <- factor(data$category, levels = c("IB2L", "IB3L", "IB4L", "IB5L", "IB6L"))
data$xpos <- c(seq(1,5,1),seq(1,4,1)+5.5,seq(1,3,1)+10,seq(1,2,1)+13.5,17)
data$color <- c("#e18cb5","#e18cb5","#e18cb5","#e18cb5","#e18cb5",
                "#6bc2af","#6bc2af","#6bc2af","#6bc2af",
                "#e6ae7f","#e6ae7f","#e6ae7f",
                "#a9a3ca","#a9a3ca",
                "#89aad7")

ggplot(data, aes(x = xpos, y = value,fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = data$color) +
  scale_x_continuous(breaks = data$xpos, labels = data$name,expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0,0.35),expand = expansion(mult = c(0.01, 0.1)))+
  geom_text(aes(label = sprintf("%0.3f", value)), vjust = -0.3) + 
  labs(x = "", y = "", title = "") +
  theme_minimal()+
  theme(panel.background = element_rect(linewidth=1, color = 'black', fill = 'transparent'),
        panel.grid = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.title.x = element_text(size=15,vjust = 1),
        axis.title.y = element_text(size=15),
        axis.text.x=element_text(hjust=0.5,size=12,color='black',vjust = 0.6),
        axis.text.y=element_text(vjust=0.5,size=12,color='black',hjust=0.5),
        axis.ticks = element_line(color = "black", size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none")
dev.off()






data <- data.frame(
  category = c(rep("IB2L",5),rep("IB3L",4),rep("IB4L",3),rep("IB5L",2),"IB6L"),
  name = c("1-2","2-3","3-4","4-5","5-6",
           "1-3","2-4","3-5","4-6",
           "1-4","2-5","3-6",
           "1-5","2-6",
           "1-6"),
  value = c(nmi_12, nmi_23,nmi_34,nmi_45,nmi_56,nmi_13,nmi_24,nmi_35,nmi_46,
            nmi_14,nmi_25,nmi_36,nmi_15,nmi_26,nmi_16)
)
data$name <- factor(data$name, levels = c("1-2","2-3","3-4","4-5","5-6",
                                          "1-3","2-4","3-5","4-6",
                                          "1-4","2-5","3-6",
                                          "1-5","2-6",
                                          "1-6"))
data$category <- factor(data$category, levels = c("IB2L", "IB3L", "IB4L", "IB5L", "IB6L"))
data$xpos <- c(seq(1,5,1),seq(1,4,1)+5.5,seq(1,3,1)+10,seq(1,2,1)+13.5,17)
data$color <- c("#e18cb5","#e18cb5","#e18cb5","#e18cb5","#e18cb5",
                "#6bc2af","#6bc2af","#6bc2af","#6bc2af",
                "#e6ae7f","#e6ae7f","#e6ae7f",
                "#a9a3ca","#a9a3ca",
                "#89aad7")

ggplot(data, aes(x = xpos, y = value,fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = data$color) +
  scale_x_continuous(breaks = data$xpos, labels = data$name,expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0,0.76),expand = expansion(mult = c(0.01, 0.1)))+
  geom_text(aes(label = sprintf("%0.3f", value)), vjust = -0.3) + 
  labs(x = "", y = "", title = "") +
  theme_minimal()+
  theme(panel.background = element_rect(linewidth=1, color = 'black', fill = 'transparent'),
        panel.grid = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.title.x = element_text(size=15,vjust = 1),
        axis.title.y = element_text(size=15),
        axis.text.x=element_text(hjust=0.5,size=12,color='black',vjust = 0.6),
        axis.text.y=element_text(vjust=0.5,size=12,color='black',hjust=0.5),
        axis.ticks = element_line(color = "black", size = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none")
dev.off()




soil_table12 <- table(as.character(module1), as.character(module2))
chi_test <- chisq.test(soil_table12)
print(chi_test)
cramer_v12 <- CramerV(soil_table12)

soil_table23 <- table(as.character(module2), as.character(module3))
chi_test <- chisq.test(soil_table23)
print(chi_test)
cramer_v23 <- CramerV(soil_table23)

soil_table34 <- table(as.character(module3), as.character(module4))
chi_test <- chisq.test(soil_table34)
print(chi_test)
cramer_v34 <- CramerV(soil_table34)

soil_table45 <- table(as.character(module4), as.character(module5))
chi_test <- chisq.test(soil_table45)
print(chi_test)
cramer_v45 <- CramerV(soil_table45)

soil_table56 <- table(as.character(module5), as.character(module6))
chi_test <- chisq.test(soil_table56)
print(chi_test)
cramer_v56 <- CramerV(soil_table56)



soil_table13 <- table(as.character(module1), as.character(module3))
chi_test <- chisq.test(soil_table13)
print(chi_test)
cramer_v13 <- CramerV(soil_table13)

soil_table24 <- table(as.character(module2), as.character(module4))
chi_test <- chisq.test(soil_table24)
print(chi_test)
cramer_v24 <- CramerV(soil_table24)

soil_table35 <- table(as.character(module3), as.character(module5))
chi_test <- chisq.test(soil_table35)
print(chi_test)
cramer_v35 <- CramerV(soil_table35)

soil_table46 <- table(as.character(module4), as.character(module6))
chi_test <- chisq.test(soil_table46)
print(chi_test)
cramer_v46 <- CramerV(soil_table46)



soil_table14 <- table(as.character(module1), as.character(module4))
chi_test <- chisq.test(soil_table14)
print(chi_test)
cramer_v14 <- CramerV(soil_table14)

soil_table25 <- table(as.character(module2), as.character(module5))
chi_test <- chisq.test(soil_table25)
print(chi_test)
cramer_v25 <- CramerV(soil_table25)

soil_table36 <- table(as.character(module3), as.character(module6))
chi_test <- chisq.test(soil_table36)
print(chi_test)
cramer_v36 <- CramerV(soil_table36)



soil_table15 <- table(as.character(module1), as.character(module5))
chi_test <- chisq.test(soil_table15)
print(chi_test)
cramer_v15 <- CramerV(soil_table15)

soil_table26 <- table(as.character(module2), as.character(module6))
chi_test <- chisq.test(soil_table26)
print(chi_test)
cramer_v26 <- CramerV(soil_table26)



soil_table16 <- table(as.character(module1), as.character(module6))
chi_test <- chisq.test(soil_table16)
print(chi_test)
cramer_v16 <- CramerV(soil_table16)


par(mfrow=c(1,5),las=1)
par(mar = c(0, 0, 0,0),oma = c(5.5, 6, 3.5, 1),mgp=c(1,1,0))

plot(1,1,type="n",ylim=c(0.25,0.53),xlim=c(0.5,6.5),ylab="",xlab="",axes=FALSE)
polygon(c(1,2,3,2),c(rep(cramer_v12,2),rep(cramer_v23,2)), 
        col = rgb(251,208,213,150,maxColorValue = 255), 
        border = rgb(251,208,213,150,maxColorValue = 255))
polygon(c(2,3,4,3),c(rep(cramer_v23,2),rep(cramer_v34,2)), 
        col = rgb(251,208,213,150,maxColorValue = 255), 
        border = rgb(251,208,213,150,maxColorValue = 255))
polygon(c(3,4,5,4),c(rep(cramer_v34,2),rep(cramer_v45,2)), 
        col = rgb(251,208,213,150,maxColorValue = 255), 
        border = rgb(251,208,213,150,maxColorValue = 255))
polygon(c(4,5,6,5),c(rep(cramer_v45,2),rep(cramer_v56,2)), 
        col = rgb(251,208,213,150,maxColorValue = 255), 
        border = rgb(251,208,213,150,maxColorValue = 255))
lines(c(1,2),rep(cramer_v12,2),lwd=2,col="#e18cb5")
lines(c(2,3),rep(cramer_v23,2),lwd=2,col="#e18cb5")
lines(c(3,4),rep(cramer_v34,2),lwd=2,col="#e18cb5")
lines(c(4,5),rep(cramer_v45,2),lwd=2,col="#e18cb5")
lines(c(5,6),rep(cramer_v56,2),lwd=2,col="#e18cb5")
box(lwd=1)
axis(1,at=1:6,labels =1:6,cex.axis=1.5)
axis(2,at = seq(0.2,0.5,length=5),
     labels =sprintf("%0.2f", seq(0.2,0.5,length=5)),
     cex.axis=1.5)
mtext("SB2L",side = 3,line = 0.3,cex = 1.2)


plot(1,1,type="n",ylim=c(0.25,0.53),xlim=c(0.5,6.5),ylab="",xlab="",axes=FALSE)
polygon(c(1,3,4,2),c(rep(cramer_v13,2),rep(cramer_v24,2)), 
        col = rgb(210,231,167,200,maxColorValue = 255), 
        border = rgb(210,231,167,200,maxColorValue = 255))
polygon(c(2,4,5,3),c(rep(cramer_v24,2),rep(cramer_v35,2)), 
        col = rgb(210,231,167,200,maxColorValue = 255), 
        border = rgb(210,231,167,200,maxColorValue = 255))
polygon(c(3,5,6,4),c(rep(cramer_v35,2),rep(cramer_v46,2)), 
        col = rgb(210,231,167,200,maxColorValue = 255), 
        border = rgb(210,231,167,200,maxColorValue = 255))
lines(c(1,3),rep(cramer_v13,2),lwd=2,col="#77ae43")
lines(c(2,4),rep(cramer_v24,2),lwd=2,col="#77ae43")
lines(c(3,5),rep(cramer_v35,2),lwd=2,col="#77ae43")
lines(c(4,6),rep(cramer_v46,2),lwd=2,col="#77ae43")
box(lwd=1)
axis(1,at=1:6,labels =1:6,cex.axis=1.5)
mtext("SB3L",side = 3,line = 0.3,cex = 1.2)


plot(1,1,type="n",ylim=c(0.25,0.53),xlim=c(0.5,6.5),ylab="",xlab="",axes=FALSE)
polygon(c(1,4,5,2),c(rep(cramer_v14,2),rep(cramer_v25,2)), 
        col = rgb(250,216,183,180,maxColorValue = 255), 
        border = rgb(250,216,183,180,maxColorValue = 255))
polygon(c(2,5,6,3),c(rep(cramer_v25,2),rep(cramer_v36,2)), 
        col = rgb(250,216,183,180,maxColorValue = 255), 
        border = rgb(250,216,183,180,maxColorValue = 255))
lines(c(1,4),rep(cramer_v14,2),lwd=2,col="#edb021")
lines(c(2,5),rep(cramer_v25,2),lwd=2,col="#edb021")
lines(c(3,6),rep(cramer_v36,2),lwd=2,col="#edb021")
box(lwd=1)
axis(1,at=1:6,labels =1:6,cex.axis=1.5)
mtext("SB4L",side = 3,line = 0.3,cex = 1.2)

plot(1,1,type="n",ylim=c(0.25,0.53),xlim=c(0.5,6.5),ylab="",xlab="",axes=FALSE)
polygon(c(1,5,6,2),c(rep(cramer_v15,2),rep(cramer_v26,2)), 
        col = rgb(221,196,237,220,maxColorValue = 255), 
        border = rgb(221,196,237,220,maxColorValue = 255))
lines(c(1,5),rep(cramer_v15,2),lwd=2,col="#7f318d")
lines(c(2,6),rep(cramer_v26,2),lwd=2,col="#7f318d")
box(lwd=1)
axis(1,at=1:6,labels =1:6,cex.axis=1.5)
mtext("SB5L",side = 3,line = 0.3,cex = 1.2)

plot(1,1,type="n",ylim=c(0.25,0.53),xlim=c(1,6),ylab="",xlab="",axes=FALSE)
lines(c(1,6),rep(cramer_v16,2),lwd=2,col="#1072BD")
box(lwd=1)
axis(1,at=1:6,labels =1:6,cex.axis=1.5)
mtext("SB6L",side = 3,line = 0.3,cex = 1.2)

mtext("Soil Layer",side = 1,line = 3.5,cex=1.4,outer = TRUE)
mtext("Cramer's V",side = 2,line = 4,cex=1.4,outer = TRUE,las=0)
dev.off()
