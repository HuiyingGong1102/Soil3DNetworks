
datt <- rbind(depth1,depth2,depth3,depth4,depth5,depth6)
Normalized_depth <- t(sweep(datt, 2, apply(datt, 2, max), FUN = "/"))

data_all <- rbind(Normalized_depth[,1:279120],Normalized_depth[,(1+279120):(2*279120)],
                  Normalized_depth[,(1+2*279120):(3*279120)],Normalized_depth[,(1+3*279120):(4*279120)],
                  Normalized_depth[,(1+4*279120):(5*279120)],Normalized_depth[,(1+5*279120):(6*279120)])
dim(data_all)

library(jsonlite)
json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer_all/json/k434-1.json")
k <- 434
module1<-json_data$max_omega_logi
length(table(module1))
table(module1)

data_module <-rep()
for (i in 1:k)
{
  index<-which(module1==names(table(module1))[i])
  if(length(index)==1){
    new_i <- data_all[,index]
  }else{
    new_i<-apply(data_all[,index],1,mean)
  }
  data_module<-rbind(data_module,new_i)
}
data_module<-t(data_module)
colnames(data_module) <- paste("M",1:k,sep = "")

total <- apply(data_module,1,sum)
index <- total[order(total)]
data_order <- data_module[order(total),]

fit <- power_equation_fit2(
  x = index,
  y = data_order,
  X_smooth = seq(index[1], index[length(index)], length.out = 90),
  thread = 24)

log_data_order <- log(data_order)
log_index <- log(index)

norder <- c(2:50)

bic<-c()
for (i in norder) {
  filename<-paste("D:/Soil3DNetworks/code/cluster/Layer_all_module/json/k",i,"-",1,".json",sep = "")
  json_data <- fromJSON(filename)
  BIC <- json_data$BIC
  bic <- c(bic,BIC)
}


filename <- "D:/Soil3DNetworks/code/cluster/Layer_all_module/json/k33-1.json"
k=33
json_data <- fromJSON(filename)
module<-json_data$max_omega_logi
# # module<-json_data$clustered
table(module)
length(table(module))
range(table(module))
k <- norder[which(bic==min(bic))]

plot(norder,bic,type = "p",xlab="",ylab = "",axes=FALSE,
     cex=1.5,col="darkgray",mgp=c(2.5,1,0),cex.lab=1.5)
lines(norder,bic,col="red",lwd=2)
abline(v=k,lty=2,col="blue",lwd=2)
box(lwd=1)
axis(2,at=seq(min(bic),max(bic),length=5),
     labels = sprintf("%0.0f", seq(min(bic),max(bic),length=5)))
axis(1,at=seq(2,20,3),labels = seq(2,20,3))
mtext("Modules", side=1, line=3,cex = 1.5)
mtext("BIC", side=2, line=3,cex = 1.5)


n_levels <- 1000
# my_palette <- colorRampPalette(c("#96C35E","#B4D780","#E3F7A0","#F0F7D3","#FFE388","#FAC463","#FE993B","#CA752F"))(n_levels)
my_palette <- colorRampPalette(c("#468c37","#96C35E","#B4D780","#E3F7A0","#F0F7D3","#FFE388","#FAC463","#FE993B","#CA752F"))(n_levels)


# 2. 将 emf2 的数值“分箱”映射到 1 到 1000 的索引上
# cut() 函数会按照数值大小，把 emf2 均匀切分成 1000 份
color_indices <- as.numeric(cut(emf2, breaks = n_levels))

# 3. 根据索引，为这 279120 个值分配颜色
final_colors <- my_palette[color_indices]

library(fields)

# 1. 正常画您的主图
plot(1, 1, col = final_colors, pch = 16, cex = 0.5, main = "",axes=FALSE,type="n",xlab="")

# 2. 生成横向彩条
image.plot(legend.only = TRUE, 
           zlim = range(emf2, na.rm = TRUE), 
           col = my_palette,
           # horizontal = TRUE,  # <--- 核心魔法：让彩条变横向
           legend.args = list(text = "EMF Value", side = 1, line = 2)) 


plot(n1,n2, col = final_colors, pch = 16,ylab = "",xlab = "",axes=FALSE, ylim=c(31.6,45),cex=0.3)


par(mfrow=c(7,10),las=1)
par(oma = c(4, 4.5, 0.5, 0.5),mgp=c(1,0.5,0))

for (i in 1:k) {
  index<-which(module==(((1:k))-1)[i])
  
  y1 <- range(fit$power_fit)
  lim2 <- y1[2]-y1[1]
  
  x1 <- range(fit$Time)
  lim1 <- x1[2]-x1[1]
  
  par(mar = c(0, 0.8, 0,0))
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.01*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  
  for (ii in index){
    # points(fit$Time,fit$original_data[,ii],col="#7BABDD",pch=16)
    points(fit$Time,fit$original_data[,ii],
           col=rgb(213,225,227,180,maxColorValue = 255),pch=16,cex=0.6)
    # lines(fit$Time,fit$power_fit[,ii],
    #        col=rgb(195,223,237,100,maxColorValue = 255),lwd=2)
  }
  
  if(length(index)==1){
    par <- power_equation_base2(fit$Time, fit$original_data[,index])
    lines(fit$Time,power.equation(par,fit$Time),col="#8fb4be",lwd=2)
  }else{
    par <- power_equation_base2(fit$Time, rowMeans(fit$original_data[,index]))
    lines(fit$Time,power.equation(par,fit$Time),col="#8fb4be",lwd=2)
  }
  box(lwd=1)
  if(i%in%((k-5+1):k)){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)),
         cex.axis=1.1,gap.axis = -1)
  }
  
  if(i%in%seq(1,k,5)){
    axis(2,at = seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4),
         labels =sprintf("%0.1f", seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4)),
         cex.axis=1.1)
  }else{
    axis(2,at = seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4),
         labels =FALSE,
         cex.axis=1.1)
  }
  
  mi <- paste("R",i," (",length(index)," modules)",sep = "")
  text(x1[1]+0.5*lim1,y1[2]-0.03*lim2,mi,cex=1.2)
  
  par(mar = c(0, 0, 0,0.8))
  plot(n1,n2, col = "gray90", pch = 16,ylab = "",xlab = "",axes=FALSE,#type="n",
       ylim=c(31.6,45),cex=0.3)
  # points(n1[indexx],n2[indexx], col = "red", pch = 16,cex=0.5)
  # indexx <- c()
  for (ii in 1:length(index)) {
    # indexx <- c(indexx,which(module1==(ii-1)))
    indexx <- which(module1==(index[ii]-1))
    # points(n1[indexx],n2[indexx], col = "#E31A1C", pch = 16,cex=0.5)
    
    points(n1[indexx],n2[indexx], col = final_colors[indexx], pch = 16,cex=0.5)
  }
  points(coords[,1],coords[,2],col="blue",cex=0.1,type = "p")
  box(lwd=1,lty=2)
  
  # text(105.9958,44,mi,cex=1.3)
  cat("完成 plot R",i,"\n")
  
}

mtext("Compartment Index ",side = 1,line = 2.8,cex=1.4,outer = TRUE)
mtext("Value of Soil Properties",side = 2,line = 2.5,cex=1.4,outer = TRUE,las=0)





