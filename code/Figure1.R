
datt <- rbind(depth1,depth2,depth3,depth4,depth5,depth6)
select_datt<- datt[,c("BD","CEC","SOC","TK","TN","TP")]

scale_datt <- apply(select_datt, 2, scale)
scale_datt[,1] <- -1*scale_datt[,1]

emf1 <- rowMeans(scale_datt)
matrix_emf1 <- matrix(emf1, nrow = 279120)

emf2 <- rowMeans(matrix_emf1)


n_levels <- 1000
my_palette <- colorRampPalette(c("#223e36","#235D3A","#397D54","#468c37","#96C35E","#B4D780","#E3F7A0","#F0F7D3","#FFE388","#FAC463","#FE993B","#CA752F","#9F5221","#91341D"))(n_levels)


color_indices <- as.numeric(cut(emf2, breaks = n_levels))
final_colors <- my_palette[color_indices]

library(fields)

plot(1, 1, col = final_colors, pch = 16, cex = 0.5, main = "",axes=FALSE,type="n",xlab="")

image.plot(legend.only = TRUE, 
           zlim = range(emf2, na.rm = TRUE), 
           col = my_palette,
           horizontal = TRUE,  
           legend.args = list(text = "EMF Value", side = 1, line = 2)) 



par(mar=c(0.1,1,0.1,1))

x <- as.numeric(names(table(n1)))
x <- x[seq(1, length(x), by = 10)]

y <- as.numeric(names(table(n2)))
y <- y[seq(1, length(y), by = 10)]


nii <- c()
for (i in 1:length(x)) {
  for(j in 1:length(y)){
    nii <- c(nii,which(n1==x[i]&n2==y[j]))
  }
}

plot(n1[nii],n2[nii], col = "gray90", pch = 15,ylab = "",xlab = "",axes=FALSE,#type="n",
     ylim=c(31.6,45),cex=0.7)
points(n1[nii],n2[nii], col = final_colors[nii], pch = 15,cex=0.7)
points(coords[,1],coords[,2],col="blue",cex=0.1,type = "p")

dev.off()



sample_units <- c(18909,81608,103132)
par(mfrow=c(3,6),las=1)
par(mar = c(0, 0, 0,0),oma = c(3.1, 3.2, 2, 2),mgp=c(1,0.5,0))

for (i in 1:3) {
  sample_i <- sample_units[i]
  
  y1 <- range(fit1$original_data,fit2$original_data,fit3$original_data,
              fit4$original_data,fit5$original_data,fit6$original_data)
  lim2 <- y1[2]-y1[1]
  
  x1 <- range(fit1$Time,fit2$Time,fit3$Time,fit4$Time,fit5$Time,fit6$Time)
  lim1 <- x1[2]-x1[1]
  
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.05*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  x<-par("usr")  
  rect(xleft=x[1],ybottom=x[3],xright=x[2],ytop=x[4],col=rgb(230,230,250,150,maxColorValue = 255))
  # abline(h=max(fit1$power_fit[,sample_i]),lty=2)
  # abline(h=max(fit2$power_fit[,sample_i]),col="#cfcfe5",lty=2)
  # abline(h=max(fit3$power_fit[,sample_i]),col="#fad9ea",lty=2)
  # abline(h=max(fit4$power_fit[,sample_i]),col="#f3eeaf",lty=2)
  # abline(h=max(fit5$power_fit[,sample_i]),col="#b8e3b2",lty=2)
  # abline(h=max(fit6$power_fit[,sample_i]),col="#fdc38d",lty=2)
  points(fit1$Time,fit1$original_data[,sample_i],col="steelblue1",pch=16,cex=1.2)
  par <- power_equation_base2(fit1$Time, fit1$original_data[,sample_i])
  lines(fit1$Time,power.equation(par,fit1$Time),col="#3568c2",lwd=2)
  # lines(fit1$Time,fit1$power_fit[,sample_i],col="#3568c2",lwd=2)
  if(i==1){
    mtext("Layer 1",side = 3,line = 0.3,cex = 0.8)
  }
  box(lwd=1)
  if(i==3){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1,gap.axis = -1)
  }
  axis(2,at = seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4),
       labels =sprintf("%0.1f", seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4)),cex.axis=1)
  
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.05*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  x<-par("usr")  
  rect(xleft=x[1],ybottom=x[3],xright=x[2],ytop=x[4],col=rgb(255,239,213,150,maxColorValue = 255))
  # abline(h=max(fit1$power_fit[,sample_i]),col="#C3DFED",lty=2)
  # abline(h=max(fit2$power_fit[,sample_i]),lty=2)
  # abline(h=max(fit3$power_fit[,sample_i]),col="#fad9ea",lty=2)
  # abline(h=max(fit4$power_fit[,sample_i]),col="#f3eeaf",lty=2)
  # abline(h=max(fit5$power_fit[,sample_i]),col="#b8e3b2",lty=2)
  # abline(h=max(fit6$power_fit[,sample_i]),col="#fdc38d",lty=2)
  points(fit2$Time,fit2$original_data[,sample_i],col="steelblue1",pch=16,cex=1.2)
  par <- power_equation_base2(fit2$Time, fit2$original_data[,sample_i])
  lines(fit2$Time,power.equation(par,fit2$Time),col="#3568c2",lwd=2)
  # lines(fit2$Time,fit2$power_fit[,sample_i],col="#3568c2",lwd=2)
  if(i==1){
    mtext("Layer 2",side = 3,line = 0.3,cex = 0.8)
  }
  box(lwd=1)
  if(i==3){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1,gap.axis = -1)
  }
  
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.05*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  x<-par("usr")  
  rect(xleft=x[1],ybottom=x[3],xright=x[2],ytop=x[4],col=rgb(230,230,250,150,maxColorValue = 255))
  # abline(h=max(fit1$power_fit[,sample_i]),col="#C3DFED",lty=2)
  # abline(h=max(fit2$power_fit[,sample_i]),col="#cfcfe5",lty=2)
  # abline(h=max(fit3$power_fit[,sample_i]),lty=2)
  # abline(h=max(fit4$power_fit[,sample_i]),col="#f3eeaf",lty=2)
  # abline(h=max(fit5$power_fit[,sample_i]),col="#b8e3b2",lty=2)
  # abline(h=max(fit6$power_fit[,sample_i]),col="#fdc38d",lty=2)
  points(fit3$Time,fit3$original_data[,sample_i],col="steelblue1",pch=16,cex=1.2)
  par <- power_equation_base2(fit3$Time, fit3$original_data[,sample_i])
  lines(fit3$Time,power.equation(par,fit3$Time),col="#3568c2",lwd=2)
  if(i==1){
    mtext("Layer 3",side = 3,line = 0.3,cex = 0.8)
  }
  box(lwd=1)
  if(i==3){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1,gap.axis = -1)
  }
  
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.05*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  x<-par("usr")  
  rect(xleft=x[1],ybottom=x[3],xright=x[2],ytop=x[4],col=rgb(255,239,213,150,maxColorValue = 255))
  # abline(h=max(fit1$power_fit[,sample_i]),col="#C3DFED",lty=2)
  # abline(h=max(fit2$power_fit[,sample_i]),col="#cfcfe5",lty=2)
  # abline(h=max(fit3$power_fit[,sample_i]),col="#fad9ea",lty=2)
  # abline(h=max(fit4$power_fit[,sample_i]),lty=2)
  # abline(h=max(fit5$power_fit[,sample_i]),col="#b8e3b2",lty=2)
  # abline(h=max(fit6$power_fit[,sample_i]),col="#fdc38d",lty=2)
  points(fit4$Time,fit4$original_data[,sample_i],col="steelblue1",pch=16,cex=1.2)
  par <- power_equation_base2(fit4$Time, fit4$original_data[,sample_i])
  lines(fit4$Time,power.equation(par,fit4$Time),col="#3568c2",lwd=2)
  # lines(fit4$Time,fit4$power_fit[,sample_i],col="#3568c2",lwd=2)
  if(i==1){
    mtext("Layer 4",side = 3,line = 0.3,cex = 0.8)
  }
  box(lwd=1)
  if(i==3){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1,gap.axis = -1)
  }
  
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.05*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  x<-par("usr") 
  rect(xleft=x[1],ybottom=x[3],xright=x[2],ytop=x[4],col=rgb(230,230,250,150,maxColorValue = 255))
  # abline(h=max(fit1$power_fit[,sample_i]),col="#C3DFED",lty=2)
  # abline(h=max(fit2$power_fit[,sample_i]),col="#cfcfe5",lty=2)
  # abline(h=max(fit3$power_fit[,sample_i]),col="#fad9ea",lty=2)
  # abline(h=max(fit4$power_fit[,sample_i]),col="#f3eeaf",lty=2)
  # abline(h=max(fit5$power_fit[,sample_i]),lty=2)
  # abline(h=max(fit6$power_fit[,sample_i]),col="#fdc38d",lty=2)
  points(fit5$Time,fit5$original_data[,sample_i],col="steelblue1",pch=16,cex=1.2)
  par <- power_equation_base2(fit5$Time, fit5$original_data[,sample_i])
  lines(fit5$Time,power.equation(par,fit5$Time),col="#3568c2",lwd=2)
  # lines(fit5$Time,fit5$power_fit[,sample_i],col="#3568c2",lwd=2)
  if(i==1){
    mtext("Layer 5",side = 3,line = 0.3,cex = 0.8)
  }
  box(lwd=1)
  if(i==3){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1,gap.axis = -1)
  }
  
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.05*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  x<-par("usr")  
  rect(xleft=x[1],ybottom=x[3],xright=x[2],ytop=x[4],col=rgb(255,239,213,150,maxColorValue = 255))
  # abline(h=max(fit1$power_fit[,sample_i]),col="#C3DFED",lty=2)
  # abline(h=max(fit2$power_fit[,sample_i]),col="#cfcfe5",lty=2)
  # abline(h=max(fit3$power_fit[,sample_i]),col="#fad9ea",lty=2)
  # abline(h=max(fit4$power_fit[,sample_i]),col="#f3eeaf",lty=2)
  # abline(h=max(fit5$power_fit[,sample_i]),col="#b8e3b2",lty=2)
  # abline(h=max(fit6$power_fit[,sample_i]),lty=2)
  points(fit6$Time,fit6$original_data[,sample_i],col="steelblue1",pch=16,cex=1.2)
  par <- power_equation_base2(fit6$Time, fit6$original_data[,sample_i])
  lines(fit6$Time,power.equation(par,fit6$Time),col="#3568c2",lwd=2)
  # lines(fit6$Time,fit6$power_fit[,sample_i],col="#3568c2",lwd=2)
  if(i==1){
    mtext("Layer 6",side = 3,line = 0.3,cex = 0.8)
  }
  box(lwd=1)
  if(i==3){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1,gap.axis = -1)
  }
  mi <- paste("grid cell",sample_i,sep = " ")
  mtext(mi,side = 4,line = 0.5,las=-0.5,cex = 0.8)
  
}

mtext("Compartment Index (×10000)",side = 1,line = 2,cex=1,outer = TRUE)
mtext("Value of Soil Properties",side = 2,line = 2,cex=1,outer = TRUE,las=0)


dev.off()
