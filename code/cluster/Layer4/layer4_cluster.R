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

norder <- c(seq(120,170,10),171:177)
bic<-c()
for (i in norder) {
  filename<-paste("D:/Soil3DNetworks/code/cluster/Layer4/json/k",i,"-",1,".json",sep = "")
  json_data <- fromJSON(filename)
  BIC <- json_data$BIC
  bic <- c(bic,BIC)
}

i=norder[which(bic==min(bic))]#172
filename<-paste("D:/Soil3DNetworks/code/cluster/Layer4/json/k",i,"-",1,".json",sep = "")
json_data <- fromJSON(filename)
module<-json_data$max_omega_logi
# module<-json_data$clustered
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


library(jsonlite)
json_data <- fromJSON("D:/Soil3DNetworks/code/cluster/Layer4/json/k172-1.json")
k <- 172
module<-json_data$max_omega_logi
length(table(module))
table(module)

par(mfrow=c(18,20),las=1)
par(oma = c(4, 4.5, 0.5, 0.5),mgp=c(1,0.5,0))

for (i in 1:k) {
  index<-which(module==(((1:k))-1)[i])
  
  y1 <- range(fit4$power_fit)
  lim2 <- y1[2]-y1[1]
  
  x1 <- range(fit4$Time)
  lim1 <- x1[2]-x1[1]
  
  par(mar = c(0, 0.8, 0,0))
  plot(1,1,type = "n",ylab = "",xlab = "",axes=FALSE,
       ylim = c(y1[1]-0.01*lim2,y1[2]+0.05*lim2),
       xlim = c(x1[1]-0.05*lim1,x1[2]+0.05*lim1))
  
  for (ii in index){
    # points(fit4$Time,fit4$original_data[,ii],col="#7BABDD",pch=16)
    points(index4,data_order4[,ii],
           col=rgb(243,238,175,100,maxColorValue = 255),pch=16,cex=0.6)
    # lines(fit4$Time,fit4$power_fit[,ii],
    #        col=rgb(195,223,237,100,maxColorValue = 255),lwd=2)
  }
  
  if(length(index)==1){
    par <- power_equation_base2(fit4$Time, fit4$original_data[,index])
    lines(fit4$Time,power.equation(par,fit4$Time),col="#eeac03",lwd=2)
  }else{
    par <- power_equation_base2(fit4$Time, rowMeans(fit4$original_data[,index]))
    lines(fit4$Time,power.equation(par,fit4$Time),col="#eeac03",lwd=2)
  }
  box(lwd=1)
  if(i%in%((k-10+1):k)){
    axis(1,at=seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3),
         labels =sprintf("%0.1f", seq(x1[1]+0.1*lim1,x1[2]-0.1*lim1,length=3)/10000),
         cex.axis=1.1,gap.axis = -1)
  }
  
  if(i%in%seq(1,k,10)){
    axis(2,at = seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4),
         labels =sprintf("%0.1f", seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4)),
         cex.axis=1.1)
  }else{
    axis(2,at = seq(y1[1]+0.05*lim2,y1[2]-0.05*lim2,length=4),
         labels =FALSE,
         cex.axis=1.1)
  }
  
  mi <- paste("M",i," (",length(index),")",sep = "")
  text(x1[1]+0.45*lim1,y1[2]-0.03*lim2,mi,cex=1.2)
  
  par(mar = c(0, 0, 0,0.8))
  plot(n1,n2, col = "gray90", pch = 16,ylab = "",xlab = "",axes=FALSE,
       ylim=c(31.6,45),cex=0.3)
  points(n1[index],n2[index], col = "red", pch = 16,cex=0.5)
  points(coords[,1],coords[,2],col="blue",cex=0.1,type = "p")
  box(lwd=1,lty=2)
  
  # text(105.9958,44,mi,cex=1.3)
  cat("完成 plot M",i,"\n")
  
}

mtext("Compartment Index (×10000)",side = 1,line = 2.8,cex=1.4,outer = TRUE)
mtext("Value of Soil Properties",side = 2,line = 2.5,cex=1.4,outer = TRUE,las=0)


dev.off()
