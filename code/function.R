power_equation_fit2 <- function(x,y,X_smooth, thread = 4) {
  data <- y
  X <- x
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterExport(cl, c("power.equation", "power_equation_base2", "data", "X","s.mle"), envir = environment())
  all_model = parLapply(cl = cl, 1:ncol(data), function(c) power_equation_base2(X, data[,c]))
  stopCluster(cl)
  
  
  names(all_model) = colnames(data)
  
  # new_x = seq(min(X), max(X), length = n)
  new_x = X
  new_x2 = X_smooth
  power_fit = sapply(1:ncol(data), function(c)power.equation(all_model[[c]],new_x))
  power_fit2 = sapply(1:ncol(data), function(c)power.equation(all_model[[c]],new_x2))
  colnames(power_fit) <- colnames(data)
  colnames(power_fit2) <- colnames(data)
  
  result = list(original_data = data, power_par = all_model, power_fit = power_fit,
                Time = X, smooth_fit = power_fit2, smooth_Time = X_smooth)
  return(result)
}

power_equation_base2 <- function(x, y){
  x <- as.numeric(x)
  y <- as.numeric(y)
  lmFit <- lm( log(y + rep(1e-10,length(y)))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  m1 <- list()
  m1[[1]] <- try(optim(c(a,b),s.mle,data=y,x=x,method = "Nelder-Mead"),TRUE)
  if (class(m1[[1]]) == "try-error"){
    m1[[1]] <- list(par=c(1,1),value=1000)
  }
  m1[[2]] <- optim(c(0.1,0.1),s.mle,data=y,x=x,method = "Nelder-Mead")
  m1[[3]] <- optim(c(0.5,0.5),s.mle,data=y,x=x,method = "Nelder-Mead")
  m1[[4]] <- optim(c(0,0),s.mle,data=y,x=x,method = "Nelder-Mead")
  min_index <- which.min(c(m1[[1]]$value,m1[[2]]$value,m1[[3]]$value,m1[[4]]$value))
  return(optim(m1[[min_index]]$par,s.mle,data=y,x=x,method = "Nelder-Mead")$par)
}


power.equation <- function(par,x){
  par[1]*x^par[2]
}

s.mle<-function(par,data,x){
  y <- power.equation(par,x)
  yi <- data
  res <- sum((yi-y)^2)
  return(res)
}

get_interaction <- function(data, col=1, depth=1){
  
  min_max_norm <- function(x) {
    return(x / (147.5))
  }
  
  dis_all <- rep(distance,each=15)
  clean_data <- data
  n <- nrow(clean_data)
  colnames(clean_data) <- as.vector(t(outer((1:6), 1:15, 
                                            FUN = function(x, y) paste(x, y, sep = "-"))))
  
  ind_no <- which(colnames(clean_data)==paste(depth, col, sep = "-"))
  dep_no <- which(as.numeric(substr(colnames(clean_data),1,1))>=depth)
  dep_no2 <- setdiff(dep_no,ind_no)
  m <- clean_data[,ind_no]
  x_matrix <- clean_data[,dep_no2]
  
  ridge1_cv <- cv.glmnet(x = x_matrix, y = m,alpha = 0)
  best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
  
  dis <- abs(dis_all[ind_no]-dis_all[dep_no2])
  dis <- min_max_norm(dis)+1
  
  # pf <- dis+(1/(best_ridge_coef))
  pf <- (dis)^0.5/(best_ridge_coef)^0.5
  # pf <- pf / mean(pf)
  
  fit_res <- cv.glmnet(x = x_matrix, y = m,alpha = 0.5,
                       penalty.factor = pf,
                       keep = TRUE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
  best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
  
  gene_list_one <- list()
  gene_list_one[[1]] <- paste(depth, col, sep = "-")
  gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
  gene_list_one[[3]] <- best_alasso_coef1@x[-1]
  
  return(gene_list_one)
}


# get_interaction <- function(data, col=1, depth=1){
#   
#   dis_all <- rep(distance,each=15)
#   clean_data <- data
#   n <- nrow(clean_data)
#   colnames(clean_data) <- as.vector(t(outer((1:6), 1:15, 
#                                    FUN = function(x, y) paste(x, y, sep = "-"))))
#   
#   ind_no <- which(colnames(clean_data)==paste(depth, col, sep = "-"))
#   m <- clean_data[,ind_no]
#   x_matrix <- clean_data[,-ind_no]
#   
#   ridge1_cv <- cv.glmnet(x = x_matrix, y = m,alpha = 0)
#   best_ridge_coef <- abs(as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1])
#   
#   dis <- abs(dis_all[ind_no]-dis_all[-ind_no])
#   dis[which(dis==0)] <- 1
#   # dis <- dis^0.5
#   
#   pf <- dis+(1/(best_ridge_coef))
#   # pf <- (dis/(best_ridge_coef))^0.5
#   pf <- pf / mean(pf)
#   
#   fit_res <- cv.glmnet(x = x_matrix, y = m,alpha = 0.5,
#                        penalty.factor = pf,
#                        keep = TRUE)
#   best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
#   
#   gene_list_one <- list()
#   gene_list_one[[1]] <- paste(depth, col, sep = "-")
#   gene_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]
#   gene_list_one[[3]] <- best_alasso_coef1@x[-1]
#   
#   return(gene_list_one)
# }


get_effect <- function(pars,time,power_par,order,y0){
  #Legendre polynomials
  LOP <-  legendre.polynomials(order, normalized=F)
  LOP_fit <-  sapply(1:length(pars),function(c) pars[c]*LOP[[c]])
  if(order==0){
    f <- function(x){LOP_fit}
  }else{
    f <- function(x){dy=do.call(sum,polynomial.values(polynomials=LOP_fit,x=x));dy}
  }
  
  max_y <- time[length(time)]
  min_y <- time[1]
  
  rescale_y <- function(yi){
    return((yi - min_y)/(max_y - min_y))
  }
  
  h <- diff(time)
  
  dy_fit <- c(y0)
  for (j in 1:(length(time) - 1)) {
    k1 <- f(rescale_y(time[j])) * power.equation(power_par,time[j])
    k2 <- f(rescale_y(time[j] + h[j] / 2)) * power.equation(power_par,(time[j] + h[j] / 2))
    k3 <- f(rescale_y(time[j] + h[j] / 2)) * power.equation(power_par,(time[j] + h[j] / 2))
    k4 <- f(rescale_y(time[j] + h[j])) * power.equation(power_par,(time[j] + h[j]))
    y <- dy_fit[j]+h[j]/6*(k1+2*(1-1/sqrt(2))*k2+2*(1+1/sqrt(2))*k3+k4)
    dy_fit <- c(dy_fit,y)
  }
  return(dy_fit)
}



#' @title calculate least-square for observed and fitted data
#' @param pars matrix of LOP parameters for ind and dep growth curve
#' @param ind the independent growth curve id
#' @param dep the dependent growth curve id
#' @param times vector of time point
#' @param data dataframe of observed data
#' @param order scalar of LOP order
#' @param effect matrix of observed data
#' @return scalar of least-square error
ode_optimize <- function(pars,ind,dep,effecti,power_par,time,smooth_time,order_ind,order_dep){
  if (all(is.na(dep))){
    ind_pars <- pars
    inital_value <- power.equation(power_par[[ind]],time[1])
    ind_effect <- get_effect(ind_pars,time,power_par[[ind]],order_ind,inital_value)
    y <- ind_effect
    ssr <- sum((effecti[,ind]-y)^2)
    smooth_ind_effect <- get_effect(ind_pars,smooth_time,power_par[[ind]],
                                    order_ind,inital_value)
    smooth_y <- smooth_ind_effect
  }else{
    ind_pars <- pars[1:(order_ind+1)]
    dep_pars <- matrix(pars[-(1:(order_ind+1))],ncol=(order_dep+1))
    # inital_value <- effecti[1,ind]
    inital_value <- power.equation(power_par[[ind]],time[1])
    ind_effect <- get_effect(ind_pars,time,power_par[[ind]],order_ind,inital_value)
    if (nrow(dep_pars)==1) {
      dep_effect <- get_effect(dep_pars,time,power_par[[dep]],order_dep,0)
      y <- ind_effect+dep_effect
    }else{
      dep_effect <- sapply(1:length(dep), function(c)
        get_effect(dep_pars[c,],time,power_par[[dep[c]]],order_dep,0))
      y <- ind_effect+rowSums(dep_effect)
    }
    ssr <- sum((effecti[,ind]-y)^2)
    
  #   smooth_ind_effect <- get_effect(ind_pars,smooth_time,power_par[[ind]],
  #                                   order_ind,inital_value)
  #   if (nrow(dep_pars)==1) {
  #     smooth_dep_effect <- get_effect(dep_pars,smooth_time,power_par[[dep]],order_dep,0)
  #     smooth_y <- smooth_ind_effect+smooth_dep_effect
  #   }else{
  #     smooth_dep_effect <- sapply(1:length(dep), function(c)
  #       get_effect(dep_pars[c,],smooth_time,power_par[[dep[c]]],order_dep,0))
  #     smooth_y <- smooth_ind_effect+rowSums(smooth_dep_effect)
  #   }
  }
  
  #add penalty
  # alpha=length(which(smooth_ind_effect<0))
  alpha=1
  alpha1=0
  # if(max(abs(pars))<1){
  #   return(ssr+alpha1*sum(pars^2)+alpha*sum(pmax(0,-smooth_ind_effect))+alpha*sum(pmax(0,-smooth_y)))
  # }else{
  #   return(ssr+alpha1*sum(pars^2)+alpha*sum(pmax(0,-smooth_ind_effect))+alpha*sum(pmax(0,-smooth_y)))+1e10
  # }
  
  penalty <- sum(pmax(0, abs(pars) - 1)^2) * 1e6
  # return(ssr + alpha1*sum(pars^2) + alpha*sum(pmax(0, -smooth_ind_effect)) + 
  #          alpha*sum(pmax(0, -smooth_y)))
  return(ssr + alpha1*sum(pars^2) + alpha*sum(pmax(0, -ind_effect)) + alpha*sum(pmax(0, -y)))
}

get_value <- function(effecti,power_par,time,smooth_time,relationship,order_ind,order_dep){
  #input
  ind <- relationship[[1]]
  dep <- relationship[[2]]
  ind_no <- as.numeric(which(colnames(effecti)==ind))
  dep_no <- as.numeric(sapply(1:length(dep), function(c) which(colnames(effecti)==dep[c])))
  if(length(dep)==0){
    init_pars <- rep(0.1,(length(ind_no)*(order_ind+1)))
  }else{
    init_pars <- rep(0.1,(length(ind_no)*(order_ind+1)+length(dep_no)*(order_dep+1)))
  }
  result <- optim(init_pars,ode_optimize,ind=ind_no,dep=dep_no,effecti=effecti,
                  power_par=power_par,time=time,smooth_time=smooth_time,
                  order_ind=order_ind,order_dep=order_dep,
                  method = "Nelder-Mead", control=list(maxit=50000,trace=T))
  return(result$par)
}



#' @title helper function to convert ODE result
#' @param relationship list contain the result of lasso-based variable election
#' @param par vector conatin LOP parameters
#' @param effect matrix of observed data
#' @param times vector of time point
#' @param order scalar of LOP order
#' @return list contain bunch of useful output result for a row
get_ode_output <- function(relationship, par, effecti, power_par,time, order_ind,order_dep){
  #effect <- t(effect)
  output <- list()
  output[[1]] <- relationship[[1]]
  output[[2]] <- relationship[[2]]
  if(length(relationship[[2]])==0){
    output[[3]] <- par
    output[[4]] <- 0
  }else{
    output[[3]] <- par[1:(order_ind+1)]
    output[[4]] <- matrix(par[-(1:(order_ind+1))],ncol=(order_dep+1))
  }
  
  ind_no <- as.numeric(which(colnames(effecti)==output[[1]]))
  dep_no <- as.numeric(sapply(1:length(output[[2]]),
                              function(c) which(colnames(effecti)==output[[2]][c])))
  # inital_value <- effecti[1,ind_no]
  inital_value <- power.equation(power_par[[ind_no]],time[1])
  ind_effect <- get_effect(as.numeric(output[[3]]),time,power_par[[ind_no]],order_ind,inital_value)
  
  if (all(is.na(dep_no))){
    dep_effect <- rep(0,length(ind_effect))
  }else{
    if (length(dep_no)==1) {
      dep_effect <- get_effect(as.numeric(output[[4]]),time,power_par[[dep_no]],order_dep,0)
    }else{
      dep_effect <- sapply(1:length(dep_no), function(c)
        get_effect(as.numeric(output[[4]][c,]),time,power_par[[dep_no[c]]],order_dep,0))
      colnames(dep_effect) <- dep_no
    }
  }
  
  
  all_effect <- cbind(ind_effect,dep_effect)
  effect_mean <- apply(all_effect,2,mean)
  output[[5]] <- effect_mean
  output[[6]] <- all_effect
  return(output)
}