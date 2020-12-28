###### Case 2: Model without smoking
source("data_generation_case2.R")
rate_list = list(rate_lung1, rate_esophagus1, rate_larynx1, rate_colrect1)

orderindex = c(1,2,3,4)
models = permutations(n=4,r=4,v=orderindex,repeats.allowed=F)

mcmc_list = list()
mcmc_list2 = list()
for(i in 1:24){
  print(i)
  rate_list1 = rate_list[models[i,]]
  
  Y1 = rate_list1[[1]]$rate[final_perm]
  Y2 = rate_list1[[2]]$rate[final_perm]
  Y3 = rate_list1[[3]]$rate[final_perm]
  Y4 = rate_list1[[4]]$rate[final_perm]
  
  X1 = as.matrix(cbind(1,rate_list1[[1]][,6:13]))[final_perm,]
  X2 = as.matrix(cbind(1,rate_list1[[2]][,6:13]))[final_perm,]
  X3 = as.matrix(cbind(1,rate_list1[[3]][,6:13]))[final_perm,]
  X4 = as.matrix(cbind(1,rate_list1[[4]][,6:13]))[final_perm,]
  
  Y = c(Y1,Y2,Y3,Y4)
  X = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))
  
  tod1=Sys.time()
  
  sink("DAGAR.txt")
  cat("
      model
      {
      for(i in 1:2){
      Tau1[i, i] <- 1
      Tau2[i, i] <- 1
      Tau3[i, i] <- 1
      Tau4[i, i] <- 1
      
      for(j in 1:k){
      B1[i,j] <- 0
      B2[i,j] <- 0
      B3[i,j] <- 0
      B4[i,j] <- 0
      }
      }
      
      for(l in 1:(k-1)){
      for(h in (l+1):k){
      Tau1[l,h] <- 0
      Tau2[l,h] <- 0
      Tau3[l,h] <- 0
      Tau4[l,h] <- 0
      }
      }
      
      Tau1[2,1] <- 0
      Tau2[2,1] <- 0
      Tau3[2,1] <- 0
      Tau4[2,1] <- 0
      
      for (i in 3:k)
      {
      Tau1[i,i] <- (1 + (ns[i-1] - 1) * rho1^2) / (1 - rho1^2)
      Tau2[i,i] <- (1 + (ns[i-1] - 1) * rho2^2) / (1 - rho2^2)
      Tau3[i,i] <- (1 + (ns[i-1] - 1) * rho3^2) / (1 - rho3^2)
      Tau4[i,i] <- (1 + (ns[i-1] - 1) * rho4^2) / (1 - rho4^2)
      for(h in 1:(i-1)){
      Tau1[i,h] <- 0
      Tau2[i,h] <- 0
      Tau3[i,h] <- 0
      Tau4[i,h] <- 0
      }
      b1[i] <- rho1 / (1 + (ns[i-1] - 1) * rho1^2)
      b2[i] <- rho2 / (1 + (ns[i-1] - 1) * rho2^2)
      b3[i] <- rho3 / (1 + (ns[i-1] - 1) * rho3^2)
      b4[i] <- rho4 / (1 + (ns[i-1] - 1) * rho4^2)
      for(j in (udnei[(cn[i-1] + 1):(cn[i-1] + ns[i-1])])){
      B1[i,j] <- b1[i]
      B2[i,j] <- b2[i]
      B3[i,j] <- b3[i]
      B4[i,j] <- b4[i]
      }
      for(j in index1[((k)*(i-3)-cn[i-1]+1) : ((k)*(i-3)-cn[i-1] + (k - ns[i-1]))]){
      B1[i,j] <- 0
      B2[i,j] <- 0
      B3[i,j] <- 0
      B4[i,j] <- 0
      }
      }
      
      Q1 <- t(I - B1) %*% Tau1 %*% (I - B1)
      Q2 <- t(I - B2) %*% Tau2 %*% (I - B2)
      Q3 <- t(I - B3) %*% Tau3 %*% (I - B3)
      Q4 <- t(I - B4) %*% Tau4 %*% (I - B4)
      
      W1[1:k] ~ dmnorm(rep(0, k), tau1*Q1)
      W2[1:k] ~ dmnorm(rep(0, k), tau2*Q2)
      W3[1:k] ~ dmnorm(rep(0, k), tau3*Q3)
      W4[1:k] ~ dmnorm(rep(0, k), tau4*Q4)
      A21 <- eta021 * I + eta121 * Minc
      A31 <- eta031 * I + eta131 * Minc
      A32 <- eta032 * I + eta132 * Minc
      A41 <- eta041 * I + eta141 * Minc
      A42 <- eta042 * I + eta142 * Minc
      A43 <- eta043 * I + eta143 * Minc
      W[1:k] <- W1
      W[(k+1):(2*k)] <- A21 %*% W1 + W2
      W[(2*k+1):(3*k)] <- A31 %*% W1 + A32 %*% W[(k+1):(2*k)] + W3
      W[(3*k+1):(4*k)] <- A41 %*% W1 + A42 %*% W[(k+1):(2*k)] + A43 %*% W[(2*k+1):(3*k)] + W4
      for (i in 1:k)
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue1)
      loglik[i] <- logdensity.norm(Y[i], mu[i], taue1)
      }
      
      for (i in (k+1):(2*k))
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue2)
      loglik[i] <- logdensity.norm(Y[i], mu[i], taue2)
      }
      
      for (i in (2*k+1):(3*k))
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue3)
      loglik[i] <- logdensity.norm(Y[i], mu[i], taue3)
      }
      for (i in (3*k+1):(4*k))
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue4)
      loglik[i] <- logdensity.norm(Y[i], mu[i], taue4)
      }
      
      rho1 ~ dunif(0, 0.999)
      rho2 ~ dunif(0, 0.999)
      rho3 ~ dunif(0, 0.999)
      rho4 ~ dunif(0, 0.999)
      tau1 ~ dgamma(2, 0.1)
      tau2 ~ dgamma(2, 0.1)
      tau3 ~ dgamma(2, 0.1)
      tau4 ~ dgamma(2, 0.1)
      eta021 ~ dnorm(0, 0.01)
      eta121 ~ dnorm(0, 0.01)
      eta031 ~ dnorm(0, 0.01)
      eta131 ~ dnorm(0, 0.01)
      eta032 ~ dnorm(0, 0.01)
      eta132 ~ dnorm(0, 0.01)
      eta041 ~ dnorm(0, 0.01)
      eta141 ~ dnorm(0, 0.01)
      eta042 ~ dnorm(0, 0.01)
      eta142 ~ dnorm(0, 0.01)
      eta043 ~ dnorm(0, 0.01)
      eta143 ~ dnorm(0, 0.01)
      
      taue1 ~ dgamma(2, 1)
      taue2 ~ dgamma(2, 1)
      taue3 ~ dgamma(2, 1)
      taue4 ~ dgamma(2, 1)
      vare1 <- 1/taue1
      vare2 <- 1/taue2
      vare3 <- 1/taue3
      vare4 <- 1/taue4
      beta[1:36] ~ dmnorm(rep(0,36), (0.001*I1))
      }
      ", fill = TRUE)
  sink()
  
  model.data <- list(k = n, index1 = index1, I = diag(n), I1 = diag(36), Minc = Minc, ns = dni, cn = c(0, cni), udnei = udnei, X = X, Y = Y)
  model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, rho3 = 0.1, tau1 = 1, tau2 = 1, tau3 = 1, eta021 = 1, 
                               eta121 = 1, eta031 = 1, eta131 = 1, eta032 = 1, eta132 = 1, eta041 = 1, 
                               eta141 = 1, eta042 = 1, eta142 = 1, eta043 = 1, eta143 = 1, taue1 = 1, taue2 = 1,
                               taue3 = 1, taue4 = 1, beta = rep(0, 36), W1 = rep(0.1, n), 
                               W2 = rep(0.1, n), W3 = rep(0.1, n), W4 = rep(0.1, n))),2)
  model.param1 <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", "eta021", "eta121",
                  "eta031", "eta131", "eta032", "eta132","eta041", "eta141",
                  "eta042", "eta142", "eta043", "eta143", "vare1", "vare2", "vare3", "vare4", "W", "loglik")
  model.param2 <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", "eta021", "eta121",
                    "eta031", "eta131", "eta032", "eta132","eta041", "eta141",
                    "eta042", "eta142", "eta043", "eta143", "vare1", "vare2", "vare3", "vare4")
  set.seed(125)
  result1 <- jags(model.data, model.inits, model.param1, "DAGAR.txt",
                  n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  result2 <- jags(model.data, model.inits, model.param2, "DAGAR.txt",
                  n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  tod2=Sys.time()
  
  mcmc_list[[i]] = result1
  mcmc_list2[[i]] = result2
}

saveRDS(mcmc_list, "mcmc_list.rds")
saveRDS(mcmc_list2, "mcmc_list2.rds")

mcmc_select <- mcmc_list[[16]]
jags_result1 = mcmc_select$BUGSoutput$sims.array[,1,c(233:281,514:525)]
jags_result2 = mcmc_select$BUGSoutput$sims.array[,2,c(233:281,514:525)]
jags_result <- rbind(jags_result1, jags_result2)

colnames(jags_result) = c("beta1", "beta2","beta3","beta4","beta5", "beta6", "beta7","beta8","beta9","beta10","beta11","beta12",
                          "beta13", "beta14","beta15","beta16","beta17", "beta18", "beta19","beta20","beta21","beta22","beta23","beta24",
                          "beta25", "beta26","beta27","beta28","beta29", "beta30", "beta31","beta32","beta33","beta34","beta35","beta36",
                          "deviance", "eta0_21", "eta0_31", "eta0_32", "eta0_41",
                          "eta0_42", "eta0_43", "eta1_21", "eta1_31", "eta1_32", 
                          "eta1_41", "eta1_42", "eta1_43", "rho1", "rho2", "rho3", "rho4",
                          "tausq1", "tausq2", "tausq3", "tausq4","sigmasq1","sigmasq2","sigmasq3","sigmasq4")
estimate1 <- round(t(apply(jags_result, 2, mysummary)),2)

table_est = paste(estimate1[,1], " (", estimate1[,4], ", ", estimate1[,5], ")", sep="")
nrow(estimate1)

table1 <- data.frame(matrix(table_est[1:36], ncol = 4, byrow=F))
table2 <- data.frame(matrix(table_est[50:61], ncol = 4, byrow=T))

table <- rbind(table1, table2)
colnames(table) <- c("Larynx", "Esophageal", "Colorectum", "Lung")


############ correlation between pairwise diseases in same region ##########
Dinv_new <- function(Rho, n, cn, ns, udnei,q){
  Tau <- list()
  B <- list()
  invD <- list()
  for(i in 1:q){
    Tau[[i]] <- diag(n)
    B[[i]] <- matrix(0, n, n)
    for(j in 3:n){
      Tau[[i]][j,j] <- (1 + (ns[j-1] - 1) * Rho[i]^2) / (1 - Rho[i]^2)
      for(k in (udnei[(cn[j-1] + 1):(cn[j-1] + ns[j-1])])){
        B[[i]][j,k] <- Rho[i] / (1 + (ns[j-1] - 1) * Rho[i]^2)
      }
    }
    invD[[i]] <- t(diag(n) - B[[i]]) %*% Tau[[i]] %*% (diag(n) - B[[i]])
  }
  return(invD)
}
cor_fun = function(x){
  #x = data.frame(param.est)
  A21_est = diag(x[["eta0_21"]], n) + x[["eta1_21"]] * Minc
  A31_est = diag(x[["eta0_31"]], n) + x[["eta1_31"]] * Minc
  A32_est = diag(x[["eta0_32"]], n) + x[["eta1_32"]] * Minc
  A41_est = diag(x[["eta0_41"]], n) + x[["eta1_41"]] * Minc
  A42_est = diag(x[["eta0_42"]], n) + x[["eta1_42"]] * Minc
  A43_est = diag(x[["eta0_43"]], n) + x[["eta1_43"]] * Minc
  
  invD_est = Dinv_new(c(x[["rho1"]],x[["rho2"]],x[["rho3"]],x[["rho4"]]), n, cn, ns, udnei,q=4)
  G_est = as.matrix(bdiag(bdiag(1/x[["tausq1"]]*solve(invD_est[[1]]), 1/x[["tausq2"]]*solve(invD_est[[2]])), 
                          bdiag(1/x[["tausq3"]]*solve(invD_est[[3]]), 1/x[["tausq4"]]*solve(invD_est[[4]]))))
  A_est = as.matrix(blockmatrix(names = c("0","A21","A31","A41","0","0","A32","A42","0","0","0","A43","0","0","0","0"), 
                                A21=A21_est, A31=A31_est, A32=A32_est, A41=A41_est, A42=A42_est, A43=A43_est, dim=c(4,4)))
  L_est = solve(diag(4*n)-A_est)
  WV_est = L_est%*%G_est%*%t(L_est)
  Sigma_est = diag(c(rep(x[["sigmasq1"]],n),rep(x[["sigmasq2"]],n),rep(x[["sigmasq3"]],n),rep(x[["sigmasq4"]],n))) 
  V_est = WV_est + Sigma_est
  
  
  cor12 = diag(V_est[1:58, (1+58):(58+58)])/sqrt(diag(V_est[1:58, 1:58])*diag(V_est[(1+58):(58+58), (1+58):(58+58)]))
  cor13 = diag(V_est[1:58, (1+58*2):(58+58*2)])/sqrt(diag(V_est[1:58, 1:58])*diag(V_est[(1+58*2):(58+58*2), (1+58*2):(58+58*2)]))
  cor14 = diag(V_est[1:58, (1+58*3):(58+58*3)])/sqrt(diag(V_est[1:58, 1:58])*diag(V_est[(1+58*3):(58+58*3), (1+58*3):(58+58*3)]))
  cor23 = diag(V_est[(1+58):(58+58), (1+58*2):(58+58*2)])/sqrt(diag(V_est[(1+58):(58+58), (1+58):(58+58)])*diag(V_est[(1+58*2):(58+58*2), (1+58*2):(58+58*2)]))
  cor24 = diag(V_est[(1+58):(58+58), (1+58*3):(58+58*3)])/sqrt(diag(V_est[(1+58):(58+58), (1+58):(58+58)])*diag(V_est[(1+58*3):(58+58*3), (1+58*3):(58+58*3)]))
  cor34 = diag(V_est[(1+58*2):(58+58*2), (1+58*3):(58+58*3)])/sqrt(diag(V_est[(1+58*2):(58+58*2), (1+58*2):(58+58*2)])*diag(V_est[(1+58*3):(58+58*3), (1+58*3):(58+58*3)]))
  
  return(c(cor12,cor13,cor14,cor23,cor24,cor34))
}
mysummary = function(invector) {
  c(mean(invector), median(invector), sd(invector), quantile(invector, .025), quantile(invector,.975))
}

cor_gibbs = apply(jags_result,1,cor_fun)
cor_gibbs_est = apply(cor_gibbs,1,mysummary)

cor_gibbs_est12 = t(cor_gibbs_est[,1:58])
cor_gibbs_est13 = t(cor_gibbs_est[,(1+58):(58+58)])
cor_gibbs_est14 = t(cor_gibbs_est[,(1+58*2):(58+58*2)])
cor_gibbs_est23 = t(cor_gibbs_est[,(1+58*3):(58+58*3)])
cor_gibbs_est24 = t(cor_gibbs_est[,(1+58*4):(58+58*4)])
cor_gibbs_est34 = t(cor_gibbs_est[,(1+58*5):(58+58*5)])

ca.poly$rate12 = cor_gibbs_est12[,1][order(final_perm)]
ca.poly$rate13 = cor_gibbs_est13[,1][order(final_perm)]
ca.poly$rate14 = cor_gibbs_est14[,1][order(final_perm)]
ca.poly$rate23 = cor_gibbs_est23[,1][order(final_perm)]
ca.poly$rate24 = cor_gibbs_est24[,1][order(final_perm)]
ca.poly$rate34 = cor_gibbs_est34[,1][order(final_perm)]
ca.poly$sig12 = 0
ca.poly$sig12[which(cor_gibbs_est12[,4][order(final_perm)] > 0 | cor_gibbs_est12[,5][order(final_perm)] < 0)] = 1
ca.poly$sig13 = 0
ca.poly$sig13[which(cor_gibbs_est13[,4][order(final_perm)] > 0 | cor_gibbs_est13[,5][order(final_perm)] < 0)] = 1
ca.poly$sig14 = 0
ca.poly$sig14[which(cor_gibbs_est14[,4][order(final_perm)] > 0 | cor_gibbs_est14[,5][order(final_perm)] < 0)] = 1
ca.poly$sig23 = 0
ca.poly$sig23[which(cor_gibbs_est23[,4][order(final_perm)] > 0 | cor_gibbs_est23[,5][order(final_perm)] < 0)] = 1
ca.poly$sig24 = 0
ca.poly$sig24[which(cor_gibbs_est24[,4][order(final_perm)] > 0 | cor_gibbs_est24[,5][order(final_perm)] < 0)] = 1
ca.poly$sig34 = 0
ca.poly$sig34[which(cor_gibbs_est34[,4][order(final_perm)] > 0 | cor_gibbs_est34[,5][order(final_perm)] < 0)] = 1


brks_fit = c(-0.5, 0.1, 0.3, 0.5, 0.7, 0.8)
color.pallete = brewer.pal(5,"Blues")
class.rate12 = classIntervals(var=ca.poly$rate12, n=5, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate13 = classIntervals(var=ca.poly$rate13, n=5, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate14 = classIntervals(var=ca.poly$rate14, n=5, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate23 = classIntervals(var=ca.poly$rate23, n=5, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate24 = classIntervals(var=ca.poly$rate24, n=5, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate34 = classIntervals(var=ca.poly$rate34, n=5, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
color.code.rate12 = findColours(class.rate12, color.pallete)
color.code.rate13 = findColours(class.rate13, color.pallete)
color.code.rate14 = findColours(class.rate14, color.pallete)
color.code.rate23 = findColours(class.rate23, color.pallete)
color.code.rate24 = findColours(class.rate24, color.pallete)
color.code.rate34 = findColours(class.rate34, color.pallete)

pdf("correlation_pairwise_case2.pdf", height = 10, width = 10)
par(mfrow=c(3,2), oma = c(4, 1, 1, 1))
leg.txt = c("< 0.1", "0.1 - 0.3", "0.3 - 0.5", "0.5 - 0.7", "0.7 - 0.9", "0.9 - 1") 

plot(ca.poly, col = color.code.rate12)
title(sub="(a) larynx cancer and esophageal cancer")
points(ca.coords[which(ca.poly$sig12==1),1], ca.coords[which(ca.poly$sig12==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate13)
title(sub="(b) larynx cancer and colorectum cancer" )
points(ca.coords[which(ca.poly$sig14==1),1], ca.coords[which(ca.poly$sig14==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate14)
title(sub="(c) larynx cancer and lung cancer" )
points(ca.coords[which(ca.poly$sig14==1),1], ca.coords[which(ca.poly$sig14==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate23)
title(sub="(d) esophageal cancer and colorectum cancer" )
points(ca.coords[which(ca.poly$sig23==1),1], ca.coords[which(ca.poly$sig23==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate24)
title(sub="(e) esophageal cancer and lung cancer" )
points(ca.coords[which(ca.poly$sig24==1),1], ca.coords[which(ca.poly$sig24==1),2], pch=10, col="yellow", cex=0.2)


plot(ca.poly, col = color.code.rate34)
title(sub="(f) colorectum cancer and lung cancer" )
points(ca.coords[which(ca.poly$sig34==1),1], ca.coords[which(ca.poly$sig34==1),2], pch=10, col="yellow", cex=0.2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("center", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

dev.off()