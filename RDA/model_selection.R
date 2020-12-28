setwd("Multivariate_DAGAR/RDA")
#Import data for case 1
source("data_generation_case1.R")

library(bridgesampling)
library(mnormt)

##################################################

orderindex = c(1,2,3,4)
models = permutations(n=4,r=4,v=orderindex,repeats.allowed=F)

Ylist = list()
Xlist = list()

for(i in 1:24){
  print(i)
  
  rate_list1 = rate_list[models[i,]]
  
  Y1 = rate_list1[[1]]$rate[final_perm]
  Y2 = rate_list1[[2]]$rate[final_perm]
  Y3 = rate_list1[[3]]$rate[final_perm]
  Y4 = rate_list1[[4]]$rate[final_perm]
  
  X1 = as.matrix(cbind(1,rate_list1[[1]][,6:14]))[final_perm,]
  X2 = as.matrix(cbind(1,rate_list1[[2]][,6:14]))[final_perm,]
  X3 = as.matrix(cbind(1,rate_list1[[3]][,6:14]))[final_perm,]
  X4 = as.matrix(cbind(1,rate_list1[[4]][,6:14]))[final_perm,]
  
  Ylist[[i]] = c(Y1,Y2,Y3,Y4)
  Xlist[[i]] = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))
}

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

### using jags output
mcmc_list_jags = readRDS("mcmc_list_new2.rds")

### Functions used for bridge sampling

#likelihood function
Likli_fun = function(samples.row, data){
  
  A21_est = diag(samples.row["eta021"], data$n) + samples.row["eta121"] * data$Minc
  A31_est = diag(samples.row["eta031"], data$n) + samples.row["eta131"] * data$Minc
  A32_est = diag(samples.row["eta032"], data$n) + samples.row["eta132"] * data$Minc
  A41_est = diag(samples.row["eta041"], data$n) + samples.row["eta141"] * data$Minc
  A42_est = diag(samples.row["eta042"], data$n) + samples.row["eta142"] * data$Minc
  A43_est = diag(samples.row["eta043"], data$n) + samples.row["eta143"] * data$Minc
  
  invD_est = Dinv_new(c(samples.row["rho1"],samples.row["rho2"],samples.row["rho3"],samples.row["rho4"]), data$n, data$cn, data$ns, data$udnei, q=4)
  G_est = as.matrix(bdiag(bdiag(1/samples.row["tau1"]*solve(invD_est[[1]]), 1/samples.row["tau2"]*solve(invD_est[[2]])), 
                          bdiag(1/samples.row["tau3"]*solve(invD_est[[3]]), 1/samples.row["tau4"]*solve(invD_est[[4]]))))
  A_est = as.matrix(blockmatrix(names = c("0","A21","A31","A41","0","0","A32","A42","0","0","0","A43","0","0","0","0"), 
                                A21=A21_est, A31=A31_est, A32=A32_est, A41=A41_est, A42=A42_est, A43=A43_est, dim=c(4,4)))
  L_est = solve(diag(4*data$n)-A_est)
  V_est = L_est%*%G_est%*%t(L_est)
  Sigma = diag(c(rep(samples.row["vare1"],data$n),rep(samples.row["vare2"],data$n),rep(samples.row[["vare3"]],data$n),rep(samples.row[["vare4"]],data$n)))
  Cov = V_est + Sigma
  Prec = solve(nearPD(Cov)$mat)
  beta = NULL
  for(par in 1:40){
    beta = c(beta, samples.row[paste("beta[", par, "]", sep="")])
  }
  mu = data$X%*%beta
  as.numeric(-2*data$n*log(2*pi) + 1/2*logdet(nearPD(Prec)$mat) - 1/2*t(data$Y-mu)%*%Prec%*%(data$Y-mu))
}

#unnormalized joint posterior function
log_posterior <- function(samples.row, data) {
  
  eta_sum = dnorm(samples.row["eta021"], 0, 10, log=T) + dnorm(samples.row["eta121"], 0, 10, log=T) +
    dnorm(samples.row["eta031"], 0, 10, log=T) + dnorm(samples.row["eta131"], 0, 10, log=T) +
    dnorm(samples.row["eta032"], 0, 10, log=T) + dnorm(samples.row["eta132"], 0, 10, log=T) +
    dnorm(samples.row["eta041"], 0, 10, log=T) + dnorm(samples.row["eta141"], 0, 10, log=T) +
    dnorm(samples.row["eta042"], 0, 10, log=T) + dnorm(samples.row["eta142"], 0, 10, log=T) +
    dnorm(samples.row["eta043"], 0, 10, log=T) + dnorm(samples.row["eta143"], 0, 10, log=T) 
  
  
  rho_sum = dunif(samples.row["rho1"], 0, 0.999, log=T) + dunif(samples.row["rho2"], 0, 0.999, log=T) +
    dunif(samples.row["rho3"], 0, 0.999, log=T) + dunif(samples.row["rho4"], 0, 0.999, log=T)
  
  tau_sum = dgamma(samples.row["tau1"], 2, rate = 0.1, log=T) + dgamma(samples.row["tau2"], 2, rate = 0.1, log=T) + 
    dgamma(samples.row["tau3"], 2, rate = 0.1, log=T) + dgamma(samples.row["tau4"], 2, rate = 0.1, log=T)
  
  sigma_sum = dinvgamma(samples.row["vare1"], 2, 1, log=T) + dinvgamma(samples.row["vare2"], 2, 1, log=T) +
    dinvgamma(samples.row["vare3"], 2, 1, log=T) + dinvgamma(samples.row["vare4"], 2, 1, log=T)
  
  beta = NULL
  for(par in 1:40){
    beta = c(beta, samples.row[paste("beta[", par, "]", sep="")])
  }
  beta_sum = dmnorm(beta, 0, diag(sqrt(1000),40), log=T)
  as.numeric(Likli_fun(samples.row, data) + eta_sum + rho_sum + tau_sum + sigma_sum + beta_sum)
}

#upper and lower bound for parameters
cname <- colnames(mcmc_list_jags[[1]]$BUGSoutput$sims.array[,1,c(1:40, 42:65)])
lb <- rep(-Inf, length(cname))
ub <- rep(Inf, length(cname))
names(lb) <- names(ub) <- cname
lb[53:64] <- 0
ub[53:56] <- 0.999

#bridge sampling for log marginal likelihood
model.bridge = list()
for(model in 1:24){
  print(model)
  set.seed(12345)
  model.bridge[[model]] <- bridge_sampler(samples = mcmc_list_jags[[model]], 
                                          data = list(Y = Ylist[[model]], X = Xlist[[model]], Minc = Minc, n = n, ns = dni, cn = c(0, cni), udnei = udnei),
                                          log_posterior = log_posterior, lb = lb,
                                          ub = ub, maxiter = 100000, silent = FALSE)
  
}
logml = rep(0, 24)
for(i in 1:24){
  logml[i] = model.bridge[[i]]$logml
}
#posterior model probabilities: model 10 gets the largest probability
round(exp(logml)/sum(exp(logml)),4)
saveRDS(model.bridge, "model_bridge_jags.rds")

