library(lattice)
library(mvtnorm)
library(matrixStats)
library(fields)
library(boot)
library(maps)
library(mapproj)
library(blockmatrix)
library(ggmcmc)
library(Matrix)
library(mcmc)
library(magic)
library(msos)
library(AICcmodavg)
library(coda)
library(invgamma)
library(MASS)
library(LaplacesDemon)
library(gtools)
library(rjags)
library(R2jags)
library(rjmcmc)
library(mnormt)
library(bridgesampling)

setwd("Multivariate_DAGAR/Simulation")
inclattice=function(m){
  n=m^2
  Minc=matrix(0,n,n)
  for(i in 1:(m-1))	for(j in 1:(m-1)) Minc[(i-1)*m+j,(i-1)*m+j+1]=Minc[(i-1)*m+j,i*m+j]=1
  for(i in 1:(m-1)) Minc[(i-1)*m+m,i*m+m]=1
  for(j in 1:(m-1)) Minc[(m-1)*m+j,(m-1)*m+j+1]=1
  Minc+t(Minc)
}

q = 3

graphnum=3

#### setting up graph: Minc, dmat and related variables #####

graph="usa"
data("state.fips")
states=setdiff(unique(state.fips$abb),"DC")

### importing the neighborhood states information
rawn <- scan("usneighbors.txt", what="", sep="\n")
rawn2 <- strsplit(rawn, "[[:space:]]+")
rawn2 <- lapply(rawn2, `[`, -1)
neighbors=lapply(rawn2,as.numeric)

n=length(neighbors)

Minc=sapply(neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Minc)=states

abborder=sapply(states,function(x) which(state.abb==x))
avglatlong=cbind(state.center$y,state.center$x)
avglatlong=avglatlong[abborder,]

### calculating minimum and maximum latitudes
### for Alber's equal area conic
latrange=round(quantile(avglatlong[,1],c(0.25,0.75)))

### Alber's equal area projection
albersproj=mapproject(avglatlong[,2],avglatlong[,1],projection = "albers",param=latrange)
projmat=cbind(albersproj$x,albersproj$y)

perm=order(albersproj$x+albersproj$y)
Minc=Minc[perm,perm]
projmat=projmat[perm,]

dmat=as.matrix(dist(projmat))
dmat=dmat/mean(dmat[which(Minc==1)])
n=nrow(Minc)
ni=rowSums(Minc)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Minc[x,]==1))
#N(i): 2:n
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n
dni=sapply(dneighbors,length)
nmax=max(dni)
cni=cumsum(dni)
dneimat=sapply(dneighbors, function(nei,nmax,n) c(nei,rep(n+1,nmax+1-length(nei))),nmax,n)
udnei=unlist(dneighbors)
for(i in 1:n) neimat[i,1:ni[i]]=neighbors[[i]]
intersectmat=matrix(0,n^2,maxn)
nijvec=rep(0,n^2)
for(i in 1:n) for(j in 1:n)
{ neighij=intersect(neighbors[[i]],neighbors[[j]])
nij=length(neighij)
nijvec[n*(i-1)+j]=nij
if(nij >0) intersectmat[n*(i-1)+j,1:nij]=neighij
}

#For neighbors without order
ni_wo = sapply(neighbors,length)
cni_wo = cumsum(ni_wo)
udnei_wo = unlist(neighbors)

region = seq(1:n)
cn = c(0, cni)
ns = dni
index = list()
for(i in 1:(n-1)){
  index[[i]] = region[-(udnei[(cn[i] + 1):(cn[i] + ns[i])])]
}
index1 = unlist(index)

# Assign true values for parameters
rho = c(0.2, 0.8, 0.5)
eta0_21 = 0.5
eta1_21 = 0.3
eta0_31 = 1
eta1_31 = 0.6
eta0_32 = 1.5
eta1_32 = 0.9
tau = c(0.25, 0.25, 0.25)

#generate each convariance matrix from Gaussian Process
G1 = rho[1]^dmat
G2 = rho[2]^dmat
G3 = rho[3]^dmat

A21 = diag(eta0_21, n) + eta1_21 * Minc
A31 = diag(eta0_31, n) + eta1_31 * Minc
A32 = diag(eta0_32, n) + eta1_32 * Minc

A = as.matrix(blockmatrix(names = c("0","A21","A31","0","0","A32","0","0","0"), 
                          A21=A21, A31=A31, A32=A32, dim=c(3,3)))
L = solve(diag(3*n)-A)
G = as.matrix(bdiag(bdiag(1/tau[1]*G1, 1/tau[2]*G2), 1/tau[3]*G3))
V = L%*%G%*%t(L)

orderindex = c(1,2,3)
models = permutations(n=3,r=3,v=orderindex,repeats.allowed=F)


Dinv_new <- function(Rho, n, cn, ns, udnei,q){
  Tau <- list()
  B <- list()
  invD <- list()
  for(i in 1:q){
    Tau[[i]] <- diag(n)
    B[[i]] <- matrix(0, n, n)
    for(j in 2:n){
      Tau[[i]][j,j] <- (1 + (ns[j-1] - 1) * Rho[i]^2) / (1 - Rho[i]^2)
      for(k in (udnei[(cn[j-1] + 1):(cn[j-1] + ns[j-1])])){
        B[[i]][j,k] <- Rho[i] / (1 + (ns[j-1] - 1) * Rho[i]^2)
      }
    }
    invD[[i]] <- t(diag(n) - B[[i]]) %*% Tau[[i]] %*% (diag(n) - B[[i]])
  }
  return(invD)
}

#Conduct simulations for different true model senarios: true_model = 1, 2, ..., 6
true_model = 1

mp_seed = list()
logml_seed = list()

##### data generation (50 datasets) #####
for(seed in 1:50){
  print(seed)
  set.seed(seed)

  W1=as.vector(rmvnorm(1,rep(0,q*n),V))
  reorder = order(rep(models[true_model,], each = 48))
  W = W1[reorder]
  
  X1 = cbind(rep(1, n), matrix(rnorm(n),ncol=1))
  X2 = cbind(rep(1, n), matrix(rnorm(2*n),ncol=2))
  X3 = cbind(rep(1, n), matrix(rnorm(2*n),ncol=2))

  X = as.matrix(bdiag(bdiag(X1, X2), X3))
  beta1 = c(1, 5)
  beta2 = c(2, 4, 5)
  beta3 = c(5, 3, 6)
  beta = c(beta1, beta2, beta3)
  
  sigma1 = sqrt(0.4)
  sigma2 = sqrt(0.4)
  sigma3 = sqrt(0.4)

  sigma = as.matrix(bdiag(bdiag((sigma1*diag(n)), (sigma2*diag(n))), (sigma3*diag(n))))
  Y = as.vector(X%*%beta+W+sigma%*%rnorm(n*q))
  Y1 = Y[1:n]
  Y2 = Y[(n+1):(2*n)]
  Y3 = Y[-(1:(2*n))]
  
  
  Y_list = list(Y1, Y2, Y3)
  X_list = list(X1, X2, X3)
  
  mcmc_list = list()
  # Results from six permutations
  for(i in 1:6){
    
    Y_list1 = Y_list[models[i,]]
    X_list1 = X_list[models[i,]]
    
    Y1 = Y_list1[[1]]
    Y2 = Y_list1[[2]]
    Y3 = Y_list1[[3]]
    
    X1 = X_list1[[1]]
    X2 = X_list1[[2]]
    X3 = X_list1[[3]]
    
    Y = c(Y1,Y2,Y3)
    X = as.matrix(bdiag(bdiag(X1, X2), X3))
    
    # Model for three diseases
    sink("DAGAR.txt")
    cat("
        model
        {
        Tau1[1, 1] <- 1
        Tau2[1, 1] <- 1
        Tau3[1, 1] <- 1
        
        for(j in 1:k){
        B1[1,j] <- 0
        B2[1,j] <- 0
        B3[1,j] <- 0
        }
        
        
        for(l in 1:(k-1)){
        for(h in (l+1):k){
        Tau1[l,h] <- 0
        Tau2[l,h] <- 0
        Tau3[l,h] <- 0
        }
        }
        
        
        for (i in 2:k)
        {
        Tau1[i,i] <- (1 + (ns[i-1] - 1) * rho1^2) / (1 - rho1^2)
        Tau2[i,i] <- (1 + (ns[i-1] - 1) * rho2^2) / (1 - rho2^2)
        Tau3[i,i] <- (1 + (ns[i-1] - 1) * rho3^2) / (1 - rho3^2)
        for(h in 1:(i-1)){
        Tau1[i,h] <- 0
        Tau2[i,h] <- 0
        Tau3[i,h] <- 0
        }
        b1[i] <- rho1 / (1 + (ns[i-1] - 1) * rho1^2)
        b2[i] <- rho2 / (1 + (ns[i-1] - 1) * rho2^2)
        b3[i] <- rho3 / (1 + (ns[i-1] - 1) * rho3^2)
        
        for(j in (udnei[(cn[i-1] + 1):(cn[i-1] + ns[i-1])])){
        B1[i,j] <- b1[i]
        B2[i,j] <- b2[i]
        B3[i,j] <- b3[i]
        }
        for(j in index1[((k)*(i-2)-cn[i-1]+1) : ((k)*(i-2)-cn[i-1] + (k - ns[i-1]))]){
        B1[i,j] <- 0
        B2[i,j] <- 0
        B3[i,j] <- 0
        }
        }
        
        Q1 <- t(I - B1) %*% Tau1 %*% (I - B1)
        Q2 <- t(I - B2) %*% Tau2 %*% (I - B2)
        Q3 <- t(I - B3) %*% Tau3 %*% (I - B3)
        
        W1[1:k] ~ dmnorm(rep(0, k), tau1*Q1)
        W2[1:k] ~ dmnorm(rep(0, k), tau2*Q2)
        W3[1:k] ~ dmnorm(rep(0, k), tau3*Q3)
        
        A21 <- eta021 * I + eta121 * Minc
        A31 <- eta031 * I + eta131 * Minc
        A32 <- eta032 * I + eta132 * Minc
        
        W[1:k] <- W1
        W[(k+1):(2*k)] <- A21 %*% W[1:k] + W2
        W[(2*k+1):(3*k)] <- A31 %*% W[1:k] + A32 %*% W[(k+1):(2*k)] + W3
        
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
        
        rho1 ~ dunif(0, 0.999)
        rho2 ~ dunif(0, 0.999)
        rho3 ~ dunif(0, 0.999)
        tau1 ~ dgamma(2, 0.1)
        tau2 ~ dgamma(2, 0.1)
        tau3 ~ dgamma(2, 0.1)
        eta021 ~ dnorm(0, 0.01)
        eta121 ~ dnorm(0, 0.01)
        eta031 ~ dnorm(0, 0.01)
        eta131 ~ dnorm(0, 0.01)
        eta032 ~ dnorm(0, 0.01)
        eta132 ~ dnorm(0, 0.01)
        
        taue1 ~ dgamma(2, 1)
        taue2 ~ dgamma(2, 1)
        taue3 ~ dgamma(2, 1)
        vare1 <- 1/taue1
        vare2 <- 1/taue2
        vare3 <- 1/taue3
        beta[1:8] ~ dmnorm(rep(0,8), (0.001*I1))
        }
        ", fill = TRUE)
    sink()
    
    model.data <- list(k = n, index1 = index1, I = diag(n), I1 = diag(8), Minc = Minc, ns = dni, cn = c(0, cni), udnei = udnei, X = X, Y = Y)
    model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, rho3 = 0.1, tau1 = 1, tau2 = 1, tau3 = 1, eta021 = 1, 
                                 eta121 = 1, eta031 = 1, eta131 = 1, eta032 = 1, eta132 = 1, taue1 = 1, taue2 = 1,
                                 taue3 = 1, beta = rep(0, 8), W1 = rep(0.1, n), 
                                 W2 = rep(0.1, n), W3 = rep(0.1, n))),2)
    model.param <- c("beta", "rho1", "rho2", "rho3", "tau1", "tau2", "tau3", "eta021", "eta121",
                     "eta031", "eta131", "eta032", "eta132", "vare1", "vare2", "vare3")
    set.seed(123)
    result1 <- jags(model.data, model.inits, model.param, "DAGAR.txt",
                    n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
    
    mcmc_list[[i]] = result1
  }
  
  Ylist = list()
  Xlist = list()
  
  for(i in 1:6){
    #print(i)
    Y_list1 = Y_list[models[i,]]
    X_list1 = X_list[models[i,]]
    
    Y1 = Y_list1[[1]]
    Y2 = Y_list1[[2]]
    Y3 = Y_list1[[3]]
    
    X1 = X_list1[[1]]
    X2 = X_list1[[2]]
    X3 = X_list1[[3]]
    
    Ylist[[i]] = c(Y1,Y2,Y3)
    Xlist[[i]] = as.matrix(bdiag(bdiag(X1, X2), X3))
  }
  
  #likelihood function
  Likli_fun = function(samples.row, data){
    A21_est = diag(samples.row["eta021"], data$n) + samples.row["eta121"] * data$Minc
    A31_est = diag(samples.row["eta031"], data$n) + samples.row["eta131"] * data$Minc
    A32_est = diag(samples.row["eta032"], data$n) + samples.row["eta132"] * data$Minc
    
    invD_est = Dinv_new(c(samples.row["rho1"],samples.row["rho2"],samples.row["rho3"]), data$n, data$cn, data$ns, data$udnei,q=3)
    G_est = as.matrix(bdiag(bdiag(1/samples.row["tau1"]*solve(invD_est[[1]]), 1/samples.row["tau2"]*solve(invD_est[[2]])), 
                            bdiag(1/samples.row["tau3"]*solve(invD_est[[3]]))))
    A_est = as.matrix(blockmatrix(names = c("0","A21","A31","0","0","A32","0","0","0"), 
                                  A21=A21, A31=A31, A32=A32, dim=c(3,3)))
    L_est = solve(diag(3*data$n)-A_est)
    V_est = L_est%*%G_est%*%t(L_est)
    Sigma = diag(c(rep(samples.row["vare1"],data$n),rep(samples.row["vare2"],data$n),rep(samples.row["vare3"],data$n)))
    Cov = V_est + Sigma
    Prec = solve(nearPD(Cov)$mat)
    beta = NULL
    for(par in 1:8){
      beta = c(beta, samples.row[paste("beta[", par, "]", sep="")])
    }
    mu = data$X%*%beta
    as.numeric(-2*data$n*log(2*pi) + 1/2*logdet(nearPD(Prec)$mat) - 1/2*t(data$Y-mu)%*%Prec%*%(data$Y-mu))
  }
  
  #unnormalized joint posterior function
  log_posterior <- function(samples.row, data) {
    
    eta_sum = dnorm(samples.row["eta021"], 0, 10, log=T) + dnorm(samples.row["eta121"], 0, 10, log=T) +
      dnorm(samples.row["eta031"], 0, 10, log=T) + dnorm(samples.row["eta131"], 0, 10, log=T) +
      dnorm(samples.row["eta032"], 0, 10, log=T) + dnorm(samples.row["eta132"], 0, 10, log=T)
    rho_sum = dunif(samples.row["rho1"], 0, 0.999, log=T) + dunif(samples.row["rho2"], 0, 0.999, log=T) +
      dunif(samples.row["rho3"], 0, 0.999, log=T) 
    tau_sum = dgamma(samples.row["tau1"], 2, rate = 0.1, log=T) + dgamma(samples.row["tau2"], 2, rate = 0.1, log=T) + 
      dgamma(samples.row["tau3"], 2, rate = 0.1, log=T) 
    sigma_sum = dinvgamma(samples.row["vare1"], 2, 1, log=T) + dinvgamma(samples.row["vare2"], 2, 1, log=T) +
      dinvgamma(samples.row["vare3"], 2, 1, log=T)
    
    beta = NULL
    for(par in 1:8){
      beta = c(beta, samples.row[paste("beta[", par, "]", sep="")])
    }
    beta_sum = dmnorm(beta, 0, diag(sqrt(1000),8), log=T)
    as.numeric(Likli_fun(samples.row, data) + eta_sum + rho_sum + tau_sum + sigma_sum + beta_sum)
  }
  
  #upper and lower bound for parameters
  beta_name = NULL
  for(par in 1:8){
    beta_name = c(beta_name, paste("beta[", par, "]", sep=""))
  }
  
  cname <- c(beta_name, "eta021", "eta121", "eta031", "eta131", "eta032", "eta132", 
             "rho1", "rho2", "rho3", "tau1", "tau2", "tau3",  "vare1", "vare2", "vare3")
  lb <- rep(-Inf, length(cname))
  ub <- rep(Inf, length(cname))
  names(lb) <- names(ub) <- cname
  lb[15:23] <- 0
  ub[15:17] <- 0.999
  
  #bridge sampling for log marginal likelihood
  model.bridge = list()
  for(model in 1:6){
    print(model)
    
    model.bridge[[model]] <- bridge_sampler(samples = mcmc_list[[model]], 
                                            data = list(Y = Ylist[[model]], X = Xlist[[model]], Minc = Minc, n = n, ns = dni, cn = c(0, cni), udnei = udnei),
                                            log_posterior = log_posterior, lb = lb,
                                            ub = ub, maxiter = 50000, silent = TRUE)
  }
  
  logml = do.call(c, lapply(model.bridge, function(x) x$logml))
  mp = exp(logml)/sum(exp(logml))
  print(mp)
  
  mp_seed[[seed]] = mp
  logml_seed[[seed]] = model.bridge
}

saveRDS(mp_seed, paste("mp_model", true_model, ".rds", sep = ""))
saveRDS(logml_seed, paste("logml_model", true_model, ".rds", sep = ""))
