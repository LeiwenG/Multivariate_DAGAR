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
library(mvtnorm)
library(invgamma)
library(MASS)
library(rjags)
library(R2jags)
library(LaplacesDemon)

inclattice=function(m){
  n=m^2
  Minc=matrix(0,n,n)
  for(i in 1:(m-1))	for(j in 1:(m-1)) Minc[(i-1)*m+j,(i-1)*m+j+1]=Minc[(i-1)*m+j,i*m+j]=1
  for(i in 1:(m-1)) Minc[(i-1)*m+m,i*m+m]=1
  for(j in 1:(m-1)) Minc[(m-1)*m+j,(m-1)*m+j+1]=1
  Minc+t(Minc)
}

mysummary = function(invector) {
  c(mean(invector), median(invector), sd(invector), quantile(invector, .025), quantile(invector,.975))
}

correlation <- function(D1, n, ni_wo, udnei_wo, cni_wo){
  Z = rep(0, n)
  
  for(j in 1:ni_wo[1]){
    Z[1] = Z[1] + D1[udnei_wo[j],1]
  }
  
  for(i in 2:n){
    for(j in 1:ni_wo[i]){
      Z[i] = Z[i] + D1[udnei_wo[(cni_wo[i-1] + j)],i]
    }
  }
  
  Z21 = as.matrix(cbind(diag(D1), Z))
  return(Z21)
}
#N = 100
q = 2

bhat = matrix(0, nrow = 5, ncol = 100)
blower = matrix(0, nrow = 5, ncol = 100)
bupper = matrix(0, nrow = 5, ncol = 100)
varehat = matrix(0, nrow = 2, ncol = 100)
varelower = matrix(0, nrow = 2, ncol = 100)
vareupper = matrix(0, nrow = 2, ncol = 100)
tauwhat = matrix(0, nrow = 2, ncol = 100)
tauwlower = matrix(0, nrow = 2, ncol = 100)
tauwupper = matrix(0, nrow = 2, ncol = 100)
rhohat = matrix(0, nrow = 2, ncol = 100)
rholower = matrix(0, nrow = 2, ncol = 100)
rhoupper = matrix(0, nrow = 2, ncol = 100)
etahat = matrix(0, nrow = 2, ncol = 100)
etalower = matrix(0, nrow = 2, ncol = 100)
etaupper = matrix(0, nrow = 2, ncol = 100)

w1hat = matrix(0, nrow = 100, ncol = 48)
w1lower = matrix(0, nrow = 100, ncol = 48)
w1upper = matrix(0, nrow = 100, ncol = 48)
w2hat = matrix(0, nrow = 100, ncol = 48)
w2lower = matrix(0, nrow = 100, ncol = 48)
w2upper = matrix(0, nrow = 100, ncol = 48)

corhat = matrix(0, nrow = 100, ncol = 48)
corlower = matrix(0, nrow = 100, ncol = 48)
corupper = matrix(0, nrow = 100, ncol = 48)

Covhat = matrix(0, nrow = 100, ncol = 48)
Covlower = matrix(0, nrow = 100, ncol = 48)
Covupper = matrix(0, nrow = 100, ncol = 48)

W1_path = matrix(0, nrow = 100, ncol = 48)
W2_path = matrix(0, nrow = 100, ncol = 48)

WAIC = rep(0,100)
g1 = rep(0,100)
D1 = rep(0,100)
P1 = rep(0,100)
g2 = rep(0,100)
D2 = rep(0,100)
P2 = rep(0,100)

KL = matrix(0, nrow = 100, ncol = 15000)

graphnum=3

# Assign true values for parameters
rho = c(0.2, 0.8)
tau = c(0.25, 0.25)
beta1 = c(1, 5)
beta2 = c(2, 4, 5)
beta = c(beta1, beta2)
sigma1 = sqrt(0.4)
sigma2 = sqrt(0.4)

# Set up for low correlation between diseases
eta0_21 = 0.05
eta1_21 = 0.1

# Change values to conduct simulation for different correlation scenarios among diseases
# low correlation: eta0_21 = 0.05, eta1_21 = 0.1
# medium correlation: eta0_21 = 0.5, eta1_21 = 0.3
# high correlation: eta0_21 = 2.5, eta1_21 = 0.5



for(seed in 1:100){
  print(seed)
  #### setting up Minc, dmat and related variables #####
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
  #colnames(Minc)[perm]
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
  
  G1 = rho[1]^dmat
  G2 = rho[2]^dmat
  
  A = diag(eta0_21, n) + eta1_21 * Minc
  L = as.matrix(blockmatrix(names = c("I1","A","0","I2"), I1 = diag((q-1)*n), 
                            A=A, I2 = diag(n), dim=c(2,2)))
  G = as.matrix(bdiag(1/tau[1]*G1, 1/tau[2]*G2))
  V = L%*%G%*%t(L)
  
  cn = c(0, cni)
  ns = dni
  #Covariance of two diseases in the same region
  Cov_region = 1/tau[1]*correlation(G1, n, ni_wo, udnei_wo, cni_wo)%*%c(eta0_21,eta1_21)
  v1_region = diag(1/tau[1]*G1)
  v2_region = diag(1/tau[1]*A%*%G1%*%t(A) + 1/tau[2]*G2)
  cor_region = Cov_region/(sqrt(v1_region*v2_region))
  
  
  ##### data generation #####
  
  set.seed(seed)
  W=as.vector(rmvnorm(1,rep(0,q*n),V))
  W1 = W[1:n]
  W2 = W[-(1:n)]
  X1 = cbind(rep(1, n), matrix(rnorm(n),ncol=1))
  X2 = cbind(rep(1, n), matrix(rnorm(2*n),ncol=2))
  
  X = as.matrix(bdiag(X1, X2))
  
  
  W1_path[seed, ] = W1
  W2_path[seed, ] = W2
  
  sigma = as.matrix(bdiag((sigma1*diag(n)), (sigma2*diag(n))))
  Y = as.vector(X%*%beta+W+sigma%*%rnorm(n*q))
  Y1 = Y[1:n]
  Y2 = Y[-(1:n)]
  
  
  D = diag(rowSums(Minc))
  
  tod1=Sys.time()
  sink("GMCAR.txt")
  cat("
      model
      {
      Q1[1:k, 1:k] <- tau1*(D - rho1*Minc)
      Q2[1:k, 1:k] <- tau2*(D - rho2*Minc)
      W1[1:k] ~ dmnorm(rep(0, k), Q1)
      W2[1:k] ~ dmnorm(rep(0, k), Q2)
      A <- eta0 * I + eta1 * Minc
      W[1:k] <- W1
      W[(k+1):(2*k)] <- A %*% W1 + W2
      
      for (i in 1:k)
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue1)
      }
      for (i in (k+1):(2*k))
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue2)
      }
      
      rho1 ~ dunif(0, 0.999)
      rho2 ~ dunif(0, 0.999)
      tau1 ~ dgamma(2, 8)
      tau2 ~ dgamma(2, 8)
      eta0 ~ dnorm(0, 0.01)
      eta1 ~ dnorm(0, 0.01)
      taue1 ~ dgamma(2, 0.4)
      taue2 ~ dgamma(2, 0.4)
      vare1 <- 1/taue1
      vare2 <- 1/taue2
      beta[1:5] ~ dmnorm(rep(0, 5), (0.001*I1))
      }
      ", fill = TRUE)
  sink()
  
  model.data <- list(k = n, I = diag(n), I1 = diag(5), Minc = Minc, D = diag(rowSums(Minc)), X = X, Y = Y)
  model.inits <- list(list(rho1 = 0.1, rho2 = 0.1, tau1 = 1, tau2 = 1, eta0 = 1, 
                           eta1 = 1, taue1 = 1, taue2 = 1, beta = rep(0, 5), W1 = rep(0.1, n), 
                           W2 = rep(0.1, n)))
  model.param <- c("beta", "rho1", "rho2", "tau1", "tau2", "eta0", "eta1", "vare1", "vare2", "W")
  set.seed(seed)
  estimate <- jags(model.data, model.inits, model.param, "GMCAR.txt",
                   n.chains = 1, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  print(estimate)
  tod2=Sys.time()
  
  
  result1 = data.frame(estimate$BUGSoutput$sims.array[,1,-c(1:96, 102)])
  colnames(result1) <- c("beta1", "beta2","beta3","beta4","beta5","eta0_21","eta1_21",
                         "rho1", "rho2", "tausq1", "tausq2", "sigmasq1", "sigmasq2")
  
  result = list()
  result[[1]] = result1
  result[[2]] = data.frame(estimate$BUGSoutput$sims.array[,1,c(1:96)])
  
  sample.mcmc = as.matrix(result1[, c(1:5, 12, 13, 10, 11, 8, 9, 6, 7)])
  estimate <- t(apply(sample.mcmc, 2, mysummary))
  
  #compute correlation between random effects of diseases in same region
  CovD = matrix(0, 1000, n)
  v1 = matrix(0, 1000, n)
  v2 = matrix(0, 1000, n)
  cor = matrix(0, 1000, n)
  for(tem in 1:1000){
    print(tem)
    invD1 = D - result[[1]]$rho1[14000+tem]*Minc
    invD2 = D - result[[1]]$rho2[14000+tem]*Minc
    D1 = solve(invD1)
    D2 = solve(invD2)
    tausq1_est = result[[1]]$tausq1[14000+tem]
    tausq2_est = result[[1]]$tausq2[14000+tem]
    eta_est = c(result[[1]]$eta0_21[14000+tem], result[[1]]$eta1_21[14000+tem])
    A21_est = eta_est[1]*diag(n) + eta_est[2]*Minc
    term = D1%*%t(A21_est)
    CovD[tem,] = diag(1/tausq1_est*term)
    v1[tem,] = diag(1/tausq1_est*D1)
    v2[tem,] = diag(1/tausq1_est*A21_est%*%term + 1/tausq2_est*D2)
    cor[tem,] = CovD[tem,] / (sqrt(v1[tem,]) * sqrt(v2[tem,]))
   }
  
  bhat[, seed] = estimate[1:5, 1]
  blower[, seed] = estimate[1:5, 4]
  bupper[, seed] = estimate[1:5, 5]
  varehat[, seed] = estimate[6:7, 1]
  varelower[, seed] = estimate[6:7, 4]
  vareupper[, seed] = estimate[6:7, 5]
  tauwhat[, seed] = estimate[8:9, 1]
  tauwlower[, seed] = estimate[8:9, 4]
  tauwupper[, seed] = estimate[8:9, 5]
  rhohat[, seed] = estimate[10:11, 1]
  rholower[, seed] = estimate[10:11, 4]
  rhoupper[, seed] = estimate[10:11, 5]
  etahat[, seed] = estimate[12:13, 1]
  etalower[, seed] = estimate[12:13, 4]
  etaupper[, seed] = estimate[12:13, 5]
  
  w1hat[seed,] = apply(result[[2]], 2, mean)
  w1lower[seed,] = apply(result[[2]], 2, quantile, 0.025)
  w1upper[seed,] = apply(result[[2]], 2, quantile, 0.975)
  
  Covhat[seed,] = apply(CovD, 2, mean)
  Covlower[seed,] = apply(CovD, 2, quantile, 0.025)
  Covupper[seed,] = apply(CovD, 2, quantile, 0.975)
  
  #correlation between diseases
  corhat[seed,] = apply(cor, 2, mean)
  corlower[seed,] = apply(cor, 2, quantile, 0.025)
  corupper[seed,] = apply(cor, 2, quantile, 0.975)
  
  W1_mcmc = result[[2]][,1:48]
  W2_mcmc = result[[2]][,49:96]
  
  #WAIC
  PL_single <- matrix(0, nrow = 2*n, ncol = 15000)
  for(i in 1:15000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:5]
    sigmasq1 <- theta[6]
    sigmasq2 <- theta[7]
    
    for(j in 1:n){
      v1 <- sigmasq1
      PL_single[j,i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v1) - 1/(2*v1)*(Y1[j]-X1[j,]%*%as.vector(beta1) - W1_mcmc[i,j])^2)
      v2 <- sigmasq2
      PL_single[(j+n),i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v2) - 1/(2*v2)*(Y2[j]-X2[j,]%*%as.vector(beta2) - W2_mcmc[i,j])^2)
    }
    
  }
  
  WAIC[seed] <- WAIC(PL_single)
  
  #D score
  Y_rep1 <- matrix(0, nrow = n, ncol = 15000)
  Y_rep2 <- matrix(0, nrow = n, ncol = 15000)
  for(i in 1:15000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:5]
    sigmasq1 <- theta[6]
    sigmasq2 <- theta[7]
    set.seed(seed)
    z1 <- rnorm(48, 0, 1)
    z2 <- rnorm(48, 0, 1)
    Y_rep1[,i] <- X1 %*% t(as.vector(beta1)) +  result[[2]]$W1[i,] + sqrt(as.numeric(sigmasq1)) * z1
    Y_rep2[,i] <- X2 %*% t(as.vector(beta2)) +  result[[2]]$W2[i,] + sqrt(as.numeric(sigmasq2)) * z2
    
  }
  mu_rep1 = rowMeans(Y_rep1)
  mu_rep2 = rowMeans(Y_rep2)
  
  var_rep1 = rowVar(Y_rep1)
  var_rep2 = rowVar(Y_rep2)
  
  G.latent1 = sum((Y1 - mu_rep1)^2)
  P.latent1 = sum(var_rep1)
  D.latent1 = G.latent1 + P.latent1
  
  G.latent2 = sum((Y2 - mu_rep2)^2)
  P.latent2 = sum(var_rep2)
  D.latent2 = G.latent2 + P.latent2
  
  g1[seed] = G.latent1
  P1[seed] = P.latent1
  D1[seed] = D.latent1
  
  g2[seed] = G.latent2
  P2[seed] = P.latent2
  D2[seed] = D.latent2
  
  #Kullback–Leibler divergence
  KL_car = NULL
  for(i in 1:nrow(estimatec)){
    rho_est = estimatec[i, 10:11]
    invD1 = diag(rowSums(Minc)) - rho_est[1]*Minc
    invD2 = diag(rowSums(Minc)) - rho_est[2]*Minc
    D1 = solve(invD1)
    D2 = solve(invD2)
    tausq1_est = estimatec[i,8]
    tausq2_est = estimatec[i,9]
    eta_est = estimatec[i,12:13]
    A21_est = eta_est[1]*diag(n) + eta_est[2]*Minc
    L_est = as.matrix(blockmatrix(names = c("I1","A","0","I2"), I1 = diag(n), 
                                  A=A21_est, I2 = diag(n), dim=c(2,2)))
    
    G_est = as.matrix(bdiag(1/tausq1_est*D1, 1/tausq2_est*D2))
    V_est = L_est%*%G_est%*%t(L_est)
    
    KL_car[i] = kl.norm(mu1 = rep(0, 2*n), S1= V, mu2 = rep(0, 2*n), S2= V_est)
  }
  
  KL[seed,] = KL_car
  
}

modelgraph = "usa"
saveRDS(bhat, paste("bhat_", modelgraph, ".rds", sep=""))
saveRDS(blower, paste("blower_", modelgraph, ".rds", sep=""))
saveRDS(bupper, paste("bupper_", modelgraph, ".rds", sep=""))

saveRDS(varehat, paste("varehat_", modelgraph, ".rds", sep=""))
saveRDS(varelower, paste("varelower_", modelgraph, ".rds", sep=""))
saveRDS(vareupper, paste("vareupper_", modelgraph, ".rds", sep=""))

saveRDS(tauwhat, paste("tauwhat_", modelgraph, ".rds", sep=""))
saveRDS(tauwlower, paste("tauwlower_", modelgraph, ".rds", sep=""))
saveRDS(tauwupper, paste("tauwupper_", modelgraph, ".rds", sep=""))

saveRDS(rhohat, paste("rhohat_", modelgraph, ".rds", sep=""))
saveRDS(rholower, paste("rholower_", modelgraph, ".rds", sep=""))
saveRDS(rhoupper, paste("rhoupper_", modelgraph, ".rds", sep=""))

saveRDS(etahat, paste("etahat_", modelgraph, ".rds", sep=""))
saveRDS(etalower, paste("etalower_", modelgraph, ".rds", sep=""))
saveRDS(etaupper, paste("etaupper_", modelgraph, ".rds", sep=""))

saveRDS(w1hat, paste("w1hat_", modelgraph, ".rds", sep=""))
saveRDS(w1lower, paste("w1lower_", modelgraph, ".rds", sep=""))
saveRDS(w1upper, paste("w1upper_", modelgraph, ".rds", sep=""))
saveRDS(w2hat, paste("w2hat_", modelgraph, ".rds", sep=""))
saveRDS(w2lower, paste("w2lower_", modelgraph, ".rds", sep=""))
saveRDS(w2upper, paste("w2upper_", modelgraph, ".rds", sep=""))

saveRDS(W1_path, paste("W1_", modelgraph, ".rds", sep=""))
saveRDS(W2_path, paste("W2_", modelgraph, ".rds", sep=""))

saveRDS(corhat, paste("corhat_", modelgraph, ".rds", sep=""))
saveRDS(corlower, paste("corlower_", modelgraph, ".rds", sep=""))
saveRDS(corupper, paste("corupper_", modelgraph, ".rds", sep=""))

saveRDS(Covhat, paste("Covhat_", modelgraph, ".rds", sep=""))
saveRDS(Covlower, paste("Covlower_", modelgraph, ".rds", sep=""))
saveRDS(Covupper, paste("Covupper_", modelgraph, ".rds", sep=""))

saveRDS(cor_region, "cor_region.usa.rds")
saveRDS(Cov_region, "Cov_region.usa.rds")

saveRDS(WAIC, paste("WAIC_", modelgraph, ".rds", sep=""))
saveRDS(KL, paste("KL_", modelgraph, ".rds", sep=""))

saveRDS(g1, paste("G1_", modelgraph, ".rds", sep=""))
saveRDS(D1, paste("D1_", modelgraph, ".rds", sep=""))
saveRDS(P1, paste("P1_", modelgraph, ".rds", sep=""))

saveRDS(g2, paste("G2_", modelgraph, ".rds", sep=""))
saveRDS(D2, paste("D2_", modelgraph, ".rds", sep=""))
saveRDS(P2, paste("P2_", modelgraph, ".rds", sep=""))