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
library(monomvn)

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



#N = 100
q = 2

bhat = matrix(0, nrow = 5, ncol = 100)
blower = matrix(0, nrow = 5, ncol = 100)
bupper = matrix(0, nrow = 5, ncol = 100)
varehat = matrix(0, nrow = 2, ncol = 100)
varelower = matrix(0, nrow = 2, ncol = 100)
vareupper = matrix(0, nrow = 2, ncol = 100)
rhohat = matrix(0, nrow = 2, ncol = 100)
rholower = matrix(0, nrow = 2, ncol = 100)
rhoupper = matrix(0, nrow = 2, ncol = 100)
Ahat = matrix(0, nrow = 3, ncol = 100)
Alower = matrix(0, nrow = 3, ncol = 100)
Aupper = matrix(0, nrow = 3, ncol = 100)

w1hat = matrix(0, nrow = 100, ncol = 96)
w1lower = matrix(0, nrow = 100, ncol = 96)
w1upper = matrix(0, nrow = 100, ncol = 96)
w2hat = matrix(0, nrow = 100, ncol = 96)
w2lower = matrix(0, nrow = 100, ncol = 96)
w2upper = matrix(0, nrow = 100, ncol = 96)


W1_path = matrix(0, nrow = 100, ncol = 48)
W2_path = matrix(0, nrow = 100, ncol = 48)

WAIC = rep(0,100)
D_score = rep(0,100)
KL = matrix(0, nrow = 100, ncol = 15000)

graphnum=3
for(seed in 51:100){
  print(seed)
  #### setting up Minc, dmat and related variables #####
  if(graphnum==1){
    graph="path"
    n=100
    Minc=matrix(0,n,n)
    for(i in 1:(n-1)) Minc[i,i+1]=Minc[i+1,i]=1
    dmat=as.matrix(dist(1:n))
  }else{
    if(graphnum==2){
      graph="grid"
      m=10
      Minc=inclattice(m)
      s=cbind(rep(1:m,m),kronecker(1:m,rep(1,m)))
      dmat=as.matrix(dist(s))
    }else{
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
    }
  }
  
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
  #rho1 = 0.2
  #rho2 = 0.9
  #tau1 = 10
  #tau2 = 10
  #Rho = c(0.2, 0.9)
  rho = c(0.2, 0.8)
  a11 = 1
  a21 = 0.7
  a22 = 1
  tau = c(0.25, 0.25)
  #tausq1 = 10
  #tausq2 = 10
  
  
  G1 = solve(diag(rowSums(Minc)) - rho[1] * Minc)
  G2 = solve(diag(rowSums(Minc)) - rho[2] * Minc)
  A = matrix(c(a11, 0, a21, a22), nrow = 2, byrow = TRUE)
  G = as.matrix(bdiag(G1, G2))
  prod = kronecker(A, diag(n))
  V = prod %*% G %*% t(prod)
  
  cn = c(0, cni)
  ns = dni
  #Covariance of two diseases in the same region
  Cov_region = diag(V[1:48, 49:96])
  v1_region = diag(V[1:48, 1:48])
  v2_region = diag(V[49:96, 49:96])
  #correlation: 0.5
  cor_region = Cov_region/(sqrt(v1_region*v2_region))
  
  #muw <- rep(0,n)
  #id <- diag(n)
  
  ##### data generation #####
  #sigs=4
  set.seed(seed)
  W=as.vector(rmvnorm(1,rep(0,q*n),V))
  W1 = W[1:n]
  W2 = W[-(1:n)]
  X1 = cbind(rep(1, n), matrix(rnorm(n),ncol=1))
  X2 = cbind(rep(1, n), matrix(rnorm(2*n),ncol=2))
  #X1 = matrix(rnorm(2*n),ncol=2)
  #X2 = matrix(rnorm(3*n),ncol=3)
  X = as.matrix(bdiag(X1, X2))
  beta1 = c(1, 5)
  beta2 = c(2, 4, 5)
  beta = c(beta1, beta2)
  #Beta1 = c(1, 5)
  #Beta2 = c(2, 4, 5)
  
  W1_path[seed, ] = W1
  W2_path[seed, ] = W2
  
  sigma1 = sqrt(0.4)
  sigma2 = sqrt(0.4)
  #sigmasq1 = 0.01
  #sigmasq2 = 0.01
  sigma = as.matrix(bdiag((sigma1*diag(n)), (sigma2*diag(n))))
  Y = as.vector(X%*%beta+W+sigma%*%rnorm(n*q))
  Y1 = Y[1:n]
  Y2 = Y[-(1:n)]
  
  D = diag(rowSums(Minc))
  
  tod1=Sys.time()
  sink("JCAR.txt")
  cat("
      model
      {
        Q1[1:k, 1:k] <- D - rho1*Minc
        Q2[1:k, 1:k] <- D - rho2*Minc
        
        f1[1:k] ~ dmnorm(rep(0, k), Q1)
        f2[1:k] ~ dmnorm(rep(0, k), Q2)
        
        W[1:k] <- A[1,1]*f1 
        W[(k+1):(2*k)] <- A[2,1] * f1 + A[2,2] * f2
      
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
      
      for(r in 1:2){
            A[r, r] ~ dlnorm(0, 16)
      }
        
      A[2, 1] ~ dnorm(0, 0.01)
      A[1, 2] <- 0
      
      rho1 ~ dunif(0, 0.999)
      rho2 ~ dunif(0, 0.999)
      
      taue1 ~ dgamma(2, 0.4)
      taue2 ~ dgamma(2, 0.4)
      vare1 <- 1/taue1
      vare2 <- 1/taue2
      beta[1:5] ~ dmnorm(rep(0, 5), (0.001*I1))
      }
      ", fill = TRUE)
  sink()
  
  model.data <- list(k = n, I1 = diag(5), 
                     Minc = Minc, D = diag(rowSums(Minc)), X = X, Y = Y)
  model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, taue1 = 1, taue2 = 1, beta = rep(0, 5))), 2)
  model.param <- c("A", "beta", "rho1", "rho2", "vare1", "vare2", "W")
  set.seed(seed)
  estimate <- jags(model.data, model.inits, model.param, "JCAR.txt",
                   n.chains = 2, n.iter = 30000, n.burnin = 15000, n.thin = 1)
  print(estimate)
  tod2=Sys.time()
  
  result1 = data.frame(estimate$BUGSoutput$sims.array[,1,-c(3, 5:100, 106)])
  colnames(result1) <- c("a11", "a21", "a22", "beta1", "beta2","beta3","beta4","beta5",
                         "rho1", "rho2", "sigmasq1", "sigmasq2")
  
  result = list()
  result[[1]] = result1
  result[[2]] = data.frame(estimate$BUGSoutput$sims.array[,1,c(5:100)])
  
  sample.mcmc = as.matrix(result1[, c(4:8, 11, 12, 9, 10, 1:3)])
  estimate <- t(apply(sample.mcmc, 2, mysummary))
  
  bhat[, seed] = estimate[1:5, 1]
  blower[, seed] = estimate[1:5, 4]
  bupper[, seed] = estimate[1:5, 5]
  varehat[, seed] = estimate[6:7, 1]
  varelower[, seed] = estimate[6:7, 4]
  vareupper[, seed] = estimate[6:7, 5]
  rhohat[, seed] = estimate[8:9, 1]
  rholower[, seed] = estimate[8:9, 4]
  rhoupper[, seed] = estimate[8:9, 5]
  Ahat[, seed] = estimate[10:12, 1]
  Alower[, seed] = estimate[10:12, 4]
  Aupper[, seed] = estimate[10:12, 5]
  
  w1hat[seed,] = apply(result[[2]][,1:48], 2, mean)
  w1lower[seed,] = apply(result[[2]][,1:48], 2, quantile, 0.025)
  w1upper[seed,] = apply(result[[2]][,1:48], 2, quantile, 0.975)
  w2hat[seed,] = apply(result[[2]][,49:96], 2, mean)
  w2lower[seed,] = apply(result[[2]][,49:96], 2, quantile, 0.025)
  w2upper[seed,] = apply(result[[2]][,49:96], 2, quantile, 0.975)
  
  
  W1_mcmc = result[[2]][,1:48]
  W2_mcmc = result[[2]][,49:96]
  
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
  
  WAIC[seed] <- WAIC(PL_single)$WAIC
  
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
    Y_rep1[,i] <- as.numeric(X1 %*% as.vector(beta1) +  W1_mcmc[i,] + sqrt(as.numeric(sigmasq1)) * z1)
    Y_rep2[,i] <- as.numeric(X2 %*% as.vector(beta2) +  W2_mcmc[i,] + sqrt(as.numeric(sigmasq2)) * z2)
    
  }
  mu_rep1 = rowMeans(Y_rep1)
  mu_rep2 = rowMeans(Y_rep2)
  
  var_rep1 = rowVars(Y_rep1)
  var_rep2 = rowVars(Y_rep2)
  
  G.latent1 = sum((Y1 - mu_rep1)^2)
  P.latent1 = sum(var_rep1)
  D.latent1 = G.latent1 + P.latent1
  
  G.latent2 = sum((Y2 - mu_rep2)^2)
  P.latent2 = sum(var_rep2)
  D.latent2 = G.latent2 + P.latent2
  
  D_score[seed] = D.latent1 + D.latent2
  
  #Kullback–Leibler divergence
  KL_dagar = NULL
  for(i in 1:15000){
    #print(i)
    theta <- as.vector(sample.mcmc[i,])
    beta1 <- as.vector(theta[1:2])
    beta2 <- as.vector(theta[3:5])
    sigmasq1 <- as.numeric(theta[6])
    sigmasq2 <- as.numeric(theta[7])
    
    KL_dagar[i] = kl.norm(mu1 = X%*%beta+W, S1 = sigma^2, mu2 = as.vector(X%*%unlist(c(beta1, beta2))) + 
                            c(as.numeric(W1_mcmc[i,]), as.numeric(W2_mcmc[i,])), 
                          S2= as.matrix(bdiag(sigmasq1*diag(n), sigmasq2*diag(n))))
  }
  
  KL[seed,] = KL_dagar
  
  modelgraph = "usa"
  num = 1
  
  saveRDS(bhat, paste("bhat_", modelgraph, num, ".rds", sep=""))
  saveRDS(blower, paste("blower_", modelgraph, num, ".rds", sep=""))
  saveRDS(bupper, paste("bupper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(varehat, paste("varehat_", modelgraph, num, ".rds", sep=""))
  saveRDS(varelower, paste("varelower_", modelgraph, num, ".rds", sep=""))
  saveRDS(vareupper, paste("vareupper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(rhohat, paste("rhohat_", modelgraph, num, ".rds", sep=""))
  saveRDS(rholower, paste("rholower_", modelgraph, num, ".rds", sep=""))
  saveRDS(rhoupper, paste("rhoupper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(w1hat, paste("w1hat_", modelgraph, num, ".rds", sep=""))
  saveRDS(w1lower, paste("w1lower_", modelgraph, num, ".rds", sep=""))
  saveRDS(w1upper, paste("w1upper_", modelgraph, num, ".rds", sep=""))
  saveRDS(w2hat, paste("w2hat_", modelgraph, num, ".rds", sep=""))
  saveRDS(w2lower, paste("w2lower_", modelgraph, num, ".rds", sep=""))
  saveRDS(w2upper, paste("w2upper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(W1_path, paste("W1_", modelgraph, num, ".rds", sep=""))
  saveRDS(W2_path, paste("W2_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(cor_region, "cor_region.usa.rds")
  saveRDS(Cov_region, "Cov_region.usa.rds")
  
  saveRDS(WAIC, paste("WAIC_", modelgraph, num, ".rds", sep=""))
  saveRDS(KL, paste("KL_", modelgraph, num, ".rds", sep=""))
  saveRDS(D_score, paste("D2_", modelgraph, num, ".rds", sep=""))
}
