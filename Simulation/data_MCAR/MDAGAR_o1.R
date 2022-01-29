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

beta_new <- function(Y1, Y2, X1, X2, W1, W2, sigmasq1, sigmasq2, sigmasq_beta){
  M1 <- solve(1/sigmasq1*t(X1)%*%X1 + 1/sigmasq_beta*diag(ncol(X1)))
  M2 <- solve(1/sigmasq2*t(X2)%*%X2 + 1/sigmasq_beta*diag(ncol(X2)))
  m1 <- 1/sigmasq1*t(X1)%*%(Y1-W1)
  m2 <- 1/sigmasq2*t(X2)%*%(Y2-W2)
  mu1 <- M1%*%m1
  mu2 <- M2%*%m2
  beta_1 <- rmvnorm(1, mean = mu1, sigma = M1)
  beta_2 <- rmvnorm(1, mean = mu2, sigma = M2)
  return(c(beta_1, beta_2))
}


Dinv_new <- function(Rho, n, cn, ns, udnei){
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

quad_new <- function(Rho, n, cn, ns, udnei, W1, W2){
  Tau1 <- NULL
  Tau2 <- NULL
  q1 <- NULL
  q2 <- NULL
  
  Tau1[1] <- 1
  Tau2[1] <- 1
  
  for(j in 2:n){
    Tau1[j] <- (1 + (ns[j-1] - 1) * Rho[1]^2) / (1 - Rho[1]^2)
    Tau2[j] <- (1 + (ns[j-1] - 1) * Rho[2]^2) / (1 - Rho[2]^2)
    #for(k in (udnei[(cn[j-1] + 1):(cn[j-1] + ns[j-1])])){
    #B1[k] <- Rho[1] / (1 + (ns[j-1] - 1) * Rho[1]^2)
    #B2[k] <- Rho[2] / (1 + (ns[j-1] - 1) * Rho[2]^2)
    b1 <- Rho[1] / (1 + (ns[j-1] - 1) * Rho[1]^2)
    b2 <- Rho[2] / (1 + (ns[j-1] - 1) * Rho[2]^2)
    q1[j-1] <- Tau1[j] * (W1[j] - b1*sum(W1[udnei[(cn[j-1] + 1):(cn[j-1] + ns[j-1])]]))^2
    q2[j-1] <- Tau2[j] * (W2[j] - b2*sum(W2[udnei[(cn[j-1] + 1):(cn[j-1] + ns[j-1])]]))^2
    
  }
  
  qf1 <- Tau1[1]*W1[1]^2 + sum(q1)
  qf2 <- Tau2[1]*W2[1]^2 + sum(q2)
  
  det1 <- sum(log(Tau1))
  det2 <- sum(log(Tau2))
  
  return(c(qf1, qf2, det1, det2))
}

W1_new <- function(Y1, X1, eta0_21, eta1_21, beta1, W2, sigmasq1, tausq1, tausq2, invD){
  A21 <- eta0_21 * diag(n) + eta1_21 * Minc
  DA <- invD[[2]]%*%A21
  M1 <- solve(tausq1*invD[[1]] + tausq2*t(A21)%*%DA + diag(1/sigmasq1, n))
  m1 <- tausq2*t(DA)%*%W2 + 1/sigmasq1*(Y1-X1%*%beta1)
  mu1 <- M1%*%m1
  return(rmvnorm(1, mean = mu1, sigma = M1))
}

W2_new <- function(Y2, X2, eta0_21, eta1_21, beta2, W1, sigmasq2, tausq2, invD){
  A21 <- eta0_21 * diag(n) + eta1_21 * Minc
  M2 <- solve(tausq2*invD[[2]] + diag(1/sigmasq2, n))
  m2 <- tausq2*invD[[2]]%*%(A21%*%as.vector(W1)) + 1/sigmasq2*(Y2-X2%*%beta2)
  mu2 <- M2%*%m2
  return(rmvnorm(1, mean = mu2, sigma = M2))
}

sigmasq_new <- function(Y1, Y2, X1, X2, W1, W2, beta1, beta2, a1, b1, n){
  #a1 = 1
  #b1=0.1
  s <- a1 + n/2
  term1 <- Y1-X1%*%beta1-as.vector(W1)
  term2 <- Y2-X2%*%beta2-as.vector(W2)
  r1 <- b1 + 1/2*t(term1)%*%(term1)
  r2 <- b1 + 1/2*t(term2)%*%(term2)
  return(c(rinvgamma(1, shape = s, r1),rinvgamma(1, shape = s, r2)))
}

tausq_new <- function(Rho, eta_21, Z21, cn, ns, udnei, W1, W2, a2, b2, n){
  #a2 = 1
  #b2 = 0.1
  shape <- n/2 + a2
  term <- as.vector(W2)-Z21%*%eta_21
  quad <- quad_new(Rho, n, cn, ns, udnei, W1, term)
  
  rate1 <- as.numeric(1/2*quad[1] + b2)
  rate2 <- as.numeric(1/2*quad[2] + b2)
  return(c(rgamma(1, shape = shape, rate = rate1),rgamma(1, shape = shape, rate = rate2)))
}

Z_matrix <- function(W1, ni_wo, udnei_wo, cni_wo){
  Z = rep(0, n)
  for(j in 1:ni_wo[1]){
    Z[1] = Z[1] + W1[udnei_wo[j]]
  }
  
  for(i in 2:n){
    for(j in 1:ni_wo[i]){
      Z[i] = Z[i] + W1[udnei_wo[(cni_wo[i-1] + j)]]
    }
  }
  Z21 = as.matrix(cbind(W1, Z))
  return(Z21)
}

eta_new <- function(Z21, W2, invD, sigmasq_eta1, sigmasq_eta2, tausq2){
  
  M <- solve(tausq2*t(Z21)%*%invD[[2]]%*%Z21 + diag(c(1/sigmasq_eta1, 1/sigmasq_eta2)))
  m <- tausq2 * t(Z21) %*% (invD[[2]] %*% W2)
  mu <- M %*% m
  return(rmvnorm(1, mean = mu, sigma = M))
}

target <- function(theta, W1, W2, eta_21, Z21, Minc, n, cn, ns, udnei, tausq1, tausq2, 
                   a1, b1, a2, b2){
  Rho <- exp(theta[1:2]) / (1 + exp(theta[1:2]))
  temp <- W2 - Z21 %*% eta_21
  quadf <- quad_new(Rho, n, cn, ns, udnei, W1, temp)
  PW <- 0.5*(quadf[4] + quadf[3]) - 0.5 * (tausq2 * quadf[2] + tausq1 * quadf[1])
  Plogitrho <- log(Rho[1]) + log(1-Rho[1]) + log(Rho[2]) + log(1-Rho[2])
  
  out <- PW + Plogitrho
  return(out)
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

gibbs_simpler <- function(Y1, Y2, X1, X2, n, q, cn, ns, udnei, sigmasq_beta, 
                          sigmasq_eta1, sigmasq_eta2, Minc, a1, b1, a2, b2,
                          ni_wo, udnei_wo, cni_wo,
                          time, burnin, thin, myscale){
  
  #Y <- c(Y1, Y2)
  #X <- as.matrix(bdiag(X1, X2))
  beta_old <- rep(0, 5)
  beta1_old <- beta_old[1:ncol(X1)]
  beta2_old <- beta_old[(ncol(X1)+1):(ncol(X1)+ncol(X2))]
  W1_old <- rep(0,n)
  W2_old <- rep(0,n)
  sigmasq1_old <- 1
  sigmasq2_old <- 1
  tausq1_old <- 1
  tausq2_old <- 1
  rho_old <- runif(2, 0, 1)
  eta_old <- c(1, 1)
  eta0_21_old <- eta_old[1]
  eta1_21_old <- eta_old[2]
  #time = 100
  #burnin = 10
  #thin=1
  t <- as.integer((time-burnin)/thin)
  beta1_t <- matrix(0, t, ncol(X1))
  beta2_t <- matrix(0, t, ncol(X2))
  W1_t <- matrix(0, t, n)
  W2_t <- matrix(0, t, n)
  sigmasq1_t <- rep(0, t)
  sigmasq2_t <- rep(0, t)
  tausq1_t <- rep(0, t)
  tausq2_t <- rep(0, t)
  rho1_t <- rep(0, t)
  rho2_t <- rep(0, t)
  eta0_21_t <- rep(0, t)
  eta1_21_t <- rep(0, t)
  
  accept_num <- 0
  myscale <- chol(myscale)
  
  invD <- Dinv_new(rho_old, n, cn, ns, udnei)
  Z21 <- as.matrix(Z_matrix(as.vector(W1_old), ni_wo, udnei_wo, cni_wo))
  
  for(i in 1:time){
    eta_old <- as.vector(eta_new(Z21, as.vector(W2_old), invD, sigmasq_eta1, sigmasq_eta2, tausq2_old))
    eta0_21_old <- eta_old[1]
    eta1_21_old <- eta_old[2]
    
    W1_old <- W1_new(Y1, X1, eta0_21_old, eta1_21_old, beta1_old, as.vector(W2_old), sigmasq1_old, tausq1_old, tausq2_old, invD)
    W2_old <- W2_new(Y2, X2, eta0_21_old, eta1_21_old, beta2_old, as.vector(W1_old), sigmasq2_old, tausq2_old, invD)
    
    sigmasq_old <- sigmasq_new(Y1, Y2, X1, X2, as.vector(W1_old), as.vector(W2_old), beta1_old, beta2_old, a1, b1, n)
    sigmasq1_old <- sigmasq_old[1]
    sigmasq2_old <- sigmasq_old[2]
    tausq_old <- tausq_new(rho_old, eta_old, Z21, cn, ns, udnei, as.vector(W1_old), as.vector(W2_old), a2, b2, n)
    tausq1_old <- tausq_old[1]
    tausq2_old <- tausq_old[2]
    #W_old <- W_new(Y, X, beta_old, sigmasq1_old, sigmasq2_old, tausq1_old, tausq2_old, invD, A21, Minc)
    #W1_old <- W_old[1:n]
    #W2_old <- W_old[-c(1:n)]
    
    beta_old <- beta_new(Y1, Y2, X1, X2, as.vector(W1_old), as.vector(W2_old), sigmasq1_old, sigmasq2_old, sigmasq_beta)
    beta1_old <- beta_old[1:ncol(X1)]
    beta2_old <- beta_old[(ncol(X1)+1):(ncol(X1)+ncol(X2))]
    
    current_x <- c(log(rho_old/(1-rho_old)))
    proposed_x = current_x + myscale%*%mvrnorm(1,rep(0,2),diag(2))
    
    A = min(1, exp(target(proposed_x, W1=as.vector(W1_old), W2=as.vector(W2_old), 
                          eta_21 = eta_old, Z21 = Z21, Minc = Minc, n=n, cn=cn, 
                          ns=ns, udnei=udnei, tausq1 = tausq1_old, tausq2 = tausq2_old,   
                          a1, b1, a2, b2) -
                     target(current_x, W1=as.vector(W1_old), W2=as.vector(W2_old),
                            eta_21 = eta_old, Z21 = Z21, Minc = Minc, n=n, cn=cn, 
                            ns=ns, udnei=udnei, tausq1 = tausq1_old, tausq2 = tausq2_old,
                            a1, b1, a2, b2)))
    if(runif(1)<A){
      out.batch = proposed_x       # accept move with probabily min(1,A)
      accept_num = accept_num + 1
    }else{
      out.batch = current_x        # otherwise "reject" move, and stay where we are
    }
    
    rho_old <- exp(out.batch[1:2])/(1+exp(out.batch[1:2]))
    
    accept_rate <- accept_num/i
    
    invD <- Dinv_new(rho_old, n, cn, ns, udnei)
    Z21 <- as.matrix(Z_matrix(as.vector(W1_old), ni_wo, udnei_wo, cni_wo))
    
    if (i>burnin && (i-burnin) %% thin== 0){
      
      tem <- (i - burnin)/thin
      beta1_t[tem,] <- beta1_old
      beta2_t[tem,] <- beta2_old
      W1_t[tem, ] <- W1_old
      W2_t[tem, ] <- W2_old
      sigmasq1_t[tem] <- sigmasq1_old
      sigmasq2_t[tem] <- sigmasq2_old
      tausq1_t[tem] <- tausq1_old
      tausq2_t[tem] <- tausq2_old
      rho1_t[tem] <- rho_old[1]
      rho2_t[tem] <- rho_old[2]
      eta0_21_t[tem] <- eta0_21_old
      eta1_21_t[tem] <- eta1_21_old
    }
  }
  print(accept_rate)
  return(list(par <- list(beta1 = beta1_t[,1], beta2 = beta1_t[,2], beta3 = beta2_t[,1], beta4 = beta2_t[,2], beta5 = beta2_t[,3],  
                          sigmasq1 = sigmasq1_t, sigmasq2 = sigmasq2_t, tausq1 = tausq1_t, tausq2 = tausq2_t, rho1 = rho1_t, rho2 = rho2_t,
                          eta0_21 = eta0_21_t, eta1_21 = eta1_21_t), W <- list(W1 = W1_t, W2 = W2_t)))
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
 
  
  set.seed(seed)
  tod1=Sys.time()  
  result <- gibbs_simpler(Y1=Y1, Y2=Y2, X1=X1, X2=X2, n=n, q=q, cn=cn, ns=ns, udnei=udnei, sigmasq_beta=1000, Minc = Minc, 
                          a1 = 2, b1 = 0.4, a2 = 2, b2 = 8, sigmasq_eta1 = 100, sigmasq_eta2 = 100,
                          ni_wo = ni_wo, udnei_wo = udnei_wo, cni_wo = cni_wo, 
                          myscale = diag(rep(6e-1),2), time = 30000, burnin = 15000, thin = 1)
  tod2=Sys.time() 
  
  
  sample.mcmc <- data.frame(matrix(unlist(result[[1]]), ncol = 13, byrow=F))
  colnames(sample.mcmc) <- c("beta1", "beta2","beta3","beta4","beta5","sigmasq1","sigmasq2",
                             "tausq1", "tausq2", "rho1", "rho2", "eta0_21", "eta1_21")
  
  estimate <- t(apply(sample.mcmc, 2, mysummary))
  
 
  
  
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
  
  w1hat[seed,] = apply(result[[2]]$W1, 2, mean)
  w1lower[seed,] = apply(result[[2]]$W1, 2, quantile, 0.025)
  w1upper[seed,] = apply(result[[2]]$W1, 2, quantile, 0.975)
  w2hat[seed,] = apply(result[[2]]$W2, 2, mean)
  w2lower[seed,] = apply(result[[2]]$W2, 2, quantile, 0.025)
  w2upper[seed,] = apply(result[[2]]$W2, 2, quantile, 0.975)
  

  
  PL_single <- matrix(0, nrow = 2*n, ncol = 15000)
  for(i in 1:15000){
    theta <- sample.mcmc[i,]
    
    beta1 <- theta[1:2]
    beta2 <- theta[3:5]
    sigmasq1 <- theta[6]
    sigmasq2 <- theta[7]
    
    for(j in 1:n){
      v1 <- sigmasq1
      PL_single[j,i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v1) - 1/(2*v1)*(Y1[j]-X1[j,]%*%t(as.vector(beta1)) - result[[2]]$W1[i,j])^2)
      v2 <- sigmasq2
      PL_single[(j+n),i] <- as.numeric(-1/2*log(2*pi) - 1/2*log(v2) - 1/(2*v2)*(Y2[j]-X2[j,]%*%t(as.vector(beta2)) - result[[2]]$W2[i,j])^2)
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
    Y_rep1[,i] <- X1 %*% t(as.vector(beta1)) +  result[[2]]$W1[i,] + sqrt(as.numeric(sigmasq1)) * z1
    Y_rep2[,i] <- X2 %*% t(as.vector(beta2)) +  result[[2]]$W2[i,] + sqrt(as.numeric(sigmasq2)) * z2
    
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
    W1vec <- result[[2]]$W1[i,]
    W2vec <- result[[2]]$W2[i,]
    beta1 <- as.vector(theta[1:2])
    beta2 <- as.vector(theta[3:5])
    sigmasq1 <- as.numeric(theta[6])
    sigmasq2 <- as.numeric(theta[7])
    
    KL_dagar[i] = kl.norm(mu1 = X%*%beta+W, S1 = sigma^2, mu2 = X%*%unlist(c(beta1, beta2)) + c(W1vec,W2vec), 
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
  
  saveRDS(tauwhat, paste("tauwhat_", modelgraph, num, ".rds", sep=""))
  saveRDS(tauwlower, paste("tauwlower_", modelgraph, num, ".rds", sep=""))
  saveRDS(tauwupper, paste("tauwupper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(rhohat, paste("rhohat_", modelgraph, num, ".rds", sep=""))
  saveRDS(rholower, paste("rholower_", modelgraph, num, ".rds", sep=""))
  saveRDS(rhoupper, paste("rhoupper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(etahat, paste("etahat_", modelgraph, num, ".rds", sep=""))
  saveRDS(etalower, paste("etalower_", modelgraph, num, ".rds", sep=""))
  saveRDS(etaupper, paste("etaupper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(w1hat, paste("w1hat_", modelgraph, num, ".rds", sep=""))
  saveRDS(w1lower, paste("w1lower_", modelgraph, num, ".rds", sep=""))
  saveRDS(w1upper, paste("w1upper_", modelgraph, num, ".rds", sep=""))
  saveRDS(w2hat, paste("w2hat_", modelgraph, num, ".rds", sep=""))
  saveRDS(w2lower, paste("w2lower_", modelgraph, num, ".rds", sep=""))
  saveRDS(w2upper, paste("w2upper_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(W1_path, paste("W1_", modelgraph, num, ".rds", sep=""))
  saveRDS(W2_path, paste("W2_", modelgraph, num, ".rds", sep=""))
  
  saveRDS(WAIC, paste("WAIC_", modelgraph, num, ".rds", sep=""))
  saveRDS(KL, paste("KL_", modelgraph, num, ".rds", sep=""))
  saveRDS(D_score, paste("D2_", modelgraph, num, ".rds", sep=""))
}

