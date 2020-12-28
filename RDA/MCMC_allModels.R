library(maps)
#Import California map
ca.county = map("county","california", fill=TRUE, plot=FALSE)

library(readr)
library(spdep)
library(maptools)
library(classInt)
library(RColorBrewer)
library(tidyr)
library(MASS)
library(Matrix)
library(mapproj)
library(rjags)
library(R2jags)
library(lattice)
library(mvtnorm)
library(matrixStats)
library(fields)
library(boot)
library(blockmatrix)
library(ggmcmc)
library(mcmc)
library(magic)
library(msos)
library(AICcmodavg)
library(coda)
library(invgamma)
library(mcmcse)
library(LaplacesDemon)
library(gtools)
library(ggmap)
setwd("/Users/Leiwen/Dropbox/Github/Multivariate_DAGAR/RDA")

#Import covariates
covariates <- read.csv("covariates.csv")
race <- read.csv("race.csv")
sex <- read.csv("sex.csv")
insurance <- read.csv("insurance.csv")
smoking <- read.csv("smoking.csv")
smoking$smoking <- as.numeric(substr(smoking$Cigarette.Smoking.Rate., 1,4))

#Import age-adjusted incidence rates for 4 cancers in California
rate_5y <- read.csv("age_adjusted.csv")
rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]

rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(extract_numeric(rate_lung$State_county)),]

rate_esophagus = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Esophagus",]
rate_esophagus = rate_esophagus[order(extract_numeric(rate_esophagus$State_county)),]

rate_larynx = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Larynx",]
rate_larynx = rate_larynx[order(extract_numeric(rate_larynx$State_county)),]

rate_colrect = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Colon and Rectum",]
rate_colrect = rate_colrect[order(extract_numeric(rate_colrect$State_county)),]

#County information
county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.poly$rate_lung = rate_lung$Age_Adjusted_Rate
ca.poly$rate_esophagus = rate_esophagus$Age_Adjusted_Rate
ca.poly$rate_larynx = rate_larynx$Age_Adjusted_Rate
ca.poly$rate_colrect = rate_colrect$Age_Adjusted_Rate
ca.poly$smoking = smoking$smoking

ca.coords = coordinates(ca.poly)


##################### Plot age-adjusted incidence rates in areal map ##################

brks_fit_lung = c(22, 41, 45, 51, 80)
brks_fit_esophagus = c(0, 3.5, 3.9, 4.5, 12)
brks_fit_larynx = c(0, 1.8, 2.1, 2.6, 5)
brks_fit_colrect = c(24, 34, 36, 38, 50)

color.pallete = rev(brewer.pal(4,"RdBu"))
class.rate_lung = classIntervals(var=ca.poly$rate_lung, n=4, style="fixed", 
                                 fixedBreaks=brks_fit_lung, dataPrecision=4)
class.rate_esophagus = classIntervals(var=ca.poly$rate_esophagus, n=4, style="fixed", 
                                      fixedBreaks=brks_fit_esophagus, dataPrecision=4)
class.rate_larynx = classIntervals(var=ca.poly$rate_larynx, n=4, style="fixed", 
                                   fixedBreaks=brks_fit_larynx, dataPrecision=4)
class.rate_colrect = classIntervals(var=ca.poly$rate_colrect, n=4, style="fixed", 
                                    fixedBreaks=brks_fit_colrect, dataPrecision=4)
color.code.rate_lung = findColours(class.rate_lung, color.pallete)
color.code.rate_esophagus = findColours(class.rate_esophagus, color.pallete)
color.code.rate_larynx = findColours(class.rate_larynx, color.pallete)
color.code.rate_colrect = findColours(class.rate_colrect, color.pallete)

pdf("incidence_rate_5y_new.pdf", height = 10, width = 10)
par(mfrow=c(2,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(ca.poly, col = color.code.rate_lung)
leg.txt1 = c("22-41", "41-45", "45-51","51-80")
legend("bottomleft", title="Lung cancer", legend=leg.txt1, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.rate_esophagus)
leg.txt2 = c("0-3.5", "3.5-3.9", "3.9-4.5", "4.5-12")
legend("bottomleft", title="Esophagus cancer",legend=leg.txt2, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.rate_larynx)
leg.txt3 = c("0-1.8", "1.8-2.1", "2.1-2.6", "2.6-5")
legend("bottomleft", title="Larynx cancer",legend=leg.txt3, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)


plot(ca.poly, col = color.code.rate_colrect)
leg.txt4 = c("24-34", "34-36", "36-38", "38-50")
legend("bottomleft",title="Colorectum cancer", legend=leg.txt4, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

dev.off()

################## Data construction and adjacency matrix ###########

## Data

# Covariates in percentage
county_attribute = covariates[substr(covariates$State_county,1,2) == "CA",]
county_attribute$state = extract_numeric(county_attribute$State_county)
county_attribute1 = data.frame(county_attribute[order(county_attribute$state),])
county_attribute1$V_Persons_age_18_ACS_2012_2016 = as.numeric(county_attribute1$V_Persons_age_18_ACS_2012_2016)/100
county_attribute1$V_Persons_age_65_ACS_2012_2016 = as.numeric(county_attribute1$V_Persons_age_65_ACS_2012_2016)/100
county_attribute1$VHighschooleducationACS2012201 = as.numeric(county_attribute1$VHighschooleducationACS2012201)/100
county_attribute1$VFamiliesbelowpovertyACS201220 = as.numeric(county_attribute1$VFamiliesbelowpovertyACS201220)/100
county_attribute1$V_Unemployed_ACS_2012_2016 = as.numeric(county_attribute1$V_Unemployed_ACS_2012_2016)/100
race1 = race[substr(race$State_county,1,2) == "CA"&race$Race_recode_White_Black_Other=="Black",]
sex1 = sex[substr(sex$State_county,1,2) == "CA"&sex$Sex=="Male",]
insurance1 = insurance[substr(insurance$State_county,1,2) == "CA"&insurance$Insurance_Recode_2007=="Uninsured",]

# Create covariates matrix for each cancer
# For Case 1, all covariates are included; But smoking is excluded in covariates for Case 2.
rate_lung1 = cbind(rate_lung, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_lung1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_esophagus1 = cbind(rate_esophagus, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_esophagus1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_larynx1 = cbind(rate_larynx, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_larynx1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_colrect1 = cbind(rate_colrect, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_colrect1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")

# Areal map for three important covariates: smoking, black and uninsured
ca.poly$black = race1$Row_Percent
ca.poly$smoke = smoking$smoking
ca.poly$uninsure = insurance1$Row_Percent

brks_fit_black = quantile(ca.poly$black)
class.black = classIntervals(var=ca.poly$black, n=4, style="fixed", 
                             fixedBreaks=brks_fit_black, dataPrecision=4)
color.code.black = findColours(class.black, color.pallete)
brks_fit_smoke = quantile(ca.poly$smoke)
class.smoke = classIntervals(var=ca.poly$smoke, n=4, style="fixed", 
                             fixedBreaks=brks_fit_smoke, dataPrecision=4)
color.code.smoke = findColours(class.smoke, color.pallete)
rks_fit_uninsure = quantile(ca.poly$uninsure)
class.uninsure = classIntervals(var=ca.poly$uninsure, n=4, style="fixed",
                             fixedBreaks=brks_fit_uninsure, dataPrecision=4)
color.code.uninsure = findColours(class.uninsure, color.pallete)

pdf("covariates.pdf", height = 6, width = 10)
par(mfrow=c(1,3), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(ca.poly, col = color.code.smoke)
leg.txt = c("6.70 - 11.50", "11.50 - 13.85", "13.85 - 16.28", "16.28 - 25.50")
legend("bottomleft", title="Smoke (%)", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE,
       fill = color.pallete)

plot(ca.poly, col = color.code.black)
leg.txt = c("0.90 - 1.90", "1.90 - 2.80", "2.80 - 5.15", "5.15 - 16.90")
legend("bottomleft", title="Black (%)", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE,
       fill = color.pallete)

plot(ca.poly, col = color.code.uninsure)
leg.txt = c("0.0 - 0.6", "0.6 - 0.9", "0.9 - 1.4", "1.4 - 3.8")
legend("bottomleft", title="Uninsured (%)", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE,
       fill = color.pallete)
dev.off()


## Create adjacency matrix and neighbor info
ca.neighbors = poly2nb(ca.poly)
n=length(ca.neighbors)

Adj=sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Adj)=county.ID
ca.coord = coordinates(ca.poly)
ca.latrange=round(quantile(ca.coord[,2],c(0.25,0.75)))
ca.albersproj=mapproject(ca.coord[,1],ca.coord[,2],projection = "albers",param=ca.latrange)

perm=order(ca.albersproj$x-ca.albersproj$y)
colnames(Adj)[perm]

Adj_new=Adj[perm,perm]

n=nrow(Adj_new)
ni=rowSums(Adj_new)
maxn=max(ni)
neimat=matrix(0,n,maxn)
neighbors=lapply(1:n,function(x) which(Adj_new[x,]==1))
#N(i): 2:n
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
#n<i: 2:n
dni=sapply(dneighbors,length)
original_perm = 1:58
index2=c(1,which(dni==0)+1)

final_perm=c(original_perm[perm][index2],
             original_perm[perm][-index2])
final_perm[order(final_perm)]

Minc = Adj[final_perm,final_perm]
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

ni_wo = sapply(neighbors,length)
cni_wo = cumsum(ni_wo)
udnei_wo = unlist(neighbors)
cn = c(0, cni)
ns = dni

region = seq(1:n)
cn = c(0, cni)
ns = dni
index = list()
for(i in 1:(n-2)){
  index[[i]] = region[-(udnei[(cn[i+1] + 1):(cn[i+1] + ns[i+1])])]
}
index1 = unlist(index)

# Implement MDAGAR for all permutations

q=4
rate_list = list(rate_lung1, rate_esophagus1, rate_larynx1, rate_colrect1)

orderindex = c(1,2,3,4)
models = permutations(n=4,r=4,v=orderindex,repeats.allowed=F)

mcmc_list = list()
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
  
  Y = c(Y1,Y2,Y3,Y4)
  X = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))
  
  nc = ncol(X)
  
  tod1=Sys.time()
  
  sink("MDAGAR.txt")
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
      beta[1:nc] ~ dmnorm(rep(0,nc), (0.001*I1))
      }
      ", fill = TRUE)
  sink()
  
  model.data <- list(k = n, index1 = index1, I = diag(n), nc = nc, I1 = diag(nc), Minc = Minc, ns = dni, cn = c(0, cni), udnei = udnei, X = X, Y = Y)
  model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, rho3 = 0.1, tau1 = 1, tau2 = 1, tau3 = 1, eta021 = 1, 
                               eta121 = 1, eta031 = 1, eta131 = 1, eta032 = 1, eta132 = 1, eta041 = 1, 
                               eta141 = 1, eta042 = 1, eta142 = 1, eta043 = 1, eta143 = 1, taue1 = 1, taue2 = 1,
                               taue3 = 1, taue4 = 1, beta = rep(0, nc), W1 = rep(0.1, n), 
                               W2 = rep(0.1, n), W3 = rep(0.1, n), W4 = rep(0.1, n))),2)
  model.param <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", "eta021", "eta121",
                   "eta031", "eta131", "eta032", "eta132","eta041", "eta141",
                   "eta042", "eta142", "eta043", "eta143", "vare1", "vare2", "vare3", "vare4", "W", "loglik")
  model.param2 <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", "eta021", "eta121",
                    "eta031", "eta131", "eta032", "eta132","eta041", "eta141",
                    "eta042", "eta142", "eta043", "eta143", "vare1", "vare2", "vare3", "vare4")

  set.seed(125)
  #Results for post-analysis
  result1 <- jags(model.data, model.inits, model.param, "MDAGAR.txt",
                  n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  #Results for model selection
  result2 <- jags(model.data, model.inits, model.param2, "DAGAR.txt",
                  n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  tod2=Sys.time()
  
  mcmc_list[[i]] = result1
  mcmc_list2[[i]] = result2
}

# Save MCMC outputs of all permutations for model selection
saveRDS(mcmc_list, "mcmc_list_new.rds")
saveRDS(mcmc_list2, "mcmc_list_new2.rds")

