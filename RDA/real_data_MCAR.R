library(maps)

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
library(Hmisc)
library("PerformanceAnalytics")
setwd("Multivariate_DAGAR/RDA")

#Import covariates
covariates <- read.csv("data/covariates.csv")
race <- read.csv("data/race.csv")
sex <- read.csv("data/sex.csv")
insurance <- read.csv("data/insurance.csv")
smoking <- read.csv("data/smoking.csv")
smoking$smoking <- as.numeric(substr(smoking$Cigarette.Smoking.Rate., 1,4))

#Import age-adjusted incidence rates for 4 cancers in California
rate_5y <- read.csv("data/age_adjusted.csv")

rate_CA = rate_5y[substr(rate_5y$State_county,1,2) == "CA",]

rate_lung = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Lung and Bronchus",]
rate_lung = rate_lung[order(readr::parse_number(as.character(rate_lung$State_county))),]

rate_esophagus = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Esophagus",]
rate_esophagus = rate_esophagus[order(readr::parse_number(as.character(rate_esophagus$State_county))),]

rate_larynx = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Larynx",]
rate_larynx = rate_larynx[order(readr::parse_number(as.character(rate_larynx$State_county))),]

rate_colrect = rate_CA[rate_CA$Site_recode_ICD_O_3_WHO_2008=="Colon and Rectum",]
rate_colrect = rate_colrect[order(readr::parse_number(as.character(rate_colrect$State_county))),]

cancer_rates = cbind(rate_lung$Age_Adjusted_Rate, rate_esophagus$Age_Adjusted_Rate, rate_larynx$Age_Adjusted_Rate, rate_colrect$Age_Adjusted_Rate)
colnames(cancer_rates) = c("Lung and Bronchus", "Esophageal", "Larynx", "Colon and Rectum")
cor_matrix <- rcorr(cancer_rates)
cor_matrix$r
cor_matrix$P
chart.Correlation(cancer_rates, histogram=TRUE, pch=19)

county.ID <- sapply(strsplit(ca.county$names, ","), function(x) x[2])
ca.poly = map2SpatialPolygons(ca.county, IDs=county.ID)
ca.poly$rate_lung = rate_lung$Age_Adjusted_Rate
ca.poly$rate_esophagus = rate_esophagus$Age_Adjusted_Rate
ca.poly$rate_larynx = rate_larynx$Age_Adjusted_Rate
ca.poly$rate_colrect = rate_colrect$Age_Adjusted_Rate
ca.poly$smoking = smoking$smoking

ca.coords = coordinates(ca.poly)


##################### Incidence rate in areal map ##################
#brks_fit_lung = c(22, quantile(ca.poly$rate_lung)[2:4], 73)
#brks_fit_esophagus = c(0, quantile(ca.poly$rate_esophagus)[2:4], 12)
#brks_fit_larynx = c(0, quantile(ca.poly$rate_larynx)[2:4], 5)
#brks_fit_colrect = c(24, quantile(ca.poly$rate_colrect)[2:4], 50)
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





################## Data construction and adjacency matrix ###########

## Data

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

rate_lung1 = cbind(rate_lung, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_lung1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_esophagus1 = cbind(rate_esophagus, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_esophagus1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_larynx1 = cbind(rate_larynx, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_larynx1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_colrect1 = cbind(rate_colrect, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_colrect1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")

#fit = lm(rate~smoking+young+old+highschool+poverty+unemployed+black+male+uninsured, data = rate_colrect1)
#fit1 = lm(rate~old, data = rate_colrect1)


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

brks_fit_uninsure = quantile(ca.poly$uninsure)
class.uninsure = classIntervals(var=ca.poly$uninsure, n=4, style="fixed", 
                                fixedBreaks=brks_fit_uninsure, dataPrecision=4)
color.code.uninsure = findColours(class.uninsure, color.pallete)


## Adjacency matrix
ca.neighbors = poly2nb(ca.poly)
n=length(ca.neighbors)

Adj=sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Adj)=county.ID
ca.coord = coordinates(ca.poly)
ca.latrange=round(quantile(ca.coord[,2],c(0.25,0.75)))
ca.albersproj=mapproject(ca.coord[,1],ca.coord[,2],projection = "albers",param=ca.latrange)
projmat=cbind(ca.albersproj$x,ca.albersproj$y)
dmat=as.matrix(dist(projmat))



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

q=4
rate_list = list(rate_lung1, rate_esophagus1, rate_larynx1, rate_colrect1)

#for(i in 1:24){
  print(i)
  rate_list1 = rate_list
  
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
  
  tod1=Sys.time()
  
  sink("JCAR.txt")
  cat("
      model
      {
      Q1[1:k, 1:k] <- D - rho1*Minc
      Q2[1:k, 1:k] <- D - rho2*Minc
      Q3[1:k, 1:k] <- D - rho3*Minc
      Q4[1:k, 1:k] <- D - rho4*Minc
      
      f1[1:k] ~ dmnorm(rep(0, k), Q1)
      f2[1:k] ~ dmnorm(rep(0, k), Q2)
      f3[1:k] ~ dmnorm(rep(0, k), Q3)
      f4[1:k] ~ dmnorm(rep(0, k), Q4)
      
      W[1:k] <- A[1,1]*f1 
      W[(k+1):(2*k)] <- A[2,1] * f1 + A[2,2] * f2
      W[(2*k+1):(3*k)] <- A[3,1] * f1 + A[3,2] * f2 + A[3,3] * f3
      W[(3*k+1):(4*k)] <- A[4,1] * f1 + A[4,2] * f2 + A[4,3] * f3 + A[4,4] * f4
      
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

      for (i in (2*k+1):(3*k))
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue3)
      }

      for (i in (3*k+1):(4*k))
      {
      mu[i] <- X[i,] %*% beta + W[i]
      Y[i] ~ dnorm(mu[i], taue4)
      }
      
      for(r in 1:4){
      A[r, r] ~ dlnorm(0, 1)
      }
      
      for(r in 2:4){
         for(l in 1:(r-1)){
             A[r, l] ~ dnorm(0, 0.01)
         }
      }

      for(r in 1:3){
         for(l in (r+1):4){
           A[r, l] <- 0
         }
      }
      
      rho1 ~ dunif(0, 0.999)
      rho2 ~ dunif(0, 0.999)
      rho3 ~ dunif(0, 0.999)
      rho4 ~ dunif(0, 0.999)
      
      taue1 ~ dgamma(2, 1)
      taue2 ~ dgamma(2, 1)
      taue3 ~ dgamma(2, 1)
      taue4 ~ dgamma(2, 1)
      vare1 <- 1/taue1
      vare2 <- 1/taue2
      vare3 <- 1/taue3
      vare4 <- 1/taue4
      
      beta[1:40] ~ dmnorm(rep(0,40), (0.001*I1))
      }
      ", fill = TRUE)
  sink()
  
  model.data <- list(k = n, I1 = diag(40), Minc = Minc, D = diag(rowSums(Minc)), X = X, Y = Y)
  model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, rho3 = 0.1, rho4 = 0.1, taue1 = 1, taue2 = 1,
                               taue3 = 1, taue4 = 1, beta = rep(0, 40))),2)
  # model.param1 <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", "eta021", "eta121",
  #                 "eta031", "eta131", "eta032", "eta132","eta041", "eta141",
  #                "eta042", "eta142", "eta043", "eta143", "vare1", "vare2", "vare3", "vare4", "W", "loglik")
  model.param2 <- c("beta", "A", "rho1", "rho2", "rho3", "rho4",
                    "vare1", "vare2", "vare3", "vare4", "W")
  set.seed(125)
  #result1 <- jags(model.data, model.inits, model.param1, "DAGAR.txt",
  #               n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  result2 <- jags(model.data, model.inits, model.param2, "JCAR.txt",
                  n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
  tod2=Sys.time()
  
#}

#saveRDS(result2, "mcmc_result_jcar.rds")


mysummary = function(invector) {
  c(mean(invector), median(invector), sd(invector), quantile(invector, .025), quantile(invector,.975))
}

mcmc_select <- result2
jags_result1 = mcmc_select$BUGSoutput$sims.array[,1,c(249:297,1:4,6:8,11,12,16)]
jags_result2 = mcmc_select$BUGSoutput$sims.array[,2,c(249:297,1:4,6:8,11,12,16)]
jags_result <- rbind(jags_result1, jags_result2)

W_mcmc1 =  mcmc_select$BUGSoutput$sims.array[,1,17:248]
W_mcmc2 =  mcmc_select$BUGSoutput$sims.array[,2,17:248]
W_mcmc = rbind(W_mcmc1, W_mcmc2)
colnames(jags_result) = c("beta1", "beta2","beta3","beta4","beta5", "beta6", "beta7","beta8","beta9","beta10","beta11","beta12",
                          "beta13", "beta14","beta15","beta16","beta17", "beta18", "beta19","beta20","beta21","beta22","beta23","beta24",
                          "beta25", "beta26","beta27","beta28","beta29", "beta30", "beta31","beta32","beta33","beta34","beta35","beta36",
                          "beta37","beta38","beta39","beta40", "deviance", "rho1", "rho2", "rho3", "rho4",
                          "sigmasq1","sigmasq2","sigmasq3","sigmasq4", "a11", "a21", "a31", "a41", "a22", "a32", "a42",
                          "a33", "a43", "a44")
estimate1 <- round(t(apply(jags_result, 2, mysummary)),2)

table_est = paste(estimate1[,1], " (", estimate1[,4], ", ", estimate1[,5], ")", sep="")
nrow(estimate1)

table1 <- data.frame(matrix(table_est[1:40], ncol = 4, byrow=F))
table2 <- data.frame(matrix(table_est[42:49], ncol = 4, byrow=T))

table <- rbind(table1, table2)
colnames(table) <- c("Lung", "Esophagus", "Larynx", "Colorectum")
write.csv(table, "table_est_jcar.csv")





W_mean = apply(W_mcmc, 2, mean)
W_mean1 = W_mean[1:58]
ca.poly$W_mean1 = W_mean1[order(final_perm)]
#brks_fit1 = quantile(W_mean1)
brks_fit1 = c(-23, -10, 0, 5, 10, 20)
color.pallete = rev(brewer.pal(5,"RdBu"))
class.W1 = classIntervals(var=ca.poly$W_mean1, n=5, style="fixed", 
                          fixedBreaks=brks_fit1, dataPrecision=4)
color.code.W1 = findColours(class.W1, color.pallete)

W_mean2 = W_mean[59:116]
ca.poly$W_mean2 = W_mean2[order(final_perm)]
#brks_fit2 = quantile(W_mean2)
brks_fit2 = c(-2, -0.5, 0, 0.1, 0.5, 1.5)
class.W2 = classIntervals(var=ca.poly$W_mean2, n=5, style="fixed", 
                          fixedBreaks=brks_fit2, dataPrecision=4)
color.code.W2 = findColours(class.W2, color.pallete)

W_mean3 = W_mean[117:174]
ca.poly$W_mean3 = W_mean3[order(final_perm)]
#brks_fit3 = quantile(W_mean3)
brks_fit3 = c(-1, -0.5, 0, 0.1, 0.5, 1.1)
class.W3 = classIntervals(var=ca.poly$W_mean3, n=5, style="fixed", 
                          fixedBreaks=brks_fit3, dataPrecision=4)
color.code.W3 = findColours(class.W3, color.pallete)

W_mean4 = W_mean[175:232]
ca.poly$W_mean4 = W_mean4[order(final_perm)]
#brks_fit4 = quantile(W_mean4)
brks_fit4 = c(-10, -3, 0, 3, 6, 12)
class.W4 = classIntervals(var=ca.poly$W_mean4, n=5, style="fixed", 
                          fixedBreaks=brks_fit4, dataPrecision=4)
color.code.W4 = findColours(class.W4, color.pallete)

# Plots for random effects
pdf("random_covariates_jcar.pdf", height = 10, width = 10)
par(mfrow=c(2,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(ca.poly, col = color.code.W1)
leg.txt1 = c("-23 - -10", "-10 - 0", "0 - 5", "5 - 10", "10 - 20")
legend("bottomleft", title="Lung cancer",legend=leg.txt1, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W2)
leg.txt2 = c("-2 - -0.5", "-0.5 - 0", "0 - 0.1", "0.1 - 0.5", "0.5 - 1.5")
legend("bottomleft", title="Esophageal cancer", legend=leg.txt2, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W3)
leg.txt3 = c("-1 - -0.5", "-0.5 - -0", "0 - 0.1", "0.1 - 0.5", "0.5 - 1.1") 
legend("bottomleft", title="Larynx cancer", legend=leg.txt3, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W4)
leg.txt4 = c("-10 - -3", "-3 - 0", "0 - 3", "3 - 6", "6 - 12")
legend("bottomleft",title="Colorectum cancer", legend=leg.txt4, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)
dev.off()


######## Posterior prediction #############
rate_list1 = rate_list

Y1 = rate_list1[[1]]$rate[final_perm]
Y2 = rate_list1[[2]]$rate[final_perm]
Y3 = rate_list1[[3]]$rate[final_perm]
Y4 = rate_list1[[4]]$rate[final_perm]

X1 = as.matrix(cbind(1, rate_list1[[1]][,6:14]))[final_perm,]
X2 = as.matrix(cbind(1, rate_list1[[2]][,6:14]))[final_perm,]
X3 = as.matrix(cbind(1, rate_list1[[3]][,6:14]))[final_perm,]
X4 = as.matrix(cbind(1, rate_list1[[4]][,6:14]))[final_perm,]

Y = c(Y1,Y2,Y3,Y4)
X = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))

y.est = matrix(0,30000,4*n)
for(iter in 1:30000){
  print(iter)
  y.est[iter,] = as.vector(X %*% as.vector(jags_result[iter, 1:40]) + W_mcmc[iter,])
}
y.mean = apply(y.est,2,mean)




ca.poly$rate_lung_est = y.mean[1:58][order(final_perm)]
ca.poly$rate_esophagus_est = y.mean[59:116][order(final_perm)]
ca.poly$rate_larynx_est = y.mean[117:174][order(final_perm)]
ca.poly$rate_colrect_est = y.mean[175:232][order(final_perm)]


brks_fit_lung = c(22, 41, 45, 51, 80)
brks_fit_esophagus = c(0, 3.5, 3.9, 4.5, 12)
brks_fit_larynx = c(0, 1.8, 2.1, 2.6, 5)
brks_fit_colrect = c(24, 34, 36, 38, 50)


color.pallete = rev(brewer.pal(4,"RdBu"))
class.rate_lung = classIntervals(var=ca.poly$rate_lung_est, n=4, style="fixed", 
                                 fixedBreaks=brks_fit_lung, dataPrecision=4)
class.rate_esophagus = classIntervals(var=ca.poly$rate_esophagus_est, n=4, style="fixed", 
                                      fixedBreaks=brks_fit_esophagus, dataPrecision=4)
class.rate_larynx = classIntervals(var=ca.poly$rate_larynx_est, n=4, style="fixed", 
                                   fixedBreaks=brks_fit_larynx, dataPrecision=4)
class.rate_colrect = classIntervals(var=ca.poly$rate_colrect_est, n=4, style="fixed", 
                                    fixedBreaks=brks_fit_colrect, dataPrecision=4)
color.code.rate_lung = findColours(class.rate_lung, color.pallete)
color.code.rate_esophagus = findColours(class.rate_esophagus, color.pallete)
color.code.rate_larynx = findColours(class.rate_larynx, color.pallete)
color.code.rate_colrect = findColours(class.rate_colrect, color.pallete)


pdf("post_incidence_jcar.pdf", height = 10, width = 10)
par(mfrow=c(2,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)


plot(ca.poly, col = color.code.rate_lung)
leg.txt1 = c("22-41", "41-45", "45-51","51-80")
legend("bottomleft", title="Lung cancer", legend=leg.txt1, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.rate_esophagus)
leg.txt2 = c("0-3.5", "3.5-3.9", "3.9-4.5", "4.5-12")
legend("bottomleft", title="Esophageal cancer",legend=leg.txt2, cex=1.25, bty="n", horiz = FALSE, 
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
