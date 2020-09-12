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
# For case 1, smoking is excluded in covariates; All covariates are included in case 2.
rate_lung1 = cbind(rate_lung, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_lung1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_esophagus1 = cbind(rate_esophagus, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_esophagus1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_larynx1 = cbind(rate_larynx, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_larynx1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")
rate_colrect1 = cbind(rate_colrect, smoking$smoking, county_attribute1[,2:6], race1$Row_Percent, sex1$Row_Percent,insurance1$Row_Percent)
colnames(rate_colrect1) = c("county", "site", "rate", "count", "population", "smoking", "young","old", "highschool", "poverty", "unemployed", "black", "male", "uninsured")

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
dneighbors=sapply(2:n,function(i) intersect(neighbors[[i]],1:(i-1)))
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
  V_est = L_est%*%G_est%*%t(L_est)
  #V_est = diag(c(rep(x[["sigmasq1"]],n),rep(x[["sigmasq2"]],n),rep(x[["sigmasq3"]],n),rep(x[["sigmasq4"]],n))) + 
  # L_est%*%G_est%*%t(L_est)
  
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

#Read in MCMC output for models
mcmc_list1 = readRDS("mcmc_list_new.rds")


# Case 2: model 15  (Case1: model 16) 
mcmc_select <- mcmc_list1[[15]]
jags_result1 = mcmc_select$BUGSoutput$sims.array[,1,c(233:281,514:529)]
jags_result2 = mcmc_select$BUGSoutput$sims.array[,2,c(233:281,514:529)]
jags_result = rbind(jags_result1, jags_result2)

W_mcmc1 =  mcmc_select$BUGSoutput$sims.array[,1,1:232]
W_mcmc2 =  mcmc_select$BUGSoutput$sims.array[,2,1:232]
W_mcmc = rbind(W_mcmc1, W_mcmc2)
colnames(jags_result) = c("beta1", "beta2","beta3","beta4","beta5", "beta6", "beta7","beta8","beta9","beta10","beta11","beta12",
                          "beta13", "beta14","beta15","beta16","beta17", "beta18", "beta19","beta20","beta21","beta22","beta23","beta24",
                          "beta25", "beta26","beta27","beta28","beta29", "beta30", "beta31","beta32","beta33","beta34","beta35","beta36",
                          "beta37","beta38","beta39","beta40", "deviance", "eta0_21", "eta0_31", "eta0_32", "eta0_41",
                          "eta0_42", "eta0_43", "eta1_21", "eta1_31", "eta1_32", "eta1_41", "eta1_42", "eta1_43", "rho1", "rho2", "rho3", "rho4",
                          "tausq1", "tausq2", "tausq3", "tausq4","sigmasq1","sigmasq2","sigmasq3","sigmasq4")
estimate1 <- round(t(apply(jags_result, 2, mysummary)),2)

# Table for coefficients
table_est = paste(estimate1[,1], " (", estimate1[,4], ", ", estimate1[,5], ")", sep="")
nrow(estimate1)

table1 <- data.frame(matrix(table_est[1:40], ncol = 4, byrow=F))
table2 <- data.frame(matrix(table_est[54:65], ncol = 4, byrow=T))

table <- rbind(table1, table2)
colnames(table) <- c("Larynx", "Esophagus", "Lung", "Colorectum")
write.csv(table, "table_est_new.csv")

#Histogram for etas
pdf("eta_final_new.pdf", height = 10, width = 10)
par(mfrow=c(3,4))
hist(jags_result[,42],xlab = "eta_021, esophagus | larynx",ylab = "", main="")
hist(jags_result[,48],xlab = "eta_121, esophagus | larynx",ylab = "", main="")
hist(jags_result[,43],xlab = "eta_031, lung | larynx",ylab = "", main="")
hist(jags_result[,49],xlab = "eta_131, lung | larynx",ylab = "", main="")
hist(jags_result[,45],xlab = "eta_041, colorectum | larynx",ylab = "", main="")
hist(jags_result[,51],xlab = "eta_141, colorectum | larynx",ylab = "", main="")
hist(jags_result[,44],xlab = "eta_032, lung | esophagus",ylab = "", main="")
hist(jags_result[,50],xlab = "eta_132, lung | esophagus",ylab = "", main="")
hist(jags_result[,46],xlab = "eta_042, colorectum | esophagus",ylab = "", main="")
hist(jags_result[,52],xlab = "eta_142, colorectum | esophagus",ylab = "", main="")
hist(jags_result[,47],xlab = "eta_043, colorectum | lung",ylab = "", main="")
hist(jags_result[,53],xlab = "eta_143, colorectum | lung",ylab = "", main="")
dev.off()


################## Correlation in each county ###############
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

## Plots for correlation in each county
brks_fit1 = c(-0.6, -0.3, -0.1, 0, 0.1, 0.25, 0.75)
color.pallete = rev(brewer.pal(6,"RdBu"))
class.W1 = classIntervals(var=ca.poly$W_mean1, n=6, style="fixed", 
                          fixedBreaks=brks_fit1, dataPrecision=4)
color.code.W1 = findColours(class.W1, color.pallete)

brks_fit = c(-0.1, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1)
color.pallete = brewer.pal(7,"Blues")
class.rate12 = classIntervals(var=ca.poly$rate12, n=7, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate13 = classIntervals(var=ca.poly$rate13, n=7, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate14 = classIntervals(var=ca.poly$rate14, n=7, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate23 = classIntervals(var=ca.poly$rate23, n=7, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate24 = classIntervals(var=ca.poly$rate24, n=7, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
class.rate34 = classIntervals(var=ca.poly$rate34, n=7, style="fixed", 
                              fixedBreaks=brks_fit, dataPrecision=4)
color.code.rate12 = findColours(class.rate12, color.pallete)
color.code.rate13 = findColours(class.rate13, color.pallete)
color.code.rate14 = findColours(class.rate14, color.pallete)
color.code.rate23 = findColours(class.rate23, color.pallete)
color.code.rate24 = findColours(class.rate24, color.pallete)
color.code.rate34 = findColours(class.rate34, color.pallete)

pdf("correlation_pairwise_15.pdf", height = 10, width = 10)
par(mfrow=c(3,2), oma = c(4, 1, 1, 1))
leg.txt = c("-0.1 - 0.1", "0.1 - 0.3", "0.3 - 0.5", "0.5 - 0.7", "0.7 - 0.8", "0.8 - 0.9", "0.9 - 1.0") 

plot(ca.poly, col = color.code.rate12)
title(sub="(a) larynx cancer and esophageal cancer")
points(ca.coords[which(ca.poly$sig12==1),1], ca.coords[which(ca.poly$sig12==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate14)
title(sub="(b) larynx cancer and colorectum cancer" )
points(ca.coords[which(ca.poly$sig14==1),1], ca.coords[which(ca.poly$sig14==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate13)
title(sub="(c) larynx cancer and lung cancer" )
points(ca.coords[which(ca.poly$sig13==1),1], ca.coords[which(ca.poly$sig13==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate24)
title(sub="(d) esophageal cancer and colorectum cancer" )
points(ca.coords[which(ca.poly$sig24==1),1], ca.coords[which(ca.poly$sig24==1),2], pch=10, col="yellow", cex=0.2)

plot(ca.poly, col = color.code.rate23)
title(sub="(e) esophageal cancer and lung cancer" )
points(ca.coords[which(ca.poly$sig23==1),1], ca.coords[which(ca.poly$sig23==1),2], pch=10, col="yellow", cex=0.2)


plot(ca.poly, col = color.code.rate34)
title(sub="(f) colorectum cancer and lung cancer" )
points(ca.coords[which(ca.poly$sig34==1),1], ca.coords[which(ca.poly$sig34==1),2], pch=10, col="yellow", cex=0.2)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("center", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

dev.off()

############# Random effects ###############
W_mean = apply(W_mcmc, 2, mean)
W_mean1 = W_mean[1:58]
ca.poly$W_mean1 = W_mean1[order(final_perm)]

W_mean2 = W_mean[59:116]
ca.poly$W_mean2 = W_mean2[order(final_perm)]
brks_fit2 = c(-2.5, -0.5, 0, 0.35, 1, 2, 3)
class.W2 = classIntervals(var=ca.poly$W_mean2, n=6, style="fixed", 
                          fixedBreaks=brks_fit2, dataPrecision=4)
color.code.W2 = findColours(class.W2, color.pallete)

W_mean3 = W_mean[117:174]
ca.poly$W_mean3 = W_mean3[order(final_perm)]
brks_fit3 = c(-21, -2, 2, 5, 10, 20, 31)
class.W3 = classIntervals(var=ca.poly$W_mean3, n=6, style="fixed", 
                          fixedBreaks=brks_fit3, dataPrecision=4)
color.code.W3 = findColours(class.W3, color.pallete)

W_mean4 = W_mean[118:232]
ca.poly$W_mean4 = W_mean4[order(final_perm)]
brks_fit4 = c(-21, -1, 1, 3, 10, 20, 31)
class.W4 = classIntervals(var=ca.poly$W_mean4, n=6, style="fixed", 
                          fixedBreaks=brks_fit4, dataPrecision=4)
color.code.W4 = findColours(class.W4, color.pallete)
# Plots for random effects
pdf("random_covariates.pdf", height = 10, width = 10)
par(mfrow=c(2,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(ca.poly, col = color.code.W3)
leg.txt3 = c("-21 - -2", "-2 - 2", "2 - 5", "5 - 10", "10 - 20", "20 - 31")
legend("bottomleft", title="Lung cancer",legend=leg.txt3, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W2)
leg.txt = c("-2.5 - -0.5", "-0.5 - 0", "0 - 0.35", "0.35 - 1", "1 - 2", "2 - 3") 
legend("bottomleft", title="Esophageal cancer", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W1)
leg.txt = c("-0.6 - -0.3", "-0.3 - -0.1", "-0.1 - 0", "0 - 0.1", "0.1 - 0.25", "0.25 - 0.75") 
legend("bottomleft", title="Larynx cancer", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W4)
leg.txt4 = c("-21 - -1", "-1 - 1", "1 - 3", "3 - 10", "10 - 20", "20 - 31")
legend("bottomleft",title="Colorectum cancer", legend=leg.txt4, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

dev.off()

######## Posterior prediction #############
rate_list1 = rate_list[models[15,]]

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

ca.poly$rate_larynx_est = y.mean[1:58][order(final_perm)]
ca.poly$rate_esophagus_est = y.mean[59:116][order(final_perm)]
ca.poly$rate_lung_est = y.mean[117:174][order(final_perm)]
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


pdf("post_incidence_new.pdf", height = 10, width = 10)
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

###### Same Model without covariates #########

i=15
rate_list1 = rate_list[models[i,]]

Y1 = rate_list1[[1]]$rate[final_perm]
Y2 = rate_list1[[2]]$rate[final_perm]
Y3 = rate_list1[[3]]$rate[final_perm]
Y4 = rate_list1[[4]]$rate[final_perm]


Y = c(Y1,Y2,Y3,Y4)
#X = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))

tod1=Sys.time()

sink("DAGAR_nocov.txt")
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
    mu[i] <- beta[1] + W[i]
    Y[i] ~ dnorm(mu[i], taue1)
    }
    
    for (i in (k+1):(2*k))
    {
    mu[i] <- beta[2] + W[i]
    Y[i] ~ dnorm(mu[i], taue2)
    }
    
    for (i in (2*k+1):(3*k))
    {
    mu[i] <- beta[3] + W[i]
    Y[i] ~ dnorm(mu[i], taue3)
    }
    for (i in (3*k+1):(4*k))
    {
    mu[i] <-beta[4] + W[i]
    Y[i] ~ dnorm(mu[i], taue4)
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
    beta[1:4] ~ dmnorm(rep(0,4), (0.001*I1))
    }
    ", fill = TRUE)
sink()

model.data <- list(k = n, index1 = index1, I = diag(n), I1 = diag(4), Minc = Minc, ns = dni, cn = c(0, cni), udnei = udnei, Y = Y)
model.inits <- rep(list(list(rho1 = 0.1, rho2 = 0.1, rho3 = 0.1, tau1 = 1, tau2 = 1, tau3 = 1, eta021 = 1, 
                             eta121 = 1, eta031 = 1, eta131 = 1, eta032 = 1, eta132 = 1, eta041 = 1, 
                             eta141 = 1, eta042 = 1, eta142 = 1, eta043 = 1, eta143 = 1, taue1 = 1, taue2 = 1,
                             taue3 = 1, taue4 = 1, beta = rep(0, 4), W1 = rep(0.1, n), 
                             W2 = rep(0.1, n), W3 = rep(0.1, n), W4 = rep(0.1, n))),2)
model.param <- c("beta", "rho1", "rho2", "rho3", "rho4", "tau1", "tau2", "tau3", "tau4", "eta021", "eta121",
                 "eta031", "eta131", "eta032", "eta132","eta041", "eta141",
                 "eta042", "eta142", "eta043", "eta143", "vare1", "vare2", "vare3", "vare4", "W")
set.seed(123)
result_nocov <- jags(model.data, model.inits, model.param, "DAGAR_nocov.txt",
                     n.chains = 2, n.iter = 30000,n.burnin = 15000, n.thin = 1)
tod2=Sys.time()

saveRDS(result_nocov, "result_nocov.rds")

mcmc_select <- result_nocov
jags_result1 = mcmc_select$BUGSoutput$sims.array[,1,c(233:261)]
jags_result2 = mcmc_select$BUGSoutput$sims.array[,2,c(233:261)]
jags_result = rbind(jags_result1, jags_result2)
#mcmc_list[[9]]$BUGSoutput$summary[c(233:281,514:525),]
W_mcmc1 =  mcmc_select$BUGSoutput$sims.array[,1,1:232]
W_mcmc2 =  mcmc_select$BUGSoutput$sims.array[,2,1:232]
W_mcmc = rbind(W_mcmc1, W_mcmc2)
colnames(jags_result) = c("beta1", "beta2","beta3","beta4", "deviance", "eta0_21", "eta0_31", "eta0_32", "eta0_41",
                          "eta0_42", "eta0_43", "eta1_21", "eta1_31", "eta1_32", "eta1_41", "eta1_42", "eta1_43", "rho1", "rho2", "rho3", "rho4",
                          "tausq1", "tausq2", "tausq3", "tausq4","sigmasq1","sigmasq2","sigmasq3","sigmasq4")
estimate1 <- round(t(apply(jags_result, 2, mysummary)),2)

#Table for coefficients
table_est = paste(estimate1[,1], " (", estimate1[,4], ", ", estimate1[,5], ")", sep="")
table1 <- data.frame(matrix(table_est[1:4], ncol = 4, byrow=F))
table2 <- data.frame(matrix(table_est[18:29], ncol = 4, byrow=T))

table <- rbind(table1, table2)
colnames(table) <- c("Larynx", "Esophagus", "Lung", "Colorectum")
write.csv(table, "table_est_new.csv")

#Histogram for etas
pdf("eta_final_new.pdf", height = 10, width = 10)
par(mfrow=c(3,4))
hist(jags_result[,42],xlab = "eta_021, esophagus | larynx",ylab = "", main="")
hist(jags_result[,48],xlab = "eta_121, esophagus | larynx",ylab = "", main="")
hist(jags_result[,43],xlab = "eta_031, lung | larynx",ylab = "", main="")
hist(jags_result[,49],xlab = "eta_131, lung | larynx",ylab = "", main="")
hist(jags_result[,45],xlab = "eta_041, colorectum | larynx",ylab = "", main="")
hist(jags_result[,51],xlab = "eta_141, colorectum | larynx",ylab = "", main="")
hist(jags_result[,44],xlab = "eta_032, lung | esophagus",ylab = "", main="")
hist(jags_result[,50],xlab = "eta_132, lung | esophagus",ylab = "", main="")
hist(jags_result[,46],xlab = "eta_042, colorectum | esophagus",ylab = "", main="")
hist(jags_result[,52],xlab = "eta_142, colorectum | esophagus",ylab = "", main="")
hist(jags_result[,47],xlab = "eta_043, colorectum | lung",ylab = "", main="")
hist(jags_result[,53],xlab = "eta_143, colorectum | lung",ylab = "", main="")
dev.off()


################## Correlation in each county ###############
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

############# Random effects ####################
W_mean = apply(W_mcmc, 2, mean)
W_mean1 = W_mean[1:58]
ca.poly$W_mean1 = W_mean1[order(final_perm)]
quantile(W_mean1)
brks_fit1 = c(-0.6, -0.3, -0.1, 0, 0.1, 0.25, 0.76)
color.pallete = rev(brewer.pal(6,"RdBu"))
class.W1 = classIntervals(var=ca.poly$W_mean1, n=6, style="fixed", 
                          fixedBreaks=brks_fit1, dataPrecision=4)


W_mean2 = W_mean[59:116]
ca.poly$W_mean2 = W_mean2[order(final_perm)]
brks_fit2 = c(-2.5, -0.5, 0, 0.35, 1, 2, 3)
class.W2 = classIntervals(var=ca.poly$W_mean2, n=6, style="fixed", 
                          fixedBreaks=brks_fit2, dataPrecision=4)
color.code.W2 = findColours(class.W2, color.pallete)

W_mean3 = W_mean[117:174]
ca.poly$W_mean3 = W_mean3[order(final_perm)]
brks_fit3 = c(-21, -2, 2, 5, 10, 20, 31)
class.W3 = classIntervals(var=ca.poly$W_mean3, n=6, style="fixed", 
                          fixedBreaks=brks_fit3, dataPrecision=4)
color.code.W3 = findColours(class.W3, color.pallete)

W_mean4 = W_mean[118:232]
ca.poly$W_mean4 = W_mean4[order(final_perm)]
brks_fit4 = c(-21, -1, 1, 3, 10, 20, 31)
class.W4 = classIntervals(var=ca.poly$W_mean4, n=6, style="fixed", 
                          fixedBreaks=brks_fit4, dataPrecision=4)
color.code.W4 = findColours(class.W4, color.pallete)



# Plots for random effects
pdf("/Users/Leiwen/Box Sync/Research/data_new/random_nocovariates.pdf", height = 10, width = 10)
par(mfrow=c(2,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(ca.poly, col = color.code.W1)
leg.txt = c("-0.6 - -0.3", "-0.3 - -0.1", "-0.1 - 0", "0 - 0.1", "0.1 - 0.25", "0.25 - 0.76") 
legend("bottomleft", title="Larynx cancer", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W2)
leg.txt = c("-2.5 - -0.5", "-0.5 - 0", "0 - 0.35", "0.35 - 1", "1 - 2", "2 - 3") 
legend("bottomleft", title="Esophageal cancer", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W3)
leg.txt3 = c("-21 - -2", "-2 - 2", "2 - 5", "5 - 10", "10 - 20", "20 - 31")
legend("bottomleft", title="Lung cancer",legend=leg.txt3, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W4)
leg.txt4 = c("-21 - -1", "-1 - 1", "1 - 3", "3 - 10", "10 - 20", "20 - 31")
legend("bottomleft",title="Colorectum cancer", legend=leg.txt4, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

dev.off()
