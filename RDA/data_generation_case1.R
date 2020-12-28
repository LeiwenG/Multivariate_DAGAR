library(maps)
ca.county = map("county","california", fill=TRUE, plot=FALSE)

library(tidyr)
library(readr)

library(spdep)
library(maptools)
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
county_attribute$state = readr::parse_number(as.character(county_attribute$State_county))
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

rate_list = list(rate_lung1, rate_esophagus1, rate_larynx1, rate_colrect1)


## Create adjacency matrix and neighbor info
ca.neighbors = poly2nb(ca.poly)
n=length(ca.neighbors)

Adj=sapply(ca.neighbors,function(x,n) {v=rep(0,n);v[x]=1;v},n)
colnames(Adj)=county.ID
ca.coord = coordinates(ca.poly)
ca.latrange=round(quantile(ca.coord[,2],c(0.25,0.75)))
ca.albersproj=mapproject(ca.coord[,1],ca.coord[,2],projection = "albers",param=ca.latrange)

perm=order(ca.albersproj$x-ca.albersproj$y)

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

q = 4
