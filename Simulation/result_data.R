#Import simulation results for all parameters in each correlation scenario (Output files from MDAGAR_gibbs_metroplis.R and GMCAR_jags.R)

#delete the 15 datasets with abnormal results for GMCAR
index = c(2,6,8,17,19,25,26,29,37,44,54,73,76,89,91)

#simulations are run in two parts parallelly with 50 datasets in each part
#results are combined for 85 datasets
#MDAGAR
bhatd = cbind(readRDS("bhat_usa1.rds")[,1:50], readRDS("bhat_usa2.rds")[,51:100])[,-index]
blowerd = cbind(readRDS("blower_usa1.rds")[,1:50], readRDS("blower_usa2.rds")[,51:100])[,-index]
bupperd = cbind(readRDS("bupper_usa1.rds")[,1:50], readRDS("bupper_usa2.rds")[,51:100])[,-index]

varehatd = cbind(readRDS("varehat_usa1.rds")[,1:50], readRDS("varehat_usa2.rds")[,51:100])[,-index]
varelowerd = cbind(readRDS("varelower_usa1.rds")[,1:50], readRDS("varelower_usa2.rds")[,51:100])[,-index]
vareupperd = cbind(readRDS("vareupper_usa1.rds")[,1:50], readRDS("vareupper_usa2.rds")[,51:100])[,-index]

tauwhatd = cbind(readRDS("tauwhat_usa1.rds")[,1:50], readRDS("tauwhat_usa2.rds")[,51:100])[,-index]
tauwlowerd = cbind(readRDS("tauwlower_usa1.rds")[,1:50], readRDS("tauwlower_usa2.rds")[,51:100])[,-index]
tauwupperd = cbind(readRDS("tauwupper_usa1.rds")[,1:50], readRDS("tauwupper_usa2.rds")[,51:100])[,-index]

etahatd = cbind(readRDS("etahat_usa1.rds")[,1:50], readRDS("etahat_usa2.rds")[,51:100])[,-index]
etalowerd = cbind(readRDS("etalower_usa1.rds")[,1:50], readRDS("etalower_usa2.rds")[,51:100])[,-index]
etaupperd = cbind(readRDS("etaupper_usa1.rds")[,1:50], readRDS("etaupper_usa2.rds")[,51:100])[,-index]

rhohatd = cbind(readRDS("rhohat_usa1.rds")[,1:50], readRDS("rhohat_usa2.rds")[,51:100])[,-index]
rholowerd = cbind(readRDS("rholower_usa1.rds")[,1:50], readRDS("rholower_usa2.rds")[,51:100])[,-index]
rhoupperd = cbind(readRDS("rhoupper_usa1.rds")[,1:50], readRDS("rhoupper_usa2.rds")[,51:100])[,-index]

w1hatd = rbind(readRDS("w1hat_usa1.rds")[1:50,], readRDS("w1hat_usa2.rds")[51:100,])[-index,]
w1lowerd = rbind(readRDS("w1lower_usa1.rds")[1:50,], readRDS("w1lower_usa2.rds")[51:100,])[-index,]
w1upperd = rbind(readRDS("w1upper_usa1.rds")[1:50,], readRDS("w1upper_usa2.rds")[51:100,])[-index,]

w2hatd = rbind(readRDS("w2hat_usa1.rds")[1:50,], readRDS("w2hat_usa2.rds")[51:100,])[-index,]
w2lowerd = rbind(readRDS("w2lower_usa1.rds")[1:50,], readRDS("w2lower_usa2.rds")[51:100,])[-index,]
w2upperd = rbind(readRDS("w2upper_usa1.rds")[1:50,], readRDS("w2upper_usa2.rds")[51:100,])[-index,]

Covhatd = rbind(readRDS("Covhat_usa1.rds")[1:50,], readRDS("Covhat_usa2.rds")[51:100,])[-index,]
Covlowerd = rbind(readRDS("Covlower_usa1.rds")[1:50,], readRDS("Covlower_usa2.rds")[51:100,])[-index,]
Covupperd = rbind(readRDS("Covupper_usa1.rds")[1:50,], readRDS("Covupper_usa2.rds")[51:100,])[-index,]

corhatd = rbind(readRDS("corhat_usa1.rds")[1:50,], readRDS("corhat_usa2.rds")[51:100,])[-index,]
corlowerd = rbind(readRDS("corlower_usa1.rds")[1:50,], readRDS("corlower_usa2.rds")[51:100,])[-index,]
corupperd = rbind(readRDS("corupper_usa1.rds")[1:50,], readRDS("corupper_usa2.rds")[51:100,])[-index,]

Cov_regiond = readRDS("Cov_region.usa.rds")
cor_regiond = readRDS("cor_region.usa.rds")

W1d = rbind(readRDS("W1_usa1.rds")[1:50,], readRDS("W1_usa2.rds")[51:100,])[-index,]
W2d = rbind(readRDS("W2_usa1.rds")[1:50,], readRDS("W2_usa2.rds")[51:100,])[-index,]

WAIC_d = c(unlist(readRDS("WAIC_usa1.rds"))[1:50], unlist(readRDS("WAIC_usa2.rds"))[51:100])[-index]

G1d = cbind(readRDS("G1_usa1.rds")[1:50], readRDS("G1_usa2.rds")[51:100])[-index]
P1d = cbind(readRDS("P1_usa1.rds")[1:50], readRDS("P1_usa2.rds")[51:100])[-index]
D1d = cbind(readRDS("D1_usa1.rds")[1:50], readRDS("D1_usa2.rds")[51:100])[-index]

G2d = cbind(readRDS("G2_usa1.rds")[1:50], readRDS("G2_usa2.rds")[51:100])[-index]
P2d = cbind(readRDS("P2_usa1.rds")[1:50], readRDS("P2_usa2.rds")[51:100])[-index]
D2d = cbind(readRDS("D2_usa1.rds")[1:50], readRDS("D2_usa2.rds")[51:100])[-index]

KL_d = rbind(readRDS("KL_usa1.rds")[1:50,],readRDS("KL_usa2.rds")[51:100,])[-index,]

#GMCAR
bhatc = cbind(readRDS("bhat_usa1.rds")[,1:50], readRDS("bhat_usa2.rds")[,51:100])[,-index]
blowerc = cbind(readRDS("blower_usa1.rds")[,1:50], readRDS("blower_usa2.rds")[,51:100])[,-index]
bupperc = cbind(readRDS("bupper_usa1.rds")[,1:50], readRDS("bupper_usa2.rds")[,51:100])[,-index]

varehatc = cbind(readRDS("varehat_usa1.rds")[,1:50], readRDS("varehat_usa2.rds")[,51:100])[,-index]
varelowerc = cbind(readRDS("varelower_usa1.rds")[,1:50], readRDS("varelower_usa2.rds")[,51:100])[,-index]
vareupperc = cbind(readRDS("vareupper_usa1.rds")[,1:50], readRDS("vareupper_usa2.rds")[,51:100])[,-index]

tauwhatc = cbind(readRDS("tauwhat_usa1.rds")[,1:50], readRDS("tauwhat_usa2.rds")[,51:100])[,-index]
tauwlowerc = cbind(readRDS("tauwlower_usa1.rds")[,1:50], readRDS("tauwlower_usa2.rds")[,51:100])[,-index]
tauwupperc = cbind(readRDS("tauwupper_usa1.rds")[,1:50], readRDS("tauwupper_usa2.rds")[,51:100])[,-index]

etahatc = cbind(readRDS("etahat_usa1.rds")[,1:50], readRDS("etahat_usa2.rds")[,51:100])[,-index]
etalowerc = cbind(readRDS("etalower_usa1.rds")[,1:50], readRDS("etalower_usa2.rds")[,51:100])[,-index]
etaupperc = cbind(readRDS("etaupper_usa1.rds")[,1:50], readRDS("etaupper_usa2.rds")[,51:100])[,-index]

rhohatc = cbind(readRDS("rhohat_usa1.rds")[,1:50], readRDS("rhohat_usa2.rds")[,51:100])[,-index]
rholowerc = cbind(readRDS("rholower_usa1.rds")[,1:50], readRDS("rholower_usa2.rds")[,51:100])[,-index]
rhoupperc = cbind(readRDS("rhoupper_usa1.rds")[,1:50], readRDS("rhoupper_usa2.rds")[,51:100])[,-index]

w1hatc = rbind(readRDS("w1hat_usa1.rds")[1:50,1:48], readRDS("w1hat_usa2.rds")[51:100,1:48])[-index,]
w1lowerc = rbind(readRDS("w1lower_usa1.rds")[1:50,1:48], readRDS("w1lower_usa2.rds")[51:100,1:48])[-index,]
w1upperc = rbind(readRDS("w1upper_usa1.rds")[1:50,1:48], readRDS("w1upper_usa2.rds")[51:100,1:48])[-index,]

w2hatc = rbind(readRDS("w1hat_usa1.rds")[1:50,49:96], readRDS("w1hat_usa2.rds")[51:100,49:96])[-index,]
w2lowerc = rbind(readRDS("w1lower_usa1.rds")[1:50,49:96], readRDS("w1lower_usa2.rds")[51:100,49:96])[-index,]
w2upperc = rbind(readRDS("w1upper_usa1.rds")[1:50,49:96], readRDS("w1upper_usa2.rds")[51:100,49:96])[-index,]

Covhatc = rbind(readRDS("Covhat_usa1.rds")[1:50,], readRDS("Covhat_usa2.rds")[51:100,])[-index,]
Covlowerc = rbind(readRDS("Covlower_usa1.rds")[1:50,], readRDS("Covlower_usa2.rds")[51:100,])[-index,]
Covupperc = rbind(readRDS("Covupper_usa1.rds")[1:50,], readRDS("Covupper_usa2.rds")[51:100,])[-index,]

corhatc = rbind(readRDS("corhat_usa1.rds")[1:50,], readRDS("corhat_usa2.rds")[51:100,])[-index,]
corlowerc = rbind(readRDS("corlower_usa1.rds")[1:50,], readRDS("corlower_usa2.rds")[51:100,])[-index,]
corupperc = rbind(readRDS("corupper_usa1.rds")[1:50,], readRDS("corupper_usa2.rds")[51:100,])[-index,]

Cov_regionc = readRDS("Cov_region.usa.rds")
cor_regionc = readRDS("cor_region.usa.rds")

W1c = rbind(readRDS("W1_usa1.rds")[1:50,], readRDS("W1_usa2.rds")[51:100,])[-index,]
W2c = rbind(readRDS("W2_usa1.rds")[1:50,], readRDS("W2_usa2.rds")[51:100,])[-index,]

WAIC_c = c(unlist(readRDS("WAIC_usa1.rds"))[1:50], unlist(readRDS("WAIC_usa2.rds"))[51:100])[-index]

G1c = cbind(readRDS("G1_usa1.rds")[1:50], readRDS("G1_usa2.rds")[51:100])[-index]
P1c = cbind(readRDS("P1_usa1.rds")[1:50], readRDS("P1_usa2.rds")[51:100])[-index]
D1c = cbind(readRDS("D1_usa1.rds")[1:50], readRDS("D1_usa2.rds")[51:100])[-index]

G2c = cbind(readRDS("G2_usa1.rds")[1:50], readRDS("G2_usa2.rds")[51:100])[-index]
P2c = cbind(readRDS("P2_usa1.rds")[1:50], readRDS("P2_usa2.rds")[51:100])[-index]
D2c = cbind(readRDS("D2_usa1.rds")[1:50], readRDS("D2_usa2.rds")[51:100])[-index]

KL_c = rbind(readRDS("KL_usa1.rds")[1:50,],readRDS("KL_usa2.rds")[51:100,])[-index,]


