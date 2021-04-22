#simulations are run in two parts parallelly with 50 datasets in each part
#results are combined for 85 datasets

#delete the 15 datasets with abnormal results for GMCAR
index = c(2,6,8,17,19,25,26,29,37,44,54,73,76,89,91)

setwd("/Users/Leiwen/Box Sync/Research/simulation/GP/DAGARl_sw")

Dd = cbind(readRDS("D2_usa1.rds")[1:50], readRDS("D2_usa2.rds")[51:100])[-index]

w1hatd = rbind(readRDS("w1hat_usa1.rds")[1:50,], readRDS("w1hat_usa2.rds")[51:100,])[-index,]
w1lowerd = rbind(readRDS("w1lower_usa1.rds")[1:50,], readRDS("w1lower_usa2.rds")[51:100,])[-index,]
w1upperd = rbind(readRDS("w1upper_usa1.rds")[1:50,], readRDS("w1upper_usa2.rds")[51:100,])[-index,]

w2hatd = rbind(readRDS("w2hat_usa1.rds")[1:50,], readRDS("w2hat_usa2.rds")[51:100,])[-index,]
w2lowerd = rbind(readRDS("w2lower_usa1.rds")[1:50,], readRDS("w2lower_usa2.rds")[51:100,])[-index,]
w2upperd = rbind(readRDS("w2upper_usa1.rds")[1:50,], readRDS("w2upper_usa2.rds")[51:100,])[-index,]

W1d = rbind(readRDS("W1_usa1.rds")[1:50,], readRDS("W1_usa2.rds")[51:100,])[-index,]
W2d = rbind(readRDS("W2_usa1.rds")[1:50,], readRDS("W2_usa2.rds")[51:100,])[-index,]

WAIC_d = c(unlist(readRDS("WAIC_usa1.rds"))[1:50], unlist(readRDS("WAIC_usa2.rds"))[51:100])[-index]
KL_d = rbind(readRDS("KL_usa1.rds")[1:50,],readRDS("KL_usa2.rds")[51:100,])[-index,]


n=48
q=2


# AMSE
amse_matrix_dagar= NULL

for (t in 1:85){
  amse_matrix_dagar[t] = sum((W1d[t,]-w1hatd[t,])^2) + sum((W2d[t,]-w2hatd[t,])^2)
}

amse_W_dagarm = sum(amse_matrix_dagar)/(85*n*q)


#SE(AMSE)
se_matrix_dagar= NULL

for (t in 1:85){
  se_matrix_dagar[t] = sum(((W1d[t,]-w1hatd[t,])^2-amse_W_dagarm)^2) + sum(((W2d[t,]-w2hatd[t,])^2-amse_W_dagarm)^2)
}

se_W_dagarm = sqrt(sum(se_matrix_dagar)/(85*n*q)/(85*n*q-1))



setwd("/Users/Leiwen/Box Sync/Research/simulation/GP/DAGARl_ne")

D_ne = cbind(readRDS("D2_usa1.rds")[1:50], readRDS("D2_usa2.rds")[51:100])[-index]

w1hatd_ne = rbind(readRDS("w1hat_usa1.rds")[1:50,], readRDS("w1hat_usa2.rds")[51:100,])[-index,]
w1lowerd_ne = rbind(readRDS("w1lower_usa1.rds")[1:50,], readRDS("w1lower_usa2.rds")[51:100,])[-index,]
w1upperd_ne = rbind(readRDS("w1upper_usa1.rds")[1:50,], readRDS("w1upper_usa2.rds")[51:100,])[-index,]

w2hatd_ne = rbind(readRDS("w2hat_usa1.rds")[1:50,], readRDS("w2hat_usa2.rds")[51:100,])[-index,]
w2lowerd_ne = rbind(readRDS("w2lower_usa1.rds")[1:50,], readRDS("w2lower_usa2.rds")[51:100,])[-index,]
w2upperd_ne = rbind(readRDS("w2upper_usa1.rds")[1:50,], readRDS("w2upper_usa2.rds")[51:100,])[-index,]

W1d_ne = rbind(readRDS("W1_usa1.rds")[1:50,], readRDS("W1_usa2.rds")[51:100,])[-index,]
W2d_ne = rbind(readRDS("W2_usa1.rds")[1:50,], readRDS("W2_usa2.rds")[51:100,])[-index,]

WAIC_dne = c(readRDS("WAIC_usa1.rds")[1:50], readRDS("WAIC_usa2.rds")[51:100])[-index]

KL_dne = rbind(readRDS("KL_usa1.rds")[1:50,],readRDS("KL_usa2.rds")[51:100,])[-index,]

amse_matrix_dagar_ne= NULL

for (t in 1:85){
  amse_matrix_dagar_ne[t] = sum((W1d_ne[t,]-w1hatd_ne[t,])^2) + sum((W2d_ne[t,]-w2hatd_ne[t,])^2)
}

amse_W_dagarm_ne = sum(amse_matrix_dagar_ne)/(85*n*q)


#SE(AMSE)
se_matrix_dagar_ne= NULL

for (t in 1:85){
  se_matrix_dagar_ne[t] = sum(((W1d_ne[t,]-w1hatd_ne[t,])^2-amse_W_dagarm_ne)^2) + sum(((W2d_ne[t,]-w2hatd_ne[t,])^2-amse_W_dagarm_ne)^2)
}

se_W_dagarm_ne = sqrt(sum(se_matrix_dagar_ne)/(85*n*q)/(85*n*q-1))

setwd("/Users/Leiwen/Box Sync/Research/simulation/GP/DAGARl_nw")

D_nw = cbind(readRDS("D2_usa1.rds")[1:50], readRDS("D2_usa2.rds")[51:100])[-index]

w1hatd_nw = rbind(readRDS("w1hat_usa1.rds")[1:50,], readRDS("w1hat_usa2.rds")[51:100,])[-index,]
w1lowerd_nw = rbind(readRDS("w1lower_usa1.rds")[1:50,], readRDS("w1lower_usa2.rds")[51:100,])[-index,]
w1upperd_nw = rbind(readRDS("w1upper_usa1.rds")[1:50,], readRDS("w1upper_usa2.rds")[51:100,])[-index,]

w2hatd_nw = rbind(readRDS("w2hat_usa1.rds")[1:50,], readRDS("w2hat_usa2.rds")[51:100,])[-index,]
w2lowerd_nw = rbind(readRDS("w2lower_usa1.rds")[1:50,], readRDS("w2lower_usa2.rds")[51:100,])[-index,]
w2upperd_nw = rbind(readRDS("w2upper_usa1.rds")[1:50,], readRDS("w2upper_usa2.rds")[51:100,])[-index,]

W1d_nw = rbind(readRDS("W1_usa1.rds")[1:50,], readRDS("W1_usa2.rds")[51:100,])[-index,]
W2d_nw = rbind(readRDS("W2_usa1.rds")[1:50,], readRDS("W2_usa2.rds")[51:100,])[-index,]

WAIC_dnw = c(readRDS("WAIC_usa1.rds")[1:50], readRDS("WAIC_usa2.rds")[51:100])[-index]

KL_dnw = rbind(readRDS("KL_usa1.rds")[1:50,],readRDS("KL_usa2.rds")[51:100,])[-index,]

amse_matrix_dagar_nw= NULL

for (t in 1:85){
  amse_matrix_dagar_nw[t] = sum((W1d_nw[t,]-w1hatd_nw[t,])^2) + sum((W2d_nw[t,]-w2hatd_nw[t,])^2)
}

amse_W_dagarm_nw = sum(amse_matrix_dagar_nw)/(85*n*q)


#SE(AMSE)
se_matrix_dagar_nw= NULL

for (t in 1:85){
  se_matrix_dagar_nw[t] = sum(((W1d_nw[t,]-w1hatd_nw[t,])^2-amse_W_dagarm_nw)^2) + sum(((W2d_nw[t,]-w2hatd_nw[t,])^2-amse_W_dagarm_nw)^2)
}

se_W_dagarm_nw = sqrt(sum(se_matrix_dagar_nw)/(85*n*q)/(85*n*q-1))

setwd("/Users/Leiwen/Box Sync/Research/simulation/GP/DAGARl_se")

D_se = cbind(readRDS("D2_usa1.rds")[1:50], readRDS("D2_usa2.rds")[51:100])[-index]

w1hatd_se = rbind(readRDS("w1hat_usa1.rds")[1:50,], readRDS("w1hat_usa2.rds")[51:100,])[-index,]
w1lowerd_se = rbind(readRDS("w1lower_usa1.rds")[1:50,], readRDS("w1lower_usa2.rds")[51:100,])[-index,]
w1upperd_se = rbind(readRDS("w1upper_usa1.rds")[1:50,], readRDS("w1upper_usa2.rds")[51:100,])[-index,]

w2hatd_se = rbind(readRDS("w2hat_usa1.rds")[1:50,], readRDS("w2hat_usa2.rds")[51:100,])[-index,]
w2lowerd_se = rbind(readRDS("w2lower_usa1.rds")[1:50,], readRDS("w2lower_usa2.rds")[51:100,])[-index,]
w2upperd_se = rbind(readRDS("w2upper_usa1.rds")[1:50,], readRDS("w2upper_usa2.rds")[51:100,])[-index,]

W1d_se = rbind(readRDS("W1_usa1.rds")[1:50,], readRDS("W1_usa2.rds")[51:100,])[-index,]
W2d_se = rbind(readRDS("W2_usa1.rds")[1:50,], readRDS("W2_usa2.rds")[51:100,])[-index,]

WAIC_dse = c(readRDS("WAIC_usa1.rds")[1:50], readRDS("WAIC_usa2.rds")[51:100])[-index]

KL_dse = rbind(readRDS("KL_usa1.rds")[1:50,],readRDS("KL_usa2.rds")[51:100,])[-index,]

amse_matrix_dagar_se= NULL

for (t in 1:85){
  amse_matrix_dagar_se[t] = sum((W1d_se[t,]-w1hatd_se[t,])^2) + sum((W2d_se[t,]-w2hatd_se[t,])^2)
}

amse_W_dagarm_se = sum(amse_matrix_dagar_se)/(85*n*q)


#SE(AMSE)
se_matrix_dagar_se= NULL

for (t in 1:85){
  se_matrix_dagar_se[t] = sum(((W1d_se[t,]-w1hatd_se[t,])^2-amse_W_dagarm_se)^2) + sum(((W2d_se[t,]-w2hatd_se[t,])^2-amse_W_dagarm_se)^2)
}

se_W_dagarm_se = sqrt(sum(se_matrix_dagar_se)/(85*n*q)/(85*n*q-1))


# D score
value = c(Dd, D_ne, D_nw, D_se)
Model = c(rep("Southwest",85), rep("Northeast",85), rep("Northwest",85), rep("Southeast",85))

df = data.frame(cbind(value, Model))
df$value = as.numeric(as.character(df$value))

library(plyr)
mu <- ddply(df, "Model", summarise, grp.mean=mean(value))

pdf("/Users/Leiwen/Box Sync/Research/data_new/D_perm.pdf", height = 5, width = 8)
ggplot(df, aes(x=value,color=Model, fill=Model)) +
  geom_density(alpha=0.4) + xlab("Score D") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))
dev.off()

#WAIC
value = c(WAIC_d, WAIC_dne, WAIC_dnw, WAIC_dse)

df = data.frame(cbind(value, Model))
df$value = as.numeric(as.character(df$value))

mu <- ddply(df, "Model", summarise, grp.mean=mean(value))

pdf("/Users/Leiwen/Box Sync/Research/data_new/WAIC_perm.pdf", height = 5, width = 8)
ggplot(df, aes(x=value,color=Model, fill=Model)) +
  geom_density(alpha=0.4) + xlab("WAIC") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))
dev.off()

mean(WAIC_d)
mean(WAIC_dne)

#KL Divergence
KL_d1 <- rowMeans(KL_d)
KL_dne1 <- rowMeans(KL_dne)
KL_dnw1 <- rowMeans(KL_dnw)
KL_dse1 <- rowMeans(KL_dse)

value = c(KL_d1, KL_dne1, KL_dnw1, KL_dse1)

df = data.frame(cbind(value, Model))
df$value = as.numeric(as.character(df$value))

mu <- ddply(df, "Model", summarise, grp.mean=mean(value))

pdf("/Users/Leiwen/Box Sync/Research/data_new/KL_perm.pdf", height = 5, width = 8)
ggplot(df, aes(x=value,color=Model, fill=Model)) +
  geom_density(alpha=0.4) + xlab("D_KL") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))
dev.off()
