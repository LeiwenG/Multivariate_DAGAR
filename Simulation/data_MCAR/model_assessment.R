#Import simulation results for model fitting parameters from MCAR, MDAGR and GMCAR models (Output files from MCAR.R, MDAGAR_o1.R, MDAGAR_o2.R, GMCAR_o1.R, GMCAR_o2.R)
# MDAGAR_o1.R and GMCAR_o1.R are for the first order, MDAGAR_o2.R and GMCAR_o2.R are for the second order

#Import results for MDAGAR (set directory to output files for MDAGAR)
w1hatd = rbind(readRDS("w1hat_usa1.rds")[51:100,], readRDS("w1hat_usa2.rds")[51:100,])
w1lowerd = rbind(readRDS("w1lower_usa1.rds")[51:100,], readRDS("w1lower_usa2.rds")[51:100,])
w1upperd = rbind(readRDS("w1upper_usa1.rds")[51:100,], readRDS("w1upper_usa2.rds")[51:100,])

w2hatd = rbind(readRDS("w2hat_usa1.rds")[51:100,], readRDS("w2hat_usa2.rds")[51:100,])
w2lowerd = rbind(readRDS("w2lower_usa1.rds")[51:100,], readRDS("w2lower_usa2.rds")[51:100,])
w2upperd = rbind(readRDS("w2upper_usa1.rds")[51:100,], readRDS("w2upper_usa2.rds")[51:100,])


W1d = rbind(readRDS("W1_usa1.rds")[51:100,], readRDS("W2_usa2.rds")[51:100,])
W2d = rbind(readRDS("W2_usa1.rds")[51:100,], rreadRDS("W1_usa2.rds")[91:100,])

WAIC_d = c(unlist(readRDS("WAIC_usa1.rds"))[51:100], unlist(readRDS("WAIC_usa2.rds"))[51:100])
Dd = c(readRDS("D2_usa1.rds")[51:100], readRDS("D2_usa2.rds")[51:100])
KL_d = rbind(readRDS("KL_usa1.rds")[51:100,], readRDS("KL_usa2.rds")[51:100,])

summary.function = function(x){
  m = mean(x)
  l = quantile(x, 0.025)
  u = quantile(x, 0.975)
  return(c(m ,l,u))
}




n=48
q=2
# AMSE
amse_matrix_dagar= NULL

for (t in 1:50){
  amse_matrix_dagar[t] = sum((W1d[t,]-w1hatd[t,])^2) + sum((W2d[t,]-w2hatd[t,])^2)
}

amse_W_dagarm = sum(amse_matrix_dagar, na.rm = TRUE)/(50*n*q)


#SE(AMSE)
se_matrix_dagar= NULL

for (t in 1:50){
  se_matrix_dagar[t] = sum(((W1d[t,]-w1hatd[t,])^2-amse_W_dagarm)^2) + sum(((W2d[t,]-w2hatd[t,])^2-amse_W_dagarm)^2)
}

se_W_dagarm = sqrt(sum(se_matrix_dagar, na.rm = TRUE)/(50*n*q)/(50*n*q-1))

amse_matrix_dagar= NULL

for (t in 51:100){
  amse_matrix_dagar[t] = sum((W1d[t,]-w1hatd[t,])^2) + sum((W2d[t,]-w2hatd[t,])^2)
}

amse_W_dagarm = sum(amse_matrix_dagar, na.rm = TRUE)/(50*n*q)


#SE(AMSE)
se_matrix_dagar= NULL

for (t in 51:100){
  se_matrix_dagar[t] = sum(((W1d[t,]-w1hatd[t,])^2-amse_W_dagarm)^2) + sum(((W2d[t,]-w2hatd[t,])^2-amse_W_dagarm)^2)
}

se_W_dagarm = sqrt(sum(se_matrix_dagar, na.rm = TRUE)/(50*n*q)/(50*n*q-1))


#Import results for GMCAR (set directory to output files for GMCAR)

w1hatc = rbind(readRDS("w1hat_usa1.rds")[51:100,1:48], readRDS("w1hat_usa2.rds")[51:100,1:48])
w1lowerc = rbind(readRDS("w1lower_usa1.rds")[51:100,1:48], readRDS("w1lower_usa2.rds")[51:100,1:48])
w1upperc = rbind(readRDS("w1upper_usa1.rds")[51:100,1:48], readRDS("w1upper_usa2.rds")[51:100,1:48])

w2hatc = rbind(readRDS("w2hat_usa1.rds")[51:100,1:48], readRDS("w2hat_usa2.rds")[51:100,1:48])
w2lowerc = rbind(readRDS("w2lower_usa1.rds")[51:100,1:48], readRDS("w2lower_usa2.rds")[51:100,1:48])
w2upperc = rbind(readRDS("w2upper_usa1.rds")[51:100,1:48], readRDS("w2upper_usa2.rds")[51:100,1:48])


W1c = rbind(readRDS("W1_usa1.rds")[51:100,], readRDS("W2_usa2.rds")[51:100,])
W2c = rbind(readRDS("W2_usa1.rds")[51:100,], readRDS("W1_usa2.rds")[51:100,])

WAIC_c = c(unlist(readRDS("WAIC_usa1.rds"))[51:100], unlist(readRDS("WAIC_usa2.rds"))[51:100])
Dc = c(readRDS("D2_usa1.rds")[51:100], readRDS("D2_usa2.rds")[51:100])
KL_c = rbind(readRDS("KL_usa1.rds")[51:100,],readRDS("KL_usa2.rds")[51:100,])

# AMSE
amse_matrix_car= NULL

for (t in 1:50){
  amse_matrix_car[t] = sum((W1c[t,]-w1hatc[t,])^2) + sum((W2c[t,]-w2hatc[t,])^2)
}

amse_W_carm = sum(amse_matrix_car, na.rm = TRUE)/(50*n*q)


#SE(AMSE)
se_matrix_car= NULL

for (t in 1:50){
  se_matrix_car[t] = sum(((W1c[t,]-w1hatc[t,])^2-amse_W_carm)^2) + sum(((W2c[t,]-w2hatc[t,])^2-amse_W_carm)^2)
}

se_W_carm = sqrt(sum(se_matrix_car, na.rm = TRUE)/(50*n*q)/(50*n*q-1))

# AMSE
amse_matrix_car= NULL

for (t in 51:100){
  amse_matrix_car[t] = sum((W1c[t,]-w1hatc[t,])^2) + sum((W2c[t,]-w2hatc[t,])^2)
}

amse_W_carm = sum(amse_matrix_car, na.rm = TRUE)/(50*n*q)


#SE(AMSE)
se_matrix_car= NULL

for (t in 51:100){
  se_matrix_car[t] = sum(((W1c[t,]-w1hatc[t,])^2-amse_W_carm)^2) + sum(((W2c[t,]-w2hatc[t,])^2-amse_W_carm)^2)
}

se_W_carm = sqrt(sum(se_matrix_car, na.rm = TRUE)/(50*n*q)/(50*n*q-1))


#Import results for MCAR (set directory to output files for MCAR)
w1hatjc = readRDS("w1hat_usa1.rds")[51:100,1:48]
w1lowerjc = readRDS("w1lower_usa1.rds")[51:100,1:48]
w1upperjc = readRDS("w1upper_usa1.rds")[51:100,1:48]

w2hatjc = readRDS("w2hat_usa1.rds")[51:100,1:48]
w2lowerjc = readRDS("w2lower_usa1.rds")[51:100,1:48]
w2upperjc = readRDS("w2upper_usa1.rds")[51:100,1:48]


W1jc = readRDS("W1_usa1.rds")[51:100,]
W2jc = readRDS("W2_usa1.rds")[51:100,]

WAIC_jc = unlist(readRDS("WAIC_usa1.rds"))[51:100]
Djc = readRDS("D2_usa1.rds")[51:100]
KL_jc = readRDS("KL_usa1.rds")[51:100,]

n=48
q=2
# AMSE
amse_matrix_jc= NULL

for (t in 1:50){
  amse_matrix_jc[t] = sum((W1jc[t,]-w1hatjc[t,])^2) + sum((W2jc[t,]-w2hatjc[t,])^2)
}

amse_W_jcm = sum(amse_matrix_jc, na.rm = TRUE)/(50*n*q)


#SE(AMSE)
se_matrix_jc= NULL

for (t in 1:50){
  se_matrix_jc[t] = sum(((W1jc[t,]-w1hatjc[t,])^2-amse_W_jcm)^2) + sum(((W2jc[t,]-w2hatjc[t,])^2-amse_W_jcm)^2)
}

se_W_jcm = sqrt(sum(se_matrix_jc)/(50*n*q)/(50*n*q-1))



##### generate plots to compare five models

# WAIC plot
WAIC_plot <- data.frame(cbind(c(WAIC_d, WAIC_c, WAIC_jc), c(rep("MDAGAR1", length(WAIC_d)/2), 
                                                   rep("MDAGAR2", length(WAIC_d)/2), 
                                                   rep("GMCAR1", length(WAIC_c)/2),
                                                   rep("GMCAR2", length(WAIC_c)/2),
                                                   rep("MCAR", length(WAIC_jc)))))
WAIC_plot$X1 <- as.numeric(as.character(WAIC_plot$X1))
colnames(WAIC_plot) <- c("WAIC", "Model")
library(plyr)
mu <- ddply(WAIC_plot, "Model", summarise, grp.mean=mean(WAIC))

pdf("WAIC_median.pdf", height = 5, width = 8)
ggplot(WAIC_plot, aes(x=WAIC, color=Model,fill=Model)) +
  geom_density(alpha=0.4, adjust = 1.2) + xlab("WAIC")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")
dev.off()

# D score plot
value = c(Dd, Dc, Djc)
Model = c(rep("MDAGAR1", length(Dd)/2),
            rep("MDAGAR2", length(Dd)/2),
            rep("GMCAR1", length(Dc)/2),
            rep("GMCAR2", length(Dc)/2),
            rep("MCAR", length(Djc)))

df = data.frame(cbind(value, Model))
df$value = as.numeric(as.character(df$value))

mu <- ddply(df, "Model", summarise, grp.mean=mean(value))

pdf("D_medium.pdf", height = 5, width = 8)

ggplot(df, aes(x=value,color=Model, fill=Model)) +
  geom_density(alpha=0.4, adjust = 1.5) + xlab("D") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")
dev.off()

# KL plot
KL_value = c(rowMeans(KL_d), rowMeans(KL_c),rowMeans(KL_jc))
KL_plot <- data.frame(KL_value, c(rep("MDAGAR1", nrow(KL_d)/2),
                                  rep("MDAGAR2", nrow(KL_d)/2),
                                  rep("GMCAR1", nrow(KL_c)/2),
                                  rep("GMCAR2", nrow(KL_c)/2),
                                  rep("MCAR", nrow(KL_jc))))



colnames(KL_plot) <- c("KL", "Model")
mu <- ddply(KL_plot, "Model", summarise, grp.mean=mean(KL), q1 = quantile(KL, 0.025), q3 = quantile(KL, 0.975))

pdf("KL_median.pdf", height = 5, width = 8)
ggplot(KL_plot, aes(x=KL, color = Model, fill=Model)) +
  geom_density(alpha=0.4, adjust = 1.2) + xlab("KL")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")
dev.off()
