
library(plyr)
library(ggplot2)

#Import simalution results for a certain correlation scenario
#Take low correlation senario as an example: we import the simulation data for low correlation senario in result_data.R
scenario = "low"
source("result_data.R")
summary.function = function(x){
  m = mean(x)
  l = quantile(x, 0.025)
  u = quantile(x, 0.975)
  return(c(m ,l,u))
}

#WAIC
summary.function(WAIC_d)
summary.function(WAIC_c)



n=48
q=2

########### MDAGAR ############
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

# Take low correlation senario as an example
eta_truem = c(0.05, 0.1)
# medium correlation: eta_truem = c(0.5, 0.3)
# high correlation: eta_truem = c(2.5, 0.5)

rho_truem = c(0.2, 0.8)
tau_truem = c(0.25, 0.25)
vare_truem = c(0.4, 0.4)
b_truem = c(1, 5, 2, 4, 5)

#Coverage probability for parameters
pc_dagarm = rep(0, length(cor_regiond))
peta_dagarm = rep(0, nrow(etahatd))
prho_dagarm = rep(0, nrow(rhohatd))
ptau_dagarm = rep(0, nrow(tauwhatd))
pvare_dagarm = rep(0, nrow(varehatd))
pb_dagarm = rep(0, nrow(bhatd))
for(t in 1:85){
  for(i in 1:length(cor_regiond)){
    if(corlowerd[t,i] <= cor_regiond[i] & corupperd[t,i] >= cor_regiond[i])
      pc_dagarm[i] = pc_dagarm[i] + 1
  }
  for(i in 1:nrow(etahatd)){
    if(etalowerd[i,t] <= eta_truem[i] & etaupperd[i,t] >= eta_truem[i])
      peta_dagarm[i] = peta_dagarm[i] + 1
  }
  for(i in 1:nrow(rhohatd)){
    if(rholowerd[i,t] <= rho_truem[i] & rhoupperd[i,t] >= rho_truem[i])
      prho_dagarm[i] = prho_dagarm[i] + 1
  }
  for(i in 1:nrow(tauwhatd)){
    if(tauwlowerd[i,t] <= tau_truem[i] & tauwupperd[i,t] >= tau_truem[i])
      ptau_dagarm[i] = ptau_dagarm[i] + 1
  }
  for(i in 1:nrow(varehatd)){
    if(varelowerd[i,t] <= vare_truem[i] & vareupperd[i,t] >= vare_truem[i])
      pvare_dagarm[i] = pvare_dagarm[i] + 1
  }
  for(i in 1:nrow(bhatd)){
    if(blowerd[i,t] <= b_truem[i] & bupperd[i,t] >= b_truem[i])
      pb_dagarm[i] = pb_dagarm[i] + 1
  }
}

pc_dagarm = pc_dagarm/85
pb_dagarm = pb_dagarm/85
peta_dagarm = peta_dagarm/85
prho_dagarm = prho_dagarm/85
ptau_dagarm = ptau_dagarm/85
pvare_dagarm = pvare_dagarm/85

################# GMCAR #######################
# AMSE
amse_matrix_car= NULL

for (t in 1:85){
  amse_matrix_car[t] = sum((W1c[t,]-w1hatc[t,])^2) + sum((W2c[t,]-w2hatc[t,])^2)
}

amse_W_carm = sum(amse_matrix_car)/(85*n*q)


#SE(AMSE)
se_matrix_car= NULL

for (t in 1:85){
  se_matrix_car[t] = sum(((W1c[t,]-w1hatc[t,])^2-amse_W_carm)^2) + sum(((W2c[t,]-w2hatc[t,])^2-amse_W_carm)^2)
}

se_W_carm = sqrt(sum(se_matrix_car)/(85*n*q)/(85*n*q-1))


#Coverage probability for parameters
pc_carm = rep(0, length(cor_regionc))
peta_carm = rep(0, nrow(etahatc))
prho_carm = rep(0, nrow(rhohatc))
ptau_carm = rep(0, nrow(tauwhatc))
pvare_carm = rep(0, nrow(varehatc))
pb_carm = rep(0, nrow(bhatc))
for(t in 1:85){
  for(i in 1:length(cor_regionc)){
    if(corlowerc[t,i] <= cor_regionc[i] & corupperc[t,i] >= cor_regionc[i])
      pc_carm[i] = pc_carm[i] + 1
  }
  for(i in 1:nrow(etahatc)){
    if(etalowerc[i,t] <= eta_truem[i] & etaupperc[i,t] >= eta_truem[i])
      peta_carm[i] = peta_carm[i] + 1
  }
  for(i in 1:nrow(rhohatc)){
    if(rholowerc[i,t] <= rho_truem[i] & rhoupperc[i,t] >= rho_truem[i])
      prho_carm[i] = prho_carm[i] + 1
  }
  for(i in 1:nrow(tauwhatc)){
    if(tauwlowerc[i,t] <= tau_truem[i] & tauwupperc[i,t] >= tau_truem[i])
      ptau_carm[i] = ptau_carm[i] + 1
  }
  for(i in 1:nrow(varehatc)){
    if(varelowerc[i,t] <= vare_truem[i] & vareupperc[i,t] >= vare_truem[i])
      pvare_carm[i] = pvare_carm[i] + 1
  }
  for(i in 1:nrow(bhatc)){
    if(blowerc[i,t] <= b_truem[i] & bupperc[i,t] >= b_truem[i])
      pb_carm[i] = pb_carm[i] + 1
  }
}

pc_carm = pc_carm/85
pb_carm = pb_carm/85
peta_carm = peta_carm/85
prho_carm = prho_carm/85
ptau_carm = ptau_carm/85
pvare_carm = pvare_carm/85

# Coverage probability for correlation between diseases in each state
pcover = c(pc_carm*100, pc_dagarm*100)
model = c( rep("GMCAR", n), rep("MDAGAR", n))
data("state.fips")
states=setdiff(unique(state.fips$abb),"DC")
plot_pc = data.frame(cbind(states,pcover, model))
colnames(plot_pc) = c("location", "value", "model")
plot_pc$value = as.numeric(as.character(plot_pc$value))

pdf(paste("cp_cor_", scenario, ".pdf", sep = ""), height = 5, width = 8)
ggplot(plot_pc, aes(x=location, y=value, group=model)) +
  #geom_line(aes(color=model, linetype = model))+
  coord_cartesian(ylim = c(0,100))+
  geom_point(aes(color=model))+ylab("Coverage probability (%)")+
  geom_hline(yintercept = 95, linetype = "dashed", color = "black")+
  scale_x_discrete(limits=states)+ xlab("States") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size=8)) 
dev.off()

# Kullback–Leibler divergence
pl = NULL
for(i in 1:85){
  pl[i] = sum(KL_d[i,]<=KL_c[i,])/15000
}

hist(pl)
quantile(pl, c(0.025,0.5,0.975))
plot(density(pl))

# D score
Dd = D1d + D2d
Dc = D1c + D2c

value = c(Dd, Dc)
Model = c(rep("MDAGAR",85), rep("GMCAR",85))

df = data.frame(cbind(value, Model))
df$value = as.numeric(as.character(df$value))

mu <- ddply(df, "Model", summarise, grp.mean=mean(value))

pdf(paste("D_", scenario, ".pdf", sep = ""), height = 5, width = 8)
ggplot(df, aes(x=value,color=Model, fill=Model)) +
  geom_density(alpha=0.4) + xlab("D") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")
dev.off()

# WAIC plot
WAIC_plot <- data.frame(cbind(c(WAIC_d, WAIC_c), c(rep("MDAGAR", length(WAIC_d)),
                                                   rep("GMCAR", length(WAIC_c)))))
WAIC_plot$X1 <- as.numeric(as.character(WAIC_plot$X1))
colnames(WAIC_plot) <- c("WAIC", "Model")

mu <- ddply(WAIC_plot, "Model", summarise, grp.mean=mean(WAIC))

pdf(paste("WAIC_", scenario, ".pdf", sep = ""), height = 5, width = 8)
ggplot(WAIC_plot, aes(x=WAIC, color=Model,fill=Model)) +
  geom_density(alpha=0.4) + xlab("WAIC")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
             linetype="dashed")
dev.off()


# 2-dimensional histogram for estimates of W against true values
whatd = cbind(w1hatd, w2hatd)
Wtrued = cbind(W1d, W2d)
Wtruedvec = as.vector(t(Wtrued))
wvecd = as.vector(t(whatd))

Wd_tab = data.frame(cbind(Wtruedvec, wvecd))
pdf(paste("W_dagar_", scenario, ".pdf", sep = ""), height = 5, width = 8)
ggplot(Wd_tab, aes(x=Wtruedvec, y=wvecd) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_abline(slope=1, intercept = 0) +
  xlab("W") + ylab("Estimate of W")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))
dev.off()

whatc = cbind(w1hatc, w2hatc)
Wtruec = cbind(W1c, W2c)
Wtruecvec = as.vector(t(Wtruec))
wvecc = as.vector(t(whatc))

Wc_tab = data.frame(cbind(Wtruecvec, wvecc))
pdf(paste("W_car_", scenario, ".pdf", sep = ""), height = 5, width = 8)
ggplot(Wc_tab, aes(x=Wtruecvec, y=wvecc) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  geom_abline(slope=1, intercept = 0) +
  xlab("W") + ylab("Estimate of W")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))
dev.off()

