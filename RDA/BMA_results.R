source("data_generation_case1.R")
library(classInt)
library(RColorBrewer)
#jags outputs for all models
mcmc_list1 = readRDS("mcmc_list_new.rds")
#model probabilities from bridge sampling
model.bridge = readRDS("model_bridge_jags.rds")
logml = rep(0, 24)
for(i in 1:24){
  logml[i] = model.bridge[[i]]$logml
}
post_prob = exp(logml)/sum(exp(logml))
round(post_prob, 4)

orderindex = c(1,2,3,4)
models = permutations(n=4,r=4,v=orderindex,repeats.allowed=F)

beta_model = data.frame(matrix(0, 40, 24))
W_model = data.frame(matrix(0, 232, 24))

for(i in 1:24){
  print(i)

  beta_result1 = mcmc_list1[[i]]$BUGSoutput$sims.array[,1,c(233:272)]
  beta_result2 = mcmc_list1[[i]]$BUGSoutput$sims.array[,2,c(233:272)]
  beta_est = data.frame(rbind(beta_result1, beta_result2))
  beta_mean = data.frame(apply(beta_est, 2, mean))
  beta_mean$cancer = rep(models[i,], each = 10)
  beta_mean1 = beta_mean[order(beta_mean$cancer),][,1]
  beta_model[,i] = beta_mean1
  
  W_mcmc1 =  mcmc_list1[[i]]$BUGSoutput$sims.array[,1,1:232]
  W_mcmc2 =  mcmc_list1[[i]]$BUGSoutput$sims.array[,2,1:232]
  W_mcmc = data.frame(rbind(W_mcmc1, W_mcmc2))
  W_mean = data.frame(apply(W_mcmc, 2, mean))
  W_mean$cancer = rep(models[i,], each = 58)
  W_mean1 = W_mean[order(W_mean$cancer),][,1]
  W_model[,i] = W_mean1
}

#BMA estimates for mean random effects and incidence rates
beta_wm = data.frame(t(apply(beta_model, 1, function(x) x * post_prob)))
beta_wmean = rowSums(beta_wm)
W_wm = data.frame(t(apply(W_model, 1, function(x) x * post_prob)))
W_wmean = rowSums(W_wm)


W_mean1 = W_wmean[1:58]
ca.poly$W_mean1 = W_mean1[order(final_perm)]
brks_fit1 = c(-20, -10, 0, 5, 10, 32)
color.pallete = rev(brewer.pal(5,"RdBu"))
class.W1 = classIntervals(var=ca.poly$W_mean1, n=5, style="fixed", 
                          fixedBreaks=brks_fit1, dataPrecision=4)
color.code.W1 = findColours(class.W1, color.pallete)


W_mean2 = W_wmean[59:116]
ca.poly$W_mean2 = W_mean2[order(final_perm)]
brks_fit2 = c(-1.2, -0.5, 0, 0.1, 0.5, 1)
class.W2 = classIntervals(var=ca.poly$W_mean2, n=5, style="fixed", 
                          fixedBreaks=brks_fit2, dataPrecision=4)
color.code.W2 = findColours(class.W2, color.pallete)

W_mean3 = W_wmean[117:174]
ca.poly$W_mean3 = W_mean3[order(final_perm)]
brks_fit3 = c(-1, -0.5, 0, 0.1, 0.5, 1.1)
class.W3 = classIntervals(var=ca.poly$W_mean3, n=5, style="fixed", 
                          fixedBreaks=brks_fit3, dataPrecision=4)
color.code.W3 = findColours(class.W3, color.pallete)

W_mean4 = W_wmean[175:232]
ca.poly$W_mean4 = W_mean4[order(final_perm)]
brks_fit4 = c(-5, -3, 0, 1, 3, 6)
class.W4 = classIntervals(var=ca.poly$W_mean4, n=5, style="fixed", 
                          fixedBreaks=brks_fit4, dataPrecision=4)
color.code.W4 = findColours(class.W4, color.pallete)

# Plots for random effects
pdf("random_covariates_ave.pdf", height = 10, width = 10)
par(mfrow=c(2,2), oma = c(0,0,4,0) + 0.1, mar = c(0,0,1,0) + 0.1)

plot(ca.poly, col = color.code.W1)
leg.txt = c("-20 - -10", "-10 - 0", "0 - 5", "5 - 10", "10 - 32")
legend("bottomleft", title="Lung cancer",legend=leg.txt, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W2)
leg.txt = c("-1.2 - -0.5", "-0.5 - 0", "0 - 0.1", "0.1 - 0.5", "0.5 - 1") 
legend("bottomleft", title="Esophageal cancer", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W3)
leg.txt = c("-1 - -0.5", "-0.5 - -0", "0 - 0.1", "0.1 - 0.5", "0.5 - 1.1") 
legend("bottomleft", title="Larynx cancer", legend=leg.txt, xpd = TRUE, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

plot(ca.poly, col = color.code.W4)
leg.txt = c("-5 - -3", "-3 - 0", "0 - 1", "1 - 3", "3 - 6")
legend("bottomleft",title="Colorectum cancer", legend=leg.txt, cex=1.25, bty="n", horiz = FALSE, 
       fill = color.pallete)

dev.off()

# Posterior mean incidence rates
Y1 = rate_list[[1]]$rate[final_perm]
Y2 = rate_list[[2]]$rate[final_perm]
Y3 = rate_list[[3]]$rate[final_perm]
Y4 = rate_list[[4]]$rate[final_perm]

X1 = as.matrix(cbind(1, rate_list[[1]][,6:14]))[final_perm,]
X2 = as.matrix(cbind(1, rate_list[[2]][,6:14]))[final_perm,]
X3 = as.matrix(cbind(1, rate_list[[3]][,6:14]))[final_perm,]
X4 = as.matrix(cbind(1, rate_list[[4]][,6:14]))[final_perm,]

Y = c(Y1,Y2,Y3,Y4)
X = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))


y.mean = as.vector(X %*% as.vector(beta_wmean) + W_wmean)

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

pdf("post_incidence_new_ave.pdf", height = 10, width = 10)
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


