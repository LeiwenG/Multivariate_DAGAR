library(rjmcmc)
library(mnormt)

setwd("/Users/Leiwen/Dropbox/Github/Multivariate_DAGAR/RDA")
mcmc_list = readRDS("mcmc_list_new.rds")

# Draw samples
jags_result = list()
draw_sample = NULL
for(i in 1:24){
  jags_result1 = mcmc_list[[i]]$BUGSoutput$sims.array[,1,c(233:285,518:529)]
  jags_result2 = mcmc_list[[i]]$BUGSoutput$sims.array[,2,c(233:285,518:529)]
  jags_result[[i]] = rbind(jags_result1, jags_result2)
  draw_sample = cbind(draw_sample, jags_result[[i]])
}

draw = function(){draw_sample[sample(dim(draw_sample)[1], 1, replace=T),
                              -which(colnames(draw_sample) == "deviance")]}


gf = list()
ginvf = list()
for(i in 1:24){
  gf[[i]] = function(psi){ psi }
  ginvf[[i]] = function(theta){theta}
}

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


# Create data for 24 models
Ylist = list()
Xlist = list()

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
  
  Ylist[[i]] = c(Y1,Y2,Y3,Y4)
  Xlist[[i]] = as.matrix(bdiag(bdiag(X1, X2), bdiag(X3,X4)))
}


# Case 2: load prior and likelihood
load("p_prior_new.RData")
load("likelihood_new.RData")

# For case 1, load("p_prior.RData"), load("likelihood.RData")


#Compare with model 1
goals_post = list()
ch_ind = NULL
for(i in 2:24){
  print(i)
  goals_post[[i]] = rjmcmcpost(post.draw = list(draw,draw), 
                               g = list(gf[[1]],gf[[i]]),
                               ginv =list(ginvf[[1]],ginvf[[i]]), 
                               likelihood = list(L[[1]],L[[i]]),
                               param.prior = list(p.prior[[1]],p.prior[[i]]),
                               model.prior = rep(1/2,2), 
                               chainlength = 2000)
  if(which(goals_post[[i]]$result$`Posterior Model Probabilities` == max(goals_post[[i]]$result$`Posterior Model Probabilities`)) == 2){
    ch_ind = c(ch_ind, i)
  }
}

saveRDS(goals_post, "goals_post_new.rds")

#Find superior models compared with model 1
ch_ind = NULL
for(i in 2:24){
  print(i)
  if(which(goals_post[[i]]$result$`Posterior Model Probabilities` == max(goals_post[[i]]$result$`Posterior Model Probabilities`)) == 2){
    ch_ind = c(ch_ind, i)
  }
}


ind_prob = do.call(rbind, lapply(goals_post, function(x) c(x$result$`Posterior Model Probabilities`, x$result$`Bayes Factors`[1,2])))
ind_prob_data = data.frame(round(ind_prob,3))
colnames(ind_prob_data) = c("prob for model 1", "prob for compared model", "BF12")
ind_prob_data$model = 2:24
ind_prob1 = ind_prob[,2]

#Find the model with largest prob compared with model 1
max_model1 = which(ind_prob1==max(ind_prob1)) + 1
ch_ind1 = ch_ind[-(which(ch_ind == max_model1))]

#Compare superior models to the model with largest prob
goals_post1 = list()
ch_ind2= NULL
iter = 0
for(i in ch_ind1){
  iter = iter + 1
  print(iter)
  print(i)
  goals_post1[[iter]] = rjmcmcpost(post.draw = list(draw,draw), 
                                   g = list(gf[[max_model1]],gf[[i]]),
                                   ginv =list(ginvf[[max_model1]],ginvf[[i]]), 
                                   likelihood = list(L[[max_model1]],L[[i]]),
                                   param.prior = list(p.prior[[max_model1]],p.prior[[i]]),
                                   model.prior = rep(1/2,2), 
                                   chainlength = 2000)
  if(which(goals_post1[[iter]]$result$`Posterior Model Probabilities` == max(goals_post1[[iter]]$result$`Posterior Model Probabilities`)) == 2){
    ch_ind2 = c(ch_ind2, i)
  }
}

saveRDS(goals_post1, "goals_post1_new.rds")

#Find superior models to the model with largest prob
ch_ind2= NULL
iter = 0
for(i in ch_ind1){
  iter = iter + 1
  print(iter)
  print(i)
  
  if(which(goals_post1[[iter]]$result$`Posterior Model Probabilities` == max(goals_post1[[iter]]$result$`Posterior Model Probabilities`)) == 2){
    ch_ind2 = c(ch_ind2, i)
  }
}


ind_prob2 = data.frame(do.call(rbind, lapply(goals_post1, function(x) c(x$result$`Posterior Model Probabilities`, x$result$`Bayes Factors`[1,2]))))
ind_prob2_data = data.frame(round(ind_prob2,3))
colnames(ind_prob2_data) = c("prob for model 1", "prob for compared model", "BF12")

ind_prob2_data$model = ch_ind1

#Save tables with posterior prob and BF
write.csv(ind_prob2_data, "prob2_data_new.csv")
write.csv(ind_prob_data, "prob_data_new.csv")


