setwd("Multivariate_DAGAR/Simulation/model_selection")

mp_model1 <- t(matrix(unlist(readRDS("mp_model1.rds")), nrow = 6))
best_model1 <- apply(mp_model1, 1, function(x) which.max(x))
prop.table(table(best_model1))

mp_model2 <- t(matrix(unlist(readRDS("mp_model2.rds")), nrow = 6))
best_model2 <- apply(mp_model2, 1, function(x) which.max(x))
prop.table(table(best_model2))

mp_model3 <- t(matrix(unlist(readRDS("mp_model3.rds")), nrow = 6))
best_model3 <- apply(mp_model3, 1, function(x) which.max(x))
prop.table(table(best_model3))

mp_model4 <- t(matrix(unlist(readRDS("mp_model4.rds")), nrow = 6))
best_model4 <- apply(mp_model4, 1, function(x) which.max(x))
prop.table(table(best_model4))

mp_model5 <- t(matrix(unlist(readRDS("mp_model5.rds")), nrow = 6))
best_model5 <- apply(mp_model5, 1, function(x) which.max(x))
prop.table(table(best_model5))

mp_model6 <- t(matrix(unlist(readRDS("mp_model6.rds")), nrow = 6))
best_model6 <- apply(mp_model6, 1, function(x) which.max(x))
prop.table(table(best_model6))
