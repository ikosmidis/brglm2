## simulation study: clotting dataset

library("brglm2")
source("glmbias_simulation_functions.R")

clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))

clot_ML <- glm(conc ~ lot*log(u), data = clotting, family = Gamma(link="log"),
               method = "brglmFit", type = "ML", maxit=1000, epsilon = 1e-8)

truepar <- c(coef(clot_ML),clot_ML$dispersion)
clotting_simulation_results <- sim_clot_log(10000, clotting, truepar, 123)
clotting_simulation_results <- summarySim(clotting_simulation_results)

bias <- clotting_simulation_results$BIAS
sd <- clotting_simulation_results$SD
rib <- clotting_simulation_results$RIB
pu <- clotting_simulation_results$PU
mae <- clotting_simulation_results$MAE
cov <- clotting_simulation_results$COV_b
mleg <-  cbind(bias[1, ], sd[1, ], rib[1, ], pu[1, ], mae[1, ], cov[1, ])
mle <-  cbind(bias[2, ], sd[2, ], rib[2, ], pu[2, ], mae[2, ], cov[2, ])
meanBR <- cbind(bias[4, ], sd[4, ], rib[4, ], pu[4, ], mae[4, ], cov[4, ])
medianBR <- cbind(bias[5, ], sd[5, ], rib[5, ], pu[5, ], mae[5, ], cov[5, ])
meanmixed <- cbind(bias[6, ], sd[6, ], rib[6, ], pu[6, ], mae[6, ], cov[6, ])

save(mle, meanBR, medianBR, meanmixed, file = "clotting_simulation_results.rda")
