sim_clot_log <- function(Nsim, clotting, truepar, seed) {
    clot_ML <- glm(conc ~ lot*log(u),data = clotting, family = Gamma("log"),
                   method = "brglmFit", type = "ML",maxit=1000,epsilon = 1e-8)
    d2afuns <- enrichwith::enrich(clot_ML$family, with = c("d1afun", "d2afun",
                                                           "d3afun", "d1variance"))$d2afun
    weights <- as.vector(clot_ML$prior.weights)

    X <- model.matrix(clot_ML)
    p <- length(truepar)
    coef0 <- truepar[-p]
    disp <- truepar[p]
    simufun <- enrichwith::get_simulate_function(clot_ML)
    datasim <- simufun(coefficients = coef0, dispersion = disp, n = Nsim,
                       seed = seed)

    mleg <- mle <- mlebc <- mlebr <- mlembr <- mlemixed <- matrix(NA,Nsim,p)
    seg <- se <- sebc <- sebr <- sembr <- semixed <- matrix(NA,Nsim,p)
    convmleg <- convmle <- convbc <- convbr <- convmbr <-
        convmixed <- rep(0,Nsim)
    indconvmleg <- indconvmle <- indconvbc <- indconvbr <-
        indconvmbr <- indconvmixed <- NULL

    errmleg <- errmle <- errbc <- errbr <- errmbr <- errmixed <-
        vector("list", Nsim)
    gradmleg <- gradmle <- gradbc <- gradbr <- gradmbr <- gradmixed <-
        vector("list", Nsim)
    inderrmleg <- inderrmle <- inderrbc <- inderrbr <- inderrmbr <-
        inderrmixed <-NULL

    for(i in seq.int(Nsim)) {
        if(i%%100 ==0) cat("Iteration: ", i,"\n")
        y <- datasim[,i]
        mlg <- try(glm(y ~ lot*log(u), data = clotting,
                       family=Gamma(link = "log"),
                       maxit=1000,epsilon = 1e-8),TRUE)
        if(inherits(mlg, "try-error")) {
            errmleg[[i]] <- mlg
            inderrmleg <- c(inderrmleg,i)
        }
        if(!inherits(mlg, "try-error")) {
            convmleg[i] <- sum(mlg$converged)
            if (sum(mlg$converged)==0) indconvmleg <- c(indconvmleg ,i)
            mlg <- summary(mlg)
            mleg[i,] <- c(mlg$coefficients[,1],mlg$dispersion)
            zetas <- -weights/mlg$dispersion
            seg[i,] <- c(mlg$coefficients[,2],
                         sqrt((2*mlg$dispersion^4)/sum(weights^2*d2afuns(zetas),
                                                       na.rm = TRUE)))
        }

        ml <- try(glm(y ~ lot*log(u), data = clotting,family=Gamma(link = "log"),
                      method = "brglmFit", type = "ML", maxit=1000,
                      epsilon = 1e-8),TRUE)
        if(inherits(ml, "try-error")) {
            errmle[[i]] <- ml
            inderrmle <- c(inderrmle,i)
        }
        if(!inherits(ml, "try-error")) {
            convmle[i] <- sum(ml$converged)
            if (sum(ml$converged)==0) indconvmle <- c(indconvmle ,i)
            gradmle[[i]]<- ml$grad
            ml <- summary(ml)
            mle[i,] <- c(ml$coefficients[,1],ml$dispersion)
            zetas <- -weights/ml$dispersion
            se[i,] <- c(ml$coefficients[,2],
                        sqrt((2*ml$dispersion^4)/sum(weights^2*d2afuns(zetas),
                                                     na.rm = TRUE)))
        }

        br <- try(glm(y ~ lot*log(u), data = clotting,family=Gamma(link = "log"),
                      method = "brglmFit", type = "AS_mean",maxit=1000,
                      epsilon = 1e-8),TRUE)
        if(inherits(br, "try-error")) {
            errbr[[i]] <- br
            inderrbr <- c(inderrbr,i)
        }
        if(!inherits(br, "try-error")) {
            convbr[i] <- sum(br$converged)
            if (sum(br$converged)==0) indconvbr <- c(indconvbr ,i)
            gradbr[[i]]<- br$grad
            br <- summary(br)
            mlebr[i,] <- c(br$coefficients[,1],br$dispersion)
            zetas <- -weights/br$dispersion
            sebr[i,] <- c(br$coefficients[,2],
                          sqrt((2*br$dispersion^4)/sum(weights^2*d2afuns(zetas),
                                                       na.rm = TRUE)))
        }

        bc <- try(glm(y ~ lot*log(u), data = clotting,family=Gamma(link = "log"),
                      method = "brglmFit", type = "correction",maxit=1000,
                      epsilon = 1e-8),TRUE)
        if(inherits(bc, "try-error")) {
            errbc[[i]] <- bc
            inderrbc <- c(inderrbc,i)
        }
        if(!inherits(bc, "try-error")) {
            convbc[i] <- sum(bc$converged)
            if (sum(bc$converged)==0) indconvbc <- c(indconvbc ,i)
            gradbc[[i]]<- bc$grad
            bc <- summary(bc)
            mlebc[i,] <- c(bc$coefficients[,1],bc$dispersion)
            zetas <- -weights/bc$dispersion
            sebc[i,] <- c(bc$coefficients[,2],
                          sqrt((2*bc$dispersion^4)/sum(weights^2*d2afuns(zetas),
                                                       na.rm = TRUE)))
        }

        mbr <- try(glm(y ~ lot*log(u), data = clotting,family=Gamma(link = "log"),
                       method = "brglmFit", type = "AS_median",maxit=1000,
                       epsilon = 1e-8),TRUE)
        if(inherits(mbr, "try-error")) {
            errmbr[[i]] <- mbr
            inderrmbr <- c(inderrmbr,i)
        }
        if(!inherits(mbr, "try-error")) {
            convmbr[i] <- sum(mbr$converged)
            if (sum(mbr$converged)==0) indconvmbr <- c(indconvmbr ,i)
            gradmbr[[i]]<- mbr$grad
            mbr <- summary(mbr)
            mlembr[i,] <- c(mbr$coefficients[,1],mbr$dispersion)
            zetas <- -weights/mbr$dispersion
            sembr[i,] <- c(mbr$coefficients[,2],
                           sqrt((2*mbr$dispersion^4)/sum(weights^2*d2afuns(zetas),
                                                         na.rm = TRUE)))
        }

        mixed <- try(glm(y ~ lot*log(u), data = clotting,
                         family=Gamma(link = "log"),method = "brglmFit",
                         type = "AS_mixed",maxit=1000,epsilon = 1e-8),TRUE)
        if(inherits(mixed, "try-error")) {
            errmixed[[i]] <- mixed
            inderrmixed <- c(inderrmixed,i)
        }
        if(!inherits(mixed, "try-error")) {
            convmixed[i] <- sum(mixed$converged)
            if (sum(mixed$converged)==0) indconvmixed <- c(indconvmixed ,i)
            gradmixed[[i]]<- mixed$grad
            mixed <- summary(mixed)
            mlemixed[i,] <- c(mixed$coefficients[,1],mixed$dispersion)
            zetas <- -weights/mixed$dispersion
            semixed[i,] <- c(mixed$coefficients[,2],sqrt((2*mixed$dispersion^4)/
                                                         sum(weights^2*d2afuns(zetas),
                                                             na.rm = TRUE)))
        }
    }
    list(mleg = mleg, mle = mle, mlebc = mlebc, mlebr = mlebr, mlembr = mlembr,
         mlemixed = mlemixed, seg = seg,se = se, sebc = sebc, sebr = sebr, sembr = sembr,
         semixed = semixed, convmleg = convmleg, convmle = convmle, convbc = convbc,
         convmixed = convmixed, convbr = convbr, convmbr = convmbr,
         indconvmleg = indconvmleg, indconvmle = indconvmle, indconvbc = indconvbc,
         indconvmixed = indconvmixed, indconvbr = indconvbr, indconvmbr = indconvmbr,
         gradmleg = gradmleg,gradmle = gradmle, gradbc = gradbc, gradbr = gradbr,
         gradmbr = gradmbr, gradmixed = gradmixed,inderrmleg = inderrmleg,
         inderrmle = inderrmle, inderrbc = inderrbc, inderrbr = inderrbr,
         inderrmbr = inderrmbr, inderrmixed = inderrmixed, errmleg = errmleg,
         errmle = errmle, errbc = errbc, errbr = errbr, errmbr = errmbr,
         errmixed = errmixed, datasim = datasim, X = X, Nsim = Nsim, truepar = truepar,
         dataset = clotting, seed = seed)
}


summarySim <- function(object) {
    convmleg <- object$convmleg
    convmle <- object$convmle
    convbc <- object$convbc
    convbr <- object$convbr
    convmbr <- object$convmbr
    convmixed <- object$convmixed
    mleg <- object$mleg
    mle <- object$mle
    mlebc <- object$mlebc
    mlebr <- object$mlebr
    mlembr <- object$mlembr
    mlemixed <- object$mlemixed
    seg <- object$seg
    se <- object$se
    sebc <- object$sebc
    sebr <- object$sebr
    sembr <- object$sembr
    semixed <- object$semixed
    truepar <- object$truepar
    p <- length(truepar)
    Nsim <- object$Nsim
    X <- object$X

    indmleg <- which(convmleg==0)
    indmle <- which(convmle==0)
    indbc <- which(convbc==0)
    indbr <- which(convbr==0)
    indmbr <- which(convmbr==0)
    indmixed <- which(convmixed==0)
    ind <- unique(c(indmleg,indmle,indbc,indbr,indmbr,indmixed))
    if (length(ind)==0){
        ind <-TRUE
    }else{
        ind <- -c(ind)
    }

    no_convergence<-rbind(mleglm=length(indmleg),mle=length(indmle),
                          mlebc=length(indbc),mlebr=length(indbr),
                          mlembr=length(indmbr),mlemixed=length(indmixed))

    trueparmat <- matrix(rep(truepar,Nsim),Nsim,p,byrow = T)
    z0.975 <- matrix(qnorm(0.975),Nsim,p)
    z0.025 <- matrix(qnorm(0.025),Nsim,p)
    z0.975 <- matrix(qnorm(0.975),Nsim,p)
    z0.975 <- z0.975[ind,]
    z0.025 <- z0.025[ind,]
    trueparmat <- trueparmat[ind,]
    mleg <- mleg[ind,]
    mle <- mle[ind,]
    mlebc <- mlebc[ind,]
    mlebr <- mlebr[ind,]
    mlembr <- mlembr[ind,]
    mlemixed <- mlemixed[ind,]
    seg <- seg[ind,]
    se <- se[ind,]
    sebc <- sebc[ind,]
    sebr <- sebr[ind,]
    sembr <- sembr[ind,]
    semixed <- semixed[ind,]

    errormleg <- apply(mleg,2,function(x) sqrt(var(x)/length(x)))
    errormle <- apply(mle,2,function(x) sqrt(var(x)/length(x)))
    errorbc <- apply(mlebc,2,function(x) sqrt(var(x)/length(x)))
    errorbr <- apply(mlebr,2,function(x) sqrt(var(x)/length(x)))
    errormbr <- apply(mlembr,2,function(x) sqrt(var(x)/length(x)))
    errormixed <- apply(mlemixed,2,function(x) sqrt(var(x)/length(x)))
    error <- rbind(errormleg,errormle,errorbc,errorbr,errormbr,errormixed)
    colnames(error) <- c(colnames(X),"dispersion")

    sdmleg <- apply(mleg,2,function(x) sqrt(var(x)))
    sdmle <- apply(mle,2,function(x) sqrt(var(x)))
    sdbc <- apply(mlebc,2,function(x) sqrt(var(x)))
    sdbr <- apply(mlebr,2,function(x) sqrt(var(x)))
    sdmbr <- apply(mlembr,2,function(x) sqrt(var(x)))
    sdmixed <- apply(mlemixed,2,function(x) sqrt(var(x)))
    sd <- rbind(sdmleg,sdmle,sdbc,sdbr,sdmbr,sdmixed)
    colnames(error) <- c(colnames(X),"dispersion")

    bmlg <- colMeans(mleg-trueparmat)
    bml <- colMeans(mle-trueparmat)
    bbc <- colMeans(mlebc-trueparmat)
    bbr <- colMeans(mlebr-trueparmat)
    bmbr <- colMeans(mlembr-trueparmat)
    bmixed <- colMeans(mlemixed-trueparmat)
    b <- round(rbind(bmlg,bml,bbc,bbr,bmbr,bmixed),6)
    dimnames(b) <-list(c("BIAS_mlg","BIAS_ml","BIAS_bc","BIAS_br",
                         "BIAS_mbr","BIAS_mixed"),c(colnames(X),"dispersion"))

    pumlg <- colMeans(mleg<trueparmat)
    puml <- colMeans(mle<trueparmat)
    pubc <- colMeans(mlebc<trueparmat)
    pubr <- colMeans(mlebr<trueparmat)
    pumbr <- colMeans(mlembr<trueparmat)
    pumixed <- colMeans(mlemixed<trueparmat)
    pu <- round(rbind(pumlg,puml,pubc,pubr,pumbr,pumixed),6)
    dimnames(pu) <- list(c("PU_mlg","PU_ml","PU_bc","PU_br",
                           "PU_mbr","PU_mixed"),c(colnames(X),"dispersion"))

    rmsemlg <- sqrt(colMeans((mleg-trueparmat)^2))
    rmseml <- sqrt(colMeans((mle-trueparmat)^2))
    rmsebc <- sqrt(colMeans((mlebc-trueparmat)^2))
    rmsebr <- sqrt(colMeans((mlebr-trueparmat)^2))
    rmsembr <- sqrt(colMeans((mlembr-trueparmat)^2))
    rmsemixed <- sqrt(colMeans((mlemixed-trueparmat)^2))
    rmse <- rbind(rmsemlg,rmseml,rmsebc,rmsebr,rmsembr,rmsemixed)
    dimnames(rmse) <- list(c("RMSE_mlg","RMSE_ml","RMSE_bc","RMSE_br",
                             "RMSE_mbr","RMSE_mixed"),c(colnames(X),"dispersion"))

    maemlg <- colMeans(abs(mleg-trueparmat))
    maeml <- colMeans(abs(mle-trueparmat))
    maebc <- colMeans(abs(mlebc-trueparmat))
    maebr <- colMeans(abs(mlebr-trueparmat))
    maembr <- colMeans(abs(mlembr-trueparmat))
    maemixed <- colMeans(abs(mlemixed-trueparmat))
    mae <- rbind(maemlg,maeml,maebc,maebr,maembr,maemixed)
    dimnames(mae) <- list(c("MAE_mlg","MAE_ml","MAE_bc","MAE_br",
                            "MAE_mbr","MAE_mixed"),c(colnames(X),"dispersion"))

    covmlg_b <- colMeans((abs((mleg-trueparmat)/seg) <= z0.975))
    covml_b <- colMeans((abs((mle-trueparmat)/se) <= z0.975))
    covbc_b <- colMeans((abs((mlebc-trueparmat)/sebc) <= z0.975))
    covbr_b <- colMeans((abs((mlebr-trueparmat)/sebr) <= z0.975))
    covmbr_b <- colMeans((abs((mlembr-trueparmat)/sembr) <= z0.975))
    covmixed_b <- colMeans((abs((mlemixed-trueparmat)/semixed) <= z0.975))
    cov_b <- round(rbind(covmlg_b,covml_b,covbc_b,covbr_b,
                         covmbr_b,covmixed_b),6)
    dimnames(cov_b) <- list(c("COV_mlg","COV_ml","COV_bc","COV_br",
                              "COV_mbr","COV_mixed"),c(colnames(X),"dispersion"))

    covmlg_l <- colMeans(((mleg-trueparmat)/seg <= z0.025))
    covml_l <- colMeans(((mle-trueparmat)/se <= z0.025))
    covbc_l <- colMeans(((mlebc-trueparmat)/sebc <= z0.025))
    covbr_l <- colMeans(((mlebr-trueparmat)/sebr <= z0.025))
    covmbr_l <- colMeans(((mlembr-trueparmat)/sembr <= z0.025))
    covmixed_l <- colMeans(((mlemixed-trueparmat)/semixed <= z0.025))
    cov_l <- rbind(covmlg_l,covml_l,covbc_l,covbr_l,covmbr_l,covmixed_l)
    dimnames(cov_l) <- list(c("COV_mlg","COV_ml","COV_bc","COV_br",
                              "COV_mbr","COV_mixed"),c(colnames(X),"dispersion"))

    covmlg_r <- colMeans(((mleg-trueparmat)/seg >= z0.975))
    covml_r <- colMeans(((mle-trueparmat)/se >= z0.975))
    covbc_r <- colMeans(((mlebc-trueparmat)/sebc >= z0.975))
    covbr_r <- colMeans(((mlebr-trueparmat)/sebr >= z0.975))
    covmbr_r <- colMeans(((mlembr-trueparmat)/sembr >= z0.975))
    covmixed_r <- colMeans(((mlemixed-trueparmat)/semixed >= z0.975))
    cov_r <- rbind(covmlg_r,covml_r,covbc_r,covbr_r,covmbr_r,covmixed_r)
    dimnames(cov_r) <- list(c("COV_mlg","COV_ml","COV_bc","COV_br",
                              "COV_mbr","COV_mixed"),c(colnames(X),"dispersion"))
    list(no_convergence=no_convergence,BIAS=round(b,6),
         PU=round(pu,6),RMSE=round(rmse,6),
         MAE=round(mae,6),COV_b=round(cov_b,6),
         COV_l=round(cov_l,6),COV_r=round(cov_r,6),
         sim_error=round(error,6),sample_size=nrow(X),
         SD=round(sd,6), RIB=round(b^2/sd^2,6),
         prior_replication=Nsim,
         samples_excluded=ifelse(sum(ind)==1,0,length(ind)))
}
