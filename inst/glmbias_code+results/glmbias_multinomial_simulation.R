######################
## Simulation study ##
######################

library("plyr")
library("doMC")
library("brglm2")
library("nnet")

## Use 18 cores for the simulation
registerDoMC(18)

all_ml <- brmultinom(foodchoice ~ size + lake , weights = freq,
                       data = alligators,
                       ref = 1,
                       type = "ML")
true_coefs <- coef(all_ml)

agresti_contrasts <- list(lake = contr.treatment(levels(alligators$lake),
                          base = 4),
                          size = contr.treatment(levels(alligators$size),
                          base = 2))

## Contrasts matrix
mat <- cbind(c(1, 1, 1, 0, 0),
             c(0, -1, 0, 0, 0),
             c(0, 0, -1, 1, 0),
             c(0, 0, -1, 0, 1),
             c(0, 0, -1, 0, 0))

## coef(all_ml_agresti) %*% mat - true_coefs


## Model specific and not at all efficient!
simulate_alligators <- function(fitted, total_factor = 1) {
    nams <- apply(alligators[c("lake", "gender", "size")], 1,
    paste0, collapse = "|")
    totals <- round(total_factor * tapply(alligators$freq, nams, sum))
    ## this is ordered
    oo <- order(nams)
    al <- alligators[oo, ]
    fitted <- fitted[oo, ]
    inds <- which(!duplicated(al[, c("lake", "gender", "size")]))
    covariate_settings <- al[inds, c("lake", "gender", "size")]
    freq <- sapply(seq.int(totals), function(j) rmultinom(1, totals[j],
                   fitted[inds[j], ]))
    out <- data.frame(foodchoice = colnames(fitted),
           covariate_settings[rep(rownames(covariate_settings),
           each = nrow(freq)),] ,
    freq = c(freq))
    out$foodchoice <- factor(out$foodchoice,
                             levels = levels(alligators$foodchoice))
    out
}

nsimu <- 10000
factors <- c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5)

res <- as.list(numeric(length(factors)))

for (k in seq_along(factors)) {
    set.seed(123)
    seeds <- sample(seq.int(nsimu*1000), nsimu, replace = FALSE)
    datasets <- llply(seeds, function(seed) {
        set.seed(seed)
        dat <- simulate_alligators(fitted(all_ml),
        total_factor = factors[k])
    }, .parallel = TRUE)
    res[[k]] <- llply(seq_along(datasets), function(j) {
        current_data <- datasets[[j]]
        ml_fit <- multinom(foodchoice ~ size + lake , weights = freq,
                           data = current_data, trace = FALSE)
        ml_c <- coef(ml_fit)
        ml_se <- summary(ml_fit)$standard.errors
        br_mean_fit <- brmultinom(foodchoice ~ size + lake ,
                                  weights = freq,
                                  data = current_data, ref = "Fish",
                                  trace = FALSE,
                                  type = "AS_mean")
        br_mean_c <- coef(br_mean_fit)
        br_mean_se <- summary(br_mean_fit)$standard.errors
        br_median_fit <- brmultinom(foodchoice ~ size + lake ,
                                    weights = freq,
                                    data = current_data, ref = "Fish",
                                    trace = FALSE,
                                    type = "AS_median")
        br_median_c <- coef(br_median_fit)
        br_median_se <- summary(br_median_fit)$standard.errors
        br_median_gamma_fit <- brmultinom(foodchoice ~ size + lake ,
                                          weights = freq,
                                          data = current_data,
                                          ref = "Fish",
                                          contrasts = agresti_contrasts,
                                          trace = FALSE,
                                          type = "AS_median")
        br_median_gamma_c <- coef(br_median_gamma_fit) %*% mat
        colnames(br_median_gamma_c) <- colnames(ml_c)
        if (j %% 10 == 0) cat(factors[k], "|", j, "\n")
        list(ml = ml_c,
             br_mean = br_mean_c,
             br_median = br_median_c,
             br_median_gamma = br_median_gamma_c,
             ml_se = ml_se, br_mean_se = br_mean_se,
             br_median_se = br_median_se)
    }, .parallel = TRUE)
}

ml_coefs <- lapply(res, function(x) lapply(x, "[[", "ml"))
ml_ses <- lapply(res, function(x) lapply(x, "[[", "ml_se"))
ml_infinite <- lapply(ml_ses, function(x)
                      sapply(x, function(y) any(is.na(y)) | any(y > 100,
                      na.rm = TRUE)))
ml_prob_separation <- sapply(ml_infinite, mean)
ml_bias <- lapply(ml_coefs, function(x) 100*(Reduce("+", x)/nsimu -
                  true_coefs)/true_coefs)
ml_pu <- lapply(ml_coefs, function(x) 100*Reduce("+", lapply(x, "<",
                true_coefs))/nsimu)
## Mean bias reduction
br_mean_coefs <- lapply(res, function(x) lapply(x, "[[", "br_mean"))
br_mean_bias <- lapply(br_mean_coefs, function(x) 100*(Reduce("+", x)/
                       nsimu - true_coefs)/true_coefs)
br_mean_pu <- lapply(br_mean_coefs, function(x) 100*Reduce("+",
                     lapply(x, "<", true_coefs))/nsimu)
## Median bias reduction
br_median_coefs <- lapply(res, function(x)
                          lapply(x, "[[", "br_median"))
br_median_bias <- lapply(br_median_coefs,
                         function(x) 100*(Reduce("+", x)/
                         nsimu - true_coefs)/true_coefs)
br_median_pu <- lapply(br_median_coefs, function(x) 100*Reduce("+",
                       lapply(x, "<", true_coefs))/nsimu)
## Median bias reduction
br_median_gamma_coefs <- lapply(res, function(x)
                                lapply(x, "[[", "br_median_gamma"))
br_median_gamma_bias <- lapply(br_median_gamma_coefs,
                               function(x) 100*(Reduce("+", x)/nsimu -
                                true_coefs)/true_coefs)
br_median_gamma_pu <- lapply(br_median_gamma_coefs,
                             function(x) 100*Reduce("+",
                             lapply(x, "<", true_coefs))/nsimu)
## Prepare for plotting
ml_bias_matrix <- do.call("rbind", ml_bias)
categories <- rownames(ml_bias_matrix)
ml_results <- stack(as.data.frame(ml_bias_matrix))
ml_results$category <- categories
ml_results$factor <- rep(factors, each = 4)
ml_results$method <- "ml"
ml_results$pu <- stack(as.data.frame(do.call("rbind", ml_pu)))$values
br_mean_bias_matrix <- do.call("rbind", br_mean_bias)
categories <- rownames(br_mean_bias_matrix)
br_mean_results <- stack(as.data.frame(br_mean_bias_matrix))
br_mean_results$category <- categories
br_mean_results$factor <- rep(factors, each = 4)
br_mean_results$method <- "mean BR"
br_mean_results$pu <- stack(as.data.frame(do.call("rbind", br_mean_pu)))$values
br_median_bias_matrix <- do.call("rbind", br_median_bias)
categories <- rownames(br_median_bias_matrix)
br_median_results <- stack(as.data.frame(br_median_bias_matrix))
br_median_results$category <- categories
br_median_results$factor <- rep(factors, each = 4)
br_median_results$method <- "median BR"
br_median_results$pu <- stack(as.data.frame(do.call("rbind", br_median_pu)))$values
br_median_gamma_bias_matrix <- do.call("rbind", br_median_gamma_bias)
categories <- rownames(br_median_gamma_bias_matrix)
br_median_gamma_results <- stack(as.data.frame(br_median_gamma_bias_matrix))
br_median_gamma_results$category <- categories
br_median_gamma_results$factor <- rep(factors, each = 4)
br_median_gamma_results$method <- "median BR gamma"
br_median_gamma_results$pu <- stack(as.data.frame(do.call("rbind", br_median_gamma_pu)))$values


save(br_mean_results, br_median_results, br_median_gamma_results, file = "alligator_simulation_results.rda")
