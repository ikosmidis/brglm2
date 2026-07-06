library(brglm2)
library(microbenchmark)

## Universal variables
DENSE_TIMES <- 14L
SPARSE_TIMES <- 4L
TYPE <- "AS_median"

cat(sprintf("Benchmarking with type = '%s'\n", TYPE))

# Data generators

#  1. Logistic regression, AR-1 predictors (baseline) 
make_logistic_ar1 <- function(n, p, rho = 0.5, seed = 42) {
    set.seed(seed)
    p_pred <- p - 1L
    Sigma  <- rho^abs(outer(seq_len(p_pred), seq_len(p_pred), "-"))
    X      <- matrix(rnorm(n * p_pred), n, p_pred) %*% chol(Sigma)
    beta   <- c(0.3, rep(c(0.5, -0.5), length.out = p_pred))
    y      <- rbinom(n, 1, plogis(cbind(1, X) %*% beta))
    list(data = data.frame(y = y, X), family = binomial("logit"),
         label = "logistic AR-1")
}

#  2. Probit regression, block-correlated predictors 
# Simulates groups of correlated biomarkers (e.g. gene expression modules)
make_probit_blocks <- function(n, p, n_blocks = 5, rho = 0.7, seed = 43) {
    set.seed(seed)
    p_pred  <- p - 1L
    blk_sz  <- ceiling(p_pred / n_blocks)
    blocks  <- lapply(seq_len(n_blocks), function(b) {
        sz <- min(blk_sz, p_pred - (b - 1) * blk_sz)
        if (sz <= 0) return(NULL)
        matrix(rho, sz, sz) + diag(1 - rho, sz)
    })
    blocks  <- Filter(Negate(is.null), blocks)
    Sigma   <- as.matrix(Matrix::bdiag(blocks))[seq_len(p_pred), seq_len(p_pred)]
    X       <- matrix(rnorm(n * p_pred), n, p_pred) %*% chol(Sigma)
    beta    <- numeric(p_pred)
    for (b in seq_len(n_blocks)) beta[(b - 1) * blk_sz + 1] <- (-1)^b * 0.6
    beta    <- c(0.1, beta)
    y       <- rbinom(n, 1, pnorm(cbind(1, X) %*% beta))
    list(data = data.frame(y = y, X), family = binomial("probit"),
         label = "probit block-corr")
}

#  3. Poisson log-linear, sparse predictors 
# Simulates rare-event count data, e.g. adverse events in a clinical trial
make_poisson_sparse <- function(n, p, prop_signal = 0.2, seed = 44) {
    set.seed(seed)
    p_pred   <- p - 1L
    X        <- matrix(rnorm(n * p_pred), n, p_pred)
    n_signal <- max(1L, round(p_pred * prop_signal))
    beta     <- numeric(p_pred)
    beta[seq_len(n_signal)] <- rep(c(0.3, -0.3), length.out = n_signal)
    beta     <- c(-1.5, beta)
    y        <- rpois(n, exp(cbind(1, X) %*% beta))
    list(data = data.frame(y = y, X), family = poisson("log"),
         label = "Poisson sparse")
}

#  4. Logistic regression, high-dimensional medical (n approx 5p) 
# Toeplitz covariance, sparse signal -  mimics omics/GWAS logistic models
make_highdim_medical <- function(n, p, rho = 0.3, seed = 45) {
    set.seed(seed)
    p_pred   <- p - 1L
    Sigma    <- rho^abs(outer(seq_len(p_pred), seq_len(p_pred), "-"))
    X        <- scale(matrix(rnorm(n * p_pred), n, p_pred) %*% chol(Sigma))
    n_signal <- max(2L, round(p_pred * 0.05))
    beta     <- numeric(p_pred)
    beta[seq_len(n_signal)] <- rep(c(0.8, -0.8), length.out = n_signal)
    beta     <- c(0.0, beta)
    y        <- rbinom(n, 1, plogis(cbind(1, X) %*% beta))
    list(data = data.frame(y = y, X), family = binomial("logit"),
         label = "logistic high-dim")
}

#  5. Complementary log-log, survival-style binary outcome 
# Discrete-time survival model; cloglog link common in epidemiology
make_cloglog_survival <- function(n, p, rho = 0.4, seed = 46) {
    set.seed(seed)
    p_pred <- p - 1L
    Sigma  <- rho^abs(outer(seq_len(p_pred), seq_len(p_pred), "-"))
    X      <- matrix(rnorm(n * p_pred), n, p_pred) %*% chol(Sigma)
    beta   <- c(-1.5, rep(c(0.4, -0.3), length.out = p_pred))
    prob   <- 1 - exp(-exp(cbind(1, X) %*% beta))
    y      <- rbinom(n, 1, pmin(pmax(prob, 1e-8), 1 - 1e-8))
    list(data = data.frame(y = y, X), family = binomial("cloglog"),
         label = "cloglog survival")
}


# Scenarios: (generator, n, p)

scenarios <- list(
    list(gen = make_logistic_ar1,     n =  200, p =   5),
    list(gen = make_logistic_ar1,     n = 1000, p =  20),
    list(gen = make_logistic_ar1,     n = 4000, p =  50),
    list(gen = make_logistic_ar1,     n = 8000, p = 100),
    list(gen = make_probit_blocks,    n =  500, p =  10),
    list(gen = make_probit_blocks,    n = 2000, p =  30),
    list(gen = make_probit_blocks,    n = 5000, p =  60),
    list(gen = make_poisson_sparse,   n =  500, p =  10),
    list(gen = make_poisson_sparse,   n = 2000, p =  25),
    list(gen = make_highdim_medical,  n =  500, p = 100),
    list(gen = make_highdim_medical,  n = 1000, p = 200),
    list(gen = make_cloglog_survival, n =  800, p =  15),
    list(gen = make_cloglog_survival, n = 3000, p =  40)
)

cat(sprintf("Benchmarking dense with %d iterations per scenario\n", DENSE_TIMES+1))


# Run

results <- lapply(scenarios, function(sc) {
    sc_data <- sc$gen(sc$n, sc$p)
    n       <- sc$n;  p <- sc$p
    label   <- sc_data$label
    fam     <- sc_data$family
    dat     <- sc_data$data
    fmla    <- as.formula(paste("y ~", paste(names(dat)[-1], collapse = " + ")))

    cat(sprintf("\n %s  n=%d  p=%d \n", label, n, p))

    mb <- microbenchmark(
        original = glm(fmla, data = dat, family = fam, type = TYPE,
                       method = "brglmFit_original", maxit = 10000),
        new     = glm(fmla, data = dat, family = fam, type = TYPE,
                       method = "brglmFit",          maxit = 10000),
        times = DENSE_TIMES,
        unit  = "ms"
    )

    cat(sprintf("Coefficient fitting\n"))

    # Coefficient agreement
    fit_orig   <- glm(fmla, data = dat, family = fam, type = TYPE,
                      method = "brglmFit_original", maxit = 10000)
    fit_new   <- glm(fmla, data = dat, family = fam, type = TYPE,
                      method = "brglmFit", maxit = 10000)

    cat(sprintf("  max |coef diff| original vs new :   %.2e\n",
                max(abs(coef(fit_orig) - coef(fit_new)),   na.rm = TRUE)))

    attr(mb, "label") <- label
    attr(mb, "n")     <- n
    attr(mb, "p")     <- p
    mb
})


# Additional sparse synthetic data scenarios


#  1. One-hot categoricals, the most common real sparse structure 
# Each factor with k levels contributes k-1 binary columns, one 1 per row.
# With enough levels fill rate = 1/(k-1) which gets sparse quickly.
make_onehot_logistic <- function(n, n_factors, levels_per_factor, seed = 101) {
    set.seed(seed)
    factors <- lapply(seq_len(n_factors), function(i)
        factor(sample(seq_len(levels_per_factor), n, replace = TRUE),
               levels = seq_len(levels_per_factor)))   # <-- factor() with explicit levels
    df <- as.data.frame(setNames(factors, paste0("f", seq_len(n_factors))))
    X  <- model.matrix(~ ., data = df)
    p  <- ncol(X)
    beta <- c(0.2, rep(c(0.5, -0.5), length.out = p - 1))
    y  <- rbinom(n, 1, plogis(X %*% beta))
    fill <- mean(X != 0)
    cat(sprintf("\n  One-hot fill rate: %.4f  (n=%d, p=%d)", fill, n, p))
    list(X = X, y = y, family = binomial("logit"),
         label = sprintf("onehot %d factors x %d levels", n_factors, levels_per_factor))
}

#  2. Interaction dummies, sparser than main effects alone 
# Two-way interactions between binary indicators; 
# each interaction column has fill = p(A=1) * p(B=1)
make_interaction_logistic <- function(n, n_base, seed = 102) {
    set.seed(seed)
    # n_base binary predictors + all two-way interactions
    X_base <- matrix(rbinom(n * n_base, 1, 0.15), n, n_base)  # 15% fill base
    pairs  <- combn(n_base, 2)
    X_int  <- apply(pairs, 2, function(ij) X_base[, ij[1]] * X_base[, ij[2]])
    X    <- cbind(X_base, X_int)
    p    <- ncol(X)
    beta <- c(0.1, rep(0.4, n_base), rep(0.0, ncol(X_int)))
    y    <- rbinom(n, 1, plogis(cbind(1, X) %*% beta))
    fill   <- mean(X != 0)
    cat(sprintf("\n  Interaction dummies fill rate: %.4f  (n=%d, p=%d)", fill, n, p))
    list(X = X, y = y, family = binomial("logit"),
         label = sprintf("interaction dummies n_base=%d", n_base))
}

make_block_diagonal <- function(n, n_categories, p_per_cat, seed = 103) {
    set.seed(seed)
    
    # Each block contributes n rows, with non-zero values only in its p_per_cat columns
    # Total columns = n_categories * p_per_cat (true block-diagonal layout)
    blocks <- lapply(seq_len(n_categories), function(k) {
        X_k <- matrix(rnorm(n * p_per_cat), n, p_per_cat)
        
        # Zero-pad LEFT and RIGHT — don't reorder, just position the block correctly
        left_zeros  <- matrix(0, n, (k - 1) * p_per_cat)
        right_zeros <- matrix(0, n, (n_categories - k) * p_per_cat)
        cbind(left_zeros, X_k, right_zeros)
    })
    
    X   <- do.call(rbind, blocks)
    p   <- ncol(X)   # = n_categories * p_per_cat
    n_r <- nrow(X)   # = n * n_categories
    
    # beta needs length p+1 (intercept + p predictors)
    beta <- c(0.2, rep(c(0.4, -0.4), length.out = p))
    y    <- rbinom(n_r, 1, plogis(cbind(1, X) %*% beta))
    
    fill <- mean(X != 0)
    cat(sprintf("\n  Block diagonal fill rate: %.4f  (n=%d, p=%d)", fill, n_r, p))
    list(X = X, y = y, family = binomial("logit"),
         label = sprintf("block diag %d cats x p=%d", n_categories, p_per_cat))
}

sparse_scenarios <- list(
    # One-hot: fill = 1/(levels-1), gets sparse with many levels
    list(gen = make_onehot_logistic, n = 2000, n_factors = 10, levels_per_factor = 20),
    # Commented out due to required time
    #list(gen = make_onehot_logistic, n = 5000, n_factors = 20, levels_per_factor = 50), 
    list(gen = make_onehot_logistic, n = 10000, n_factors = 5,  levels_per_factor = 100),
    # Interaction dummies: fill = 0.15^2 = 0.0225 for interactions
    list(gen = make_interaction_logistic, n = 2000, n_base = 15),
    list(gen = make_interaction_logistic, n = 5000, n_base = 20),
    # Block diagonal
    list(gen = make_block_diagonal, n = 500,  n_categories = 10, p_per_cat = 8),
    list(gen = make_block_diagonal, n = 1000, n_categories = 20, p_per_cat = 10)
)

cat(sprintf("Benchmarking sparse with %d iterations per scenario\n", SPARSE_TIMES+1))

`%||%` <- function(x, y) if (!is.null(x)) x else y

results_sparse <- lapply(sparse_scenarios, function(sc) {
    sc_data <- sc$gen(
        sc$n,
        sc$n_factors %||% sc$n_base %||% sc$n_categories,
        sc$levels_per_factor %||% sc$p_per_cat %||% NULL
    )

    # Build a data.frame + formula, same as the original benchmark
    X   <- sc_data$X
    y   <- sc_data$y
    fam <- sc_data$family
    dat <- data.frame(y = y, X)
    fmla <- as.formula(paste("y ~", paste(names(dat)[-1], collapse = " + ")))

    cat(sprintf("\n %s \n", sc_data$label))

    mb <- microbenchmark(
        original = glm(fmla, data = dat, family = fam, type = TYPE,
                       method = "brglmFit_original", maxit = 10000),
        new      = glm(fmla, data = dat, family = fam, type = TYPE,
                       method = "brglmFit",          maxit = 10000),
        times = SPARSE_TIMES,
        unit  = "ms"
    )

    cat(sprintf("Coefficient fitting\n"))

    # Coefficient agreement check
    fit_orig <- glm(fmla, data = dat, family = fam, type = TYPE,
                    method = "brglmFit_original", maxit = 10000)
    fit_new  <- glm(fmla, data = dat, family = fam, type = TYPE,
                    method = "brglmFit",          maxit = 10000)

    cat(sprintf("  max |coef diff| original vs new: %.2e\n",
                max(abs(coef(fit_orig) - coef(fit_new)), na.rm = TRUE)))

    attr(mb, "label") <- sc_data$label
    attr(mb, "n")     <- nrow(X)
    attr(mb, "p")     <- ncol(X)
    mb
})


# Summary table

cat("\n\n== Timing summary (ms) ==\n")
cat(sprintf("%-40s  %-6s %-4s  %-9s  %12s  %12s\n",
            "scenario", "n", "p", "stat", "original", "new"))
cat(strrep("-", 100), "\n")

for (r in results) {
    label <- attr(r, "label")
    n     <- attr(r, "n")
    p     <- attr(r, "p")
    
    # r here is a raw microbenchmark data.frame: columns expr, time (nanoseconds)
    stats <- tapply(r$time / 1e6, r$expr, function(t)
        c(min = min(t), median = median(t), max = max(t)))
    
    for (stat in c("min", "median", "max")) {
        cat(sprintf("%-40s  %-6d %-4d  %-9s  %12.2f  %12.2f\n",
                    if (stat == "min") label else "",
                    n, p, stat,
                    stats[["original"]][[stat]],
                    stats[["new"]][[stat]]))
    }
    cat(strrep("-", 100), "\n")
}

cat("\n\n== Sparse timing summary (s) ==\n")
cat(sprintf("%-40s  %-6s %-4s  %-9s  %12s  %12s\n",
            "scenario", "n", "p", "stat", "original", "new"))
cat(strrep("-", 100), "\n")

for (r in results_sparse) {
    label <- attr(r, "label")
    n     <- attr(r, "n")
    p     <- attr(r, "p")
    
    # r here is a raw microbenchmark data.frame: columns expr, time (nanoseconds)
    stats <- tapply(r$time / 1e9, r$expr, function(t)
        c(min = min(t), median = median(t), max = max(t)))
    
    for (stat in c("min", "median", "max")) {
        cat(sprintf("%-40s  %-6d %-4d  %-9s  %12.2f  %12.2f\n",
                    if (stat == "min") label else "",
                    n, p, stat,
                    stats[["original"]][[stat]],
                    stats[["new"]][[stat]]))
    }
    cat(strrep("-", 100), "\n")
}