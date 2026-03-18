## benchmark_largep.R
## Compares brglmFit_original, brglmFit (base/CG default), brglmFit_dogleg,
## and brglmFit_cg across multiple synthetic data scenarios covering different
## families, link functions, sparsity structures, and high-dim settings.

library(brglm2)
library(microbenchmark)

# ══════════════════════════════════════════════════════════════════════════════
# Data generators
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. Logistic regression, AR-1 predictors (baseline) ───────────────────────
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

# ── 2. Probit regression, block-correlated predictors ────────────────────────
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

# ── 3. Poisson log-linear, sparse predictors ─────────────────────────────────
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

# ── 4. Logistic regression, high-dimensional medical (n ~ 5p) ────────────────
# Toeplitz covariance, sparse signal — mimics omics/GWAS logistic models
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

# ── 5. Complementary log-log, survival-style binary outcome ──────────────────
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

# ══════════════════════════════════════════════════════════════════════════════
# Scenarios: (generator, n, p)
# ══════════════════════════════════════════════════════════════════════════════
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

TIMES <- 20L

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════
results <- lapply(scenarios, function(sc) {
    sc_data <- sc$gen(sc$n, sc$p)
    n       <- sc$n;  p <- sc$p
    label   <- sc_data$label
    fam     <- sc_data$family
    dat     <- sc_data$data
    fmla    <- as.formula(paste("y ~", paste(names(dat)[-1], collapse = " + ")))

    cat(sprintf("\n── %s  n=%d  p=%d ──\n", label, n, p))

    mb <- microbenchmark(
        original = glm(fmla, data = dat, family = fam,
                       method = "brglmFit_original", maxit = 10000),
        new     = glm(fmla, data = dat, family = fam,
                       method = "brglmFit",          maxit = 10000),
        times = TIMES,
        unit  = "ms"
    )

    # Coefficient agreement
    fit_orig   <- glm(fmla, data = dat, family = fam,
                      method = "brglmFit_original", maxit = 10000)
    fit_new   <- glm(fmla, data = dat, family = fam,
                      method = "brglmFit", maxit = 10000)

    cat(sprintf("  max |coef diff| original vs new :   %.2e\n",
                max(abs(coef(fit_orig) - coef(fit_new)),   na.rm = TRUE)))

    smry         <- summary(mb)
    smry$n       <- n
    smry$p       <- p
    smry$label   <- label
    smry
})

# ══════════════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════════════
cat("\n\n== Timing summary (ms) ==\n")
cat(sprintf("%-28s  %-6s %-4s  %-9s  %12s  %12s \n",
            "scenario", "n", "p", "stat", "original", "new"))
cat(strrep("-", 95), "\n")

for (r in results) {
    stats <- setNames(
        lapply(split(r, r$expr), function(d)
            c(min = d$min, median = d$median, max = d$max)),
        r$expr[!duplicated(r$expr)]
    )
    for (stat in c("min", "median", "max")) {
        cat(sprintf("%-28s  %-6d %-4d  %-9s  %12.2f  %12.2f\n",
                    if (stat == "min") r$label[1] else "",
                    r$n[1], r$p[1], stat,
                    stats[["original"]][[stat]],
                    stats[["new"]][[stat]]))
    }
    cat(strrep("-", 95), "\n")
}