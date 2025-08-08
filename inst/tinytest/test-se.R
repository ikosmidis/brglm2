library("brglm2")

kappa0 <- 0.2
gamma0 <- 5
alpha0 <- 0.88
theta0 <- 1

iota0 <- 2
mu0 <- 0.7
b0 <- 1.2
sigma0 <- 2.3

expect_equal(se0(mu0, b0, sigma0, kappa0, gamma0, alpha0),
             se1(mu0, b0, sigma0, 0, kappa0, gamma0, alpha0, 0)[1:3])

## 3 par
soln0 <- solve_se(kappa0, gamma0, alpha0, start = c(0.5, 1, 1),
                  init_iter = 23, init_method = "BFGS")
expect_equal(attr(soln0, "funcs"), rep(0, 3))

opt_str <- c("optim", "23", "BFGS", "nleqslv")
for (st in opt_str) {
    expect_true(grepl(st, attr(soln0, "optimization-chain")))
}

nu0 <- sqrt(soln0[1]^2 * gamma0^2 + kappa0 * soln0[3]^2)
soln0_c <- solve_se(kappa0, nu0, alpha0, start = soln0, corrupted = TRUE,
                    init_iter = "only", init_method = "Nelder-Mead")
expect_equal(soln0_c, soln0, check.attributes = FALSE, tolerance = 1e-07)

## 4 par
soln1 <- solve_se(kappa0, gamma0, alpha0, theta0, start = c(0.5, 1, 1, 1),
                  init_iter = 44, init_method = "Nelder-Mead")
expect_equal(attr(soln1, "funcs"), rep(0, 4))

soln1t <- solve_se(kappa0, gamma0, alpha0, theta0, start = c(0.5, 1, 1, 1),
                   init_iter = 44, init_method = "Nelder-Mead", transform = FALSE)
expect_equal(attr(soln1t, "funcs"), rep(0, 4))
expect_equal(soln1t, soln1, tolerance = 1e-07)

opt_str <- c("optim", "44", "Nelder-Mead", "nleqslv")
for (st in opt_str) {
    expect_true(grepl(st, attr(soln1, "optimization-chain")))
}

nu0 <- sqrt(soln1[1]^2 * gamma0^2 + kappa0 * soln1[3]^2)
soln1_c <- solve_se(kappa0, nu0, alpha0, soln1[4], start = soln1, corrupted = TRUE,
                    init_iter = 44, init_method = "Nelder-Mead")
expect_equal(soln1_c[1:3], soln1[1:3], check.attributes = FALSE, tolerance = 1e-07)
expect_equal(soln1_c[4], theta0, check.attributes = FALSE, tolerance = 1e-07)


## Reproducing Table 13 of Zhao et al. (2022)
## at https://doi.org/10.3150/21-BEJ1401

thetas <- c(0, 0.5, 1, 2, 2.5)
gamma0 <- 5
pars3 <- matrix(NA, length(thetas), 3)
pars4 <- matrix(NA, length(thetas), 4)
colnames(pars4) <- c("I_mu", "I_b", "I_sigma", "I_iota")
colnames(pars3) <- c("I_mu", "I_b", "I_sigma")
for (i in seq_along(thetas)) {
    pars3[i, ] <- solve_se(kappa = 0.2, ss = sqrt(5 + thetas[i]^2), alpha = 1, start = c(0.5, 1, 1), init_iter = 0)
    pars4[i, ] <- solve_se(kappa = 0.2, ss = sqrt(5), alpha = 1, intercept = thetas[i], start = c(pars3[i, ], thetas[i]), init_iter = 0)
}

cbind(pars4, pars3)
Zhaoetal_t13_I <- matrix(c(1.5, 1.51, 1.56, 1.83, 2.31,
                           3.03, 3.13, 3.45, 5.47, 8.96,
                           4.75, 4.84, 5.16, 7.01, 10.00,
                           0, 0.76, 1.559, 3.68, 5.8), ncol = 4)
Zhaoetal_t13_II <- matrix(c(1.5, 1.51, 1.55, 1.75, 1.95,
                           3.03, 3.12, 3.42, 4.83, 6.26,
                           4.75, 4.84, 5.13, 6.45, 7.73), ncol = 3)
expect_equal(Zhaoetal_t13_II, pars3, tolerance = 1e-02, check.attributes = FALSE)
expect_equal(Zhaoetal_t13_I, pars4, tolerance = 1e-01, check.attributes = FALSE)
