data("stemcell", package = "brglm2")

# Adjacent category logit (non-proportional odds)
fit_bracl <- bracl(research ~ religion + gender, weights = frequency,
                   data = stemcell, type = "ML")
# Adjacent category logit (proportional odds)
fit_bracl_p <- bracl(research ~ religion + gender, weights = frequency,
                    data = stemcell, type = "ML", parallel = TRUE)

tol <- 1e-08
## ordinal_superiority agrees with fitted probabilities for main-effects models
gamma_fun <- function(p0, p1) {
    out <- outer(p0, p1, "*")
    sum(out[upper.tri(out)]) + sum(diag(out)) / 2
}
expected_ordsup <- function(fit) {
    nd0 <- data.frame(religion = factor(levels(stemcell$religion),
                                        levels = levels(stemcell$religion),
                                        ordered = is.ordered(stemcell$religion)),
                      gender = factor("male", levels = levels(stemcell$gender)))
    nd1 <- nd0
    nd1$gender <- factor("female", levels = levels(stemcell$gender))
    p0 <- predict(fit, nd0, type = "probs")
    p1 <- predict(fit, nd1, type = "probs")
    vapply(seq_len(nrow(p0)), function(i) gamma_fun(p0[i, ], p1[i, ]), numeric(1))
}
os_np <- ordinal_superiority(fit_bracl, ~ gender)
os_p <- ordinal_superiority(fit_bracl_p, ~ gender)

expect_equal(unname(os_np[, "gamma"]), expected_ordsup(fit_bracl), tolerance = tol)
expect_equal(unname(os_p[, "gamma"]), expected_ordsup(fit_bracl_p), tolerance = tol)
expect_equal(unname(ordinal_superiority(fit_bracl, ~ gender, measure = "Delta")[, "Delta"]),
             2 * unname(os_np[, "gamma"]) - 1, tolerance = tol)
expect_equal(unname(ordinal_superiority(fit_bracl_p, ~ gender, measure = "Delta")[, "Delta"]),
             2 * unname(os_p[, "gamma"]) - 1, tolerance = tol)

## data argument works with predictor-only data
pred_data <- stemcell[c("religion", "gender")]
os_np_data <- ordinal_superiority(fit_bracl, ~ gender, data = pred_data)
os_p_data <- ordinal_superiority(fit_bracl_p, ~ gender, data = pred_data)
expect_equal(os_np_data, os_np, tolerance = tol, check.attributes = FALSE)
expect_equal(os_p_data, os_p, tolerance = tol, check.attributes = FALSE)

# Adjacent category logit models with interaction
fit_bracl_int <- bracl(research ~ religion * gender, weights = frequency,
                       data = stemcell, type = "ML")
fit_bracl_p_int <- bracl(research ~ religion * gender, weights = frequency,
                         data = stemcell, type = "ML", parallel = TRUE)

expected_ordsup_int <- function(fit) {
    nd0 <- data.frame(religion = factor(levels(stemcell$religion),
                                        levels = levels(stemcell$religion),
                                        ordered = is.ordered(stemcell$religion)),
                      gender = factor("male", levels = levels(stemcell$gender)))
    nd1 <- nd0
    nd1$gender <- factor("female", levels = levels(stemcell$gender))
    p0 <- predict(fit, nd0, type = "probs")
    p1 <- predict(fit, nd1, type = "probs")
    vapply(seq_len(nrow(p0)), function(i) gamma_fun(p0[i, ], p1[i, ]), numeric(1))
}

os_np_int <- ordinal_superiority(fit_bracl_int, ~ gender)
os_p_int <- ordinal_superiority(fit_bracl_p_int, ~ gender)

expect_identical(nrow(os_np_int), 3L)
expect_identical(nrow(os_p_int), 3L)
expect_equal(unname(os_np_int[, "gamma"]), expected_ordsup_int(fit_bracl_int), tolerance = tol)
expect_equal(unname(os_p_int[, "gamma"]), expected_ordsup_int(fit_bracl_p_int), tolerance = tol)

pred_data_int <- stemcell[c("religion", "gender")]
os_np_int_data <- ordinal_superiority(fit_bracl_int, ~ gender, data = pred_data_int)
os_p_int_data <- ordinal_superiority(fit_bracl_p_int, ~ gender, data = pred_data_int)
expect_equal(os_np_int_data, os_np_int, tolerance = tol, check.attributes = FALSE)
expect_equal(os_p_int_data, os_p_int, tolerance = tol, check.attributes = FALSE)

# Adjacent category logit models with transformed covariates in formula
fit_bracl_tf <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
                      data = stemcell, type = "ML")
fit_bracl_p_tf <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
                        data = stemcell, type = "ML", parallel = TRUE)

os_np_tf <- ordinal_superiority(fit_bracl_tf, ~ gender)
os_p_tf <- ordinal_superiority(fit_bracl_p_tf, ~ gender)

expect_equal(unname(os_np_tf[, "gamma"]), expected_ordsup(fit_bracl_tf), tolerance = tol)
expect_equal(unname(os_p_tf[, "gamma"]), expected_ordsup(fit_bracl_p_tf), tolerance = tol)

pred_data_tf <- stemcell[c("religion", "gender")]
os_np_tf_data <- ordinal_superiority(fit_bracl_tf, ~ gender, data = pred_data_tf)
os_p_tf_data <- ordinal_superiority(fit_bracl_p_tf, ~ gender, data = pred_data_tf)
expect_equal(os_np_tf_data, os_np_tf, tolerance = tol, check.attributes = FALSE)
expect_equal(os_p_tf_data, os_p_tf, tolerance = tol, check.attributes = FALSE)

# Something odd
fit_bracl_tf <- bracl(research ~ sqrt(as.numeric(religion)) + gender, weights = frequency,
                      data = stemcell, type = "ML")
fit_bracl_p_tf <- bracl(research ~ sqrt(as.numeric(religion)) + gender, weights = frequency,
                        data = stemcell, type = "ML", parallel = TRUE)

os_np_tf <- ordinal_superiority(fit_bracl_tf, ~ gender)
os_p_tf <- ordinal_superiority(fit_bracl_p_tf, ~ gender)

expect_equal(unname(os_np_tf[, "gamma"]), expected_ordsup(fit_bracl_tf), tolerance = tol)
expect_equal(unname(os_p_tf[, "gamma"]), expected_ordsup(fit_bracl_p_tf), tolerance = tol)

pred_data_tf <- stemcell[c("religion", "gender")]
os_np_tf_data <- ordinal_superiority(fit_bracl_tf, ~ gender, data = pred_data_tf)
os_p_tf_data <- ordinal_superiority(fit_bracl_p_tf, ~ gender, data = pred_data_tf)
expect_equal(os_np_tf_data, os_np_tf, tolerance = tol, check.attributes = FALSE)
expect_equal(os_p_tf_data, os_p_tf, tolerance = tol, check.attributes = FALSE)
