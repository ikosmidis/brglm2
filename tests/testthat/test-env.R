
library("lobstr")
## Here xx is not copied

xx <- matrix(rnorm(1000), 500, 2)
out <- set_brglmFit_env(
    x = xx,
    y = rbinom(500, 1, 0.3),
    weights = 1,
    offset = NULL,
    family = binomial(),
    fixed_totals = NULL,
    control = brglm_control(),
    singular.ok = FALSE,
    start = c(1,2,3),
    env = TRUE)
expect_identical(ref(xx), ref(out$x))


xx <- matrix(rnorm(1000), 500, 2)
out2 <- set_brglmFit_env(
    x = xx,
    y = rbinom(500, 1, 0.3),
    weights = 1,
    offset = NULL,
    family = binomial(),
    fixed_totals = NULL,
    control = brglm_control(),
    singular.ok = FALSE,
    start = c(1,2,3),
    env = FALSE)
expect_identical(ref(xx), ref(out$x))

## Here x will be copied
xx <- as.data.frame(matrix(rnorm(1000), 500, 2))
out <- create_brglmFit_env(
    x = xx,
    y = rbinom(500, 1, 0.3),
    weights = 1,
    offset = NULL,
    family = binomial(),
    fixed_totals = NULL,
    control = brglm_control(),
    singular.ok = FALSE)
expect_false(identical(ref(xx), ref(out$x)))
