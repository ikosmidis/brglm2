b_control <- brglmControl()

## the object brglmControl() returns is as expected
expect_identical(b_control$epsilon, 1e-06)
expect_identical(b_control$maxit, 100)
expect_false(b_control$trace)
expect_null(b_control$response_adjustment)
expect_equal(b_control$Trans, expression(dispersion))
expect_equal(b_control$inverseTrans, expression(transformed_dispersion))
expect_equal(b_control$transformation, "identity")
expect_equal(b_control$slowit, 1)
expect_equal(b_control$max_step_factor, 12)

b_control <- brglmControl(epsilon = 1e-02, trace = TRUE, response_adjustment = c(0.3, 0.2))
## the object brglmControl returns with defaults is as expected
expect_identical(b_control$epsilon, 1e-02)
expect_identical(b_control$maxit, 100)
expect_true(b_control$trace)
expect_identical(b_control$response_adjustment, c(0.3, 0.2))
expect_equal(b_control$Trans, expression(dispersion))
expect_equal(b_control$inverseTrans, expression(transformed_dispersion))
expect_equal(b_control$transformation, "identity")
expect_equal(b_control$slowit, 1)
expect_equal(b_control$max_step_factor, 12)

data("coalition", package = "brglm2")
mod <- glm(duration ~ fract + numst2, family = Gamma, data = coalition, method = "brglmFit")
expect_warning(mod <- update(mod, maxit = 7, transformation = "sqrt"), "brglmFit: algorithm did not converge")
expect_equal(mod$iter, 7)
expect_equal(mod$transformation, "sqrt")

b_control <- brglmControl(epsilon = 1e-02, ABCDEFG123 = 1, response_adjustment = c(0.3, 0.2), trace = TRUE, )
## the object brglmControl returns with defaults is as expected
expect_identical(b_control$epsilon, 1e-02)
expect_identical(b_control$maxit, 100)
expect_true(b_control$trace)
expect_identical(b_control$response_adjustment, c(0.3, 0.2))
expect_equal(b_control$Trans, expression(dispersion))
expect_equal(b_control$inverseTrans, expression(transformed_dispersion))
expect_equal(b_control$transformation, "identity")
expect_equal(b_control$slowit, 1)
expect_equal(b_control$max_step_factor, 12)

## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))

clot_formula <- conc ~ lot*log(u) + I(2*log(u))

## check_aliasing option works as expected
expect_error(glm(clot_formula, data = clotting, family = Gamma(), method = brglm_fit, check_aliasing = FALSE),
             pattern = "NA/NaN/Inf in foreign")
mod <- glm(clot_formula, data = clotting, family = Gamma(), method = brglm_fit, check_aliasing = TRUE)
expect_true(is.na(coef(mod)["I(2 * log(u))"]))
## check defaults
mod <- glm(clot_formula, data = clotting, family = Gamma(), method = brglm_fit)
expect_true(is.na(coef(mod)["I(2 * log(u))"]))

## data("sex2", package = "logistf")
## ## brglmControl arguments can be passed directly from the brglmFit call
## expect_warning(glm(case ~ dia, data = sex2, family = binomial(), method = brglm_fit))
## ## Set response_adjustment to 2 just to avoid non-integer successes warning
## m1 <- glm(case ~ dia, data = sex2, family = binomial(), method = brglm_fit, response_adjustment = 2)
## expect_equal(m1$iter, 7, tolerance = .Machine$double.eps/2)
