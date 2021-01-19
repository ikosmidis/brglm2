context("brglmControl")

b_control <- brglmControl()

test_that("the object brglmControl() returns is as expected", {
    expect_identical(b_control$epsilon, 1e-06)
    expect_identical(b_control$maxit, 100)
    expect_false(b_control$trace)
    expect_null(b_control$response_adjustment)
    expect_equal(b_control$Trans, expression(dispersion))
    expect_equal(b_control$inverseTrans, expression(transformed_dispersion))
    expect_equal(b_control$transformation, "identity")
    expect_equal(b_control$slowit, 1)
    expect_equal(b_control$max_step_factor, 12)
})

b_control <- brglmControl(epsilon = 1e-02, trace = TRUE, response_adjustment = c(0.3, 0.2))
test_that("the object brglmControl returns with defaults is as expected", {
    expect_identical(b_control$epsilon, 1e-02)
    expect_identical(b_control$maxit, 100)
    expect_true(b_control$trace)
    expect_identical(b_control$response_adjustment, c(0.3, 0.2))
    expect_equal(b_control$Trans, expression(dispersion))
    expect_equal(b_control$inverseTrans, expression(transformed_dispersion))
    expect_equal(b_control$transformation, "identity")
    expect_equal(b_control$slowit, 1)
    expect_equal(b_control$max_step_factor, 12)
})

data("sex2", package = "logistf")

test_that("brglmControl arguments can be passed directly from the brglmFit call", {
    expect_warning(glm(case ~ dia, data = sex2, family = binomial(), method = brglm_fit))
    ## Set response_adjustment to 2 just to avoid non-integer successes warning
    m1 <- glm(case ~ dia, data = sex2, family = binomial(), method = brglm_fit, response_adjustment = 2)
    expect_equal(m1$iter, 7, tolerance = .Machine$double.eps/2)
})



b_control <- brglmControl(epsilon = 1e-02, ABCDEFG123 = 1, response_adjustment = c(0.3, 0.2), trace = TRUE, )
test_that("the object brglmControl returns with defaults is as expected", {
    expect_identical(b_control$epsilon, 1e-02)
    expect_identical(b_control$maxit, 100)
    expect_true(b_control$trace)
    expect_identical(b_control$response_adjustment, c(0.3, 0.2))
    expect_equal(b_control$Trans, expression(dispersion))
    expect_equal(b_control$inverseTrans, expression(transformed_dispersion))
    expect_equal(b_control$transformation, "identity")
    expect_equal(b_control$slowit, 1)
    expect_equal(b_control$max_step_factor, 12)
})
