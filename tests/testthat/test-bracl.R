context("tests for bracl and its agreement with other methods")


## Stem cell research data from Agresti (2017, Analysis of ordinal data, Table 4.1)
freq <- c(34, 41, 58, 21, 30, 64,
          67, 83, 63, 52, 52, 50,
          30, 23, 15, 24, 18, 16,
          25, 14, 12, 15, 11, 11)
fund_research <- factor(rep(c("definitely", "probably", "probably not", "definitely not"), each = 6),
                        levels = c("definitely", "probably", "probably not", "definitely not"),
                        ordered = TRUE)
gender <- factor(rep(rep(c("female", "male"), each = 3), 4), levels = c("male", "female"))
## gender <- rep(rep(c(1, 0), each = 3), 4)
religion <- factor(rep(c("fundamentalist", "moderate", "liberal"), 8),
                   levels = c("fundamentalist", "moderate", "liberal"),
                   ordered = TRUE)
stemcell <- data.frame(frequency = freq, research = fund_research, gender = gender, religion = religion)

ff <- matrix(freq, nrow = 4, byrow = TRUE)
re <- c(1, 2, 3, 1, 2, 3)
ge <- c(1, 1, 1, 0, 0, 0)

library("VGAM")

## With proportional odds
fit_vgam_p <- vglm(cbind(ff[1, ], ff[2, ], ff[3, ], ff[4, ]) ~ re + ge, family = acat(reverse = TRUE, parallel = TRUE))
expect_warning(
    fit_bracl_p <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell, type = "ML", parallel = TRUE)
)

## Without proportional odds
fit_vgam <- vglm(cbind(ff[1, ], ff[2, ], ff[3, ], ff[4, ]) ~ re + ge, family = acat(reverse = TRUE, parallel = FALSE))
expect_warning(
    fit_bracl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell, type = "ML")
)

tol <- 1e-06
test_that("VGAM::vglm and bracl return the same coefficients", {
    expect_equal(unname(coef(fit_vgam)), unname(coef(fit_bracl)), tolerance = tol)
    expect_equal(unname(coef(fit_vgam_p)), unname(coef(fit_bracl_p)), tolerance = tol)
})

test_that("Difference in deviance is the same with VGAM::vglm and bracl", {
    expect_equal(logLik(fit_vgam) - logLik(fit_vgam_p),
                 unclass(logLik(fit_bracl) - logLik(fit_bracl_p)), tolerance = tol,
                 check.attributes = FALSE)
})

test_that("logLik returns the correct df for cracl", {
    expect_identical(attr(logLik(fit_bracl), "df"), as.integer(9))
    expect_identical(attr(logLik(fit_bracl_p), "df"), as.integer(5))
})

test_that("bracl returns the correct fitted values", {
    expect_equal(fitted(fit_bracl)[1:6, ], fit_vgam@fitted.values, tolerance = tol,
                 check.attributes = FALSE)
    expect_equal(fitted(fit_bracl_p)[1:6, ], fit_vgam_p@fitted.values, tolerance = tol,
                 check.attributes = FALSE)
})

shu <- function(dat)  dat[sample(seq.int(nrow(dat)), nrow(dat)), ]

test_that("bracl results is invariance to shuffling of the data", {
    for (j in 1:10) {

        expect_warning(
            fit_bracl_p_r <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = shu(stemcell), type = "ML", parallel = TRUE)
        )

        expect_warning(
            fit_bracl_r <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = shu(stemcell), type = "ML", parallel = FALSE)
        )

        expect_equal(coef(fit_bracl), coef(fit_bracl_r), tolerance = tol)
        expect_equal(coef(fit_bracl_p), coef(fit_bracl_p_r), tolerance = tol)
    }

})

tol  <-  1e-03
test_that("vcov method for bracl returns the correct vcov matrix", {
    expect_equal(vcov(fit_vgam_p), vcov(fit_bracl_p), tolerance = tol, check.attributes = FALSE)
    expect_equal(vcov(fit_vgam), vcov(fit_bracl), tolerance = tol, check.attributes = FALSE)
})


s1 <- summary(fit_bracl)
s1p <- summary(fit_bracl_p)
s2 <- summary(fit_vgam)
s2p <- summary(fit_vgam_p)

tol <- 1e-06
test_that("summary method for bracl returns the correct coef mat", {
    expect_equal(coef(s1), VGAM::coef(s2), tolerance = tol, check.attributes = FALSE)
    expect_equal(coef(s1p), VGAM::coef(s2p), tolerance = tol, check.attributes = FALSE)
})
