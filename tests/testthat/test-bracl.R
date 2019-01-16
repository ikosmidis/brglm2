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
fit_adj_p <- vglm(cbind(ff[1, ], ff[2, ], ff[3, ], ff[4, ]) ~ re + ge, family = acat(reverse = TRUE, parallel = TRUE))
expect_warning(
    fit_acl_p <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell, type = "ML", parallel = TRUE)
)

## Without proportional odds
fit_adj <- vglm(cbind(ff[1, ], ff[2, ], ff[3, ], ff[4, ]) ~ re + ge, family = acat(reverse = TRUE, parallel = FALSE))
expect_warning(
    fit_acl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell, type = "ML")
)

tol <- 1e-06
test_that("VGAM::vglm and bracl return the same coefficients", {
    expect_equal(unname(coef(fit_adj)), unname(coef(fit_acl)), tolerance = tol)
    expect_equal(unname(coef(fit_adj_p)), unname(coef(fit_acl_p)), tolerance = tol)
})

test_that("Difference in deviance is the same with VGAM::vglm and bracl", {
    expect_equal(logLik(fit_adj) - logLik(fit_adj_p),
                 unclass(logLik(fit_acl) - logLik(fit_acl_p)), tolerance = tol,
                 check.attributes = FALSE)
})

test_that("logLik returns the correct df for cracl", {
    expect_identical(attr(logLik(fit_acl), "df"), as.integer(9))
    expect_identical(attr(logLik(fit_acl_p), "df"), as.integer(5))
})

test_that("bracl returns the correct fitted values", {
    expect_equal(fitted(fit_acl)[1:6, ], fit_adj@fitted.values, tolerance = tol,
                 check.attributes = FALSE)
    expect_equal(fitted(fit_acl_p)[1:6, ], fit_adj_p@fitted.values, tolerance = tol,
                 check.attributes = FALSE)
})

shu <- function(dat)  dat[sample(seq.int(nrow(dat)), nrow(dat)), ]

test_that("bracl results is invariance to shuffling of the data", {
    for (j in 1:10) {

        expect_warning(
            fit_acl_p_r <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = shu(stemcell), type = "ML", parallel = TRUE)
        )

        expect_warning(
            fit_acl_r <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = shu(stemcell), type = "ML", parallel = FALSE)
        )

        expect_equal(coef(fit_acl), coef(fit_acl_r), tolerance = tol)
        expect_equal(coef(fit_acl_p), coef(fit_acl_p_r), tolerance = tol)
    }

})

tol  <-  1e-03
test_that("vcov method for bracl returns the correct vcov matrix", {
    expect_equal(vcov(fit_adj_p), vcov(fit_acl_p), tolerance = tol, check.attributes = FALSE)
    inds <- c(1, 4, 7, 2, 5, 8, 3,6, 9)
    expect_equal(vcov(fit_adj), vcov(fit_acl)[inds, inds], tolerance = tol, check.attributes = FALSE)
})
