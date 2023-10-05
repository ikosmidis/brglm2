library("VGAM")

### tests for bracl and its agreement with other methods

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

## With proportional odds
suppressWarnings(fit_vgam_p <- vglm(cbind(ff[1, ], ff[2, ], ff[3, ], ff[4, ]) ~ re + ge, family = acat(reverse = TRUE, parallel = TRUE)))


expect_warning(
    fit_bracl_p <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell, type = "ML", parallel = TRUE)
)


fit_vgam <- vglm(cbind(ff[1, ], ff[2, ], ff[3, ], ff[4, ]) ~ re + ge, family = acat(reverse = TRUE, parallel = FALSE))
## Without proportional odds
expect_warning(
    fit_bracl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency, data = stemcell, type = "ML")
)

tol <- 1e-06
## "VGAM::vglm and bracl return the same coefficients
expect_equal(unname(coef(fit_vgam)), unname(coef(fit_bracl)), tolerance = tol)
expect_equal(unname(coef(fit_vgam_p)), unname(coef(fit_bracl_p)), tolerance = tol)

## Difference in deviance is the same with VGAM::vglm and bracl
expect_equal(logLik(fit_vgam) - logLik(fit_vgam_p),
             unclass(logLik(fit_bracl) - logLik(fit_bracl_p)), tolerance = tol,
             check.attributes = FALSE)


## logLik returns the correct df for cracl
expect_identical(attr(logLik(fit_bracl), "df"), as.integer(9))
expect_identical(attr(logLik(fit_bracl_p), "df"), as.integer(5))

## bracl returns the correct fitted values
expect_equal(fitted(fit_bracl)[1:6, ], fit_vgam@fitted.values, tolerance = tol,
             check.attributes = FALSE)
expect_equal(fitted(fit_bracl_p)[1:6, ], fit_vgam_p@fitted.values, tolerance = tol,
             check.attributes = FALSE)

shu <- function(dat)  dat[sample(seq.int(nrow(dat)), nrow(dat)), ]

## bracl results is invariance to shuffling of the data
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

tol  <-  1e-02
## vcov method for bracl returns the correct vcov matrix
expect_equal(vcov(fit_vgam_p), vcov(fit_bracl_p), tolerance = tol, check.attributes = FALSE)
expect_equal(vcov(fit_vgam), vcov(fit_bracl), tolerance = tol, check.attributes = FALSE)

s1 <- summary(fit_bracl)
s1p <- summary(fit_bracl_p)
s2 <- summary(fit_vgam)
s2p <- summary(fit_vgam_p)

tol <- 1e-06
## summary method for bracl returns the correct coef mat
expect_equal(coef(s1), s2@coef3, tolerance = tol, check.attributes = FALSE)
expect_equal(coef(s1p), s2p@coef3, tolerance = tol, check.attributes = FALSE)

newdata <- expand.grid(gender = c("male", "female"),  religion = c("moderate", "fundamentalist"))
## predict.bracl works as expected
pp <- predict(fit_bracl_p, newdata = stemcell, type = "probs")
p <- predict(fit_bracl, newdata = stemcell, type = "probs")
expect_equal(predict(fit_vgam_p, type = "response"),
             pp[19:24, ],
             tolerance = 1e-06, check.attributes = FALSE)
expect_equal(predict(fit_vgam, type = "response"),
             p[19:24, ],
             tolerance = 1e-06, check.attributes = FALSE)

## no intercept returns warning
expect_warning(
    fit_bracl_p_r <- bracl(research ~ -1 + as.numeric(religion) + gender, weights = frequency, data = shu(stemcell), type = "ML", parallel = FALSE)
)

## prediction with NAs works
newd <- newdata
newd[3, 2] <- NA
expect_true(is.na(predict(fit_bracl_p, newd, "class")[3]))
expect_true(all(is.na(predict(fit_bracl_p, newd, "probs")[3, ])))



## simulate from bracl objects returns a data frame with expected characteristics
simu_df <- simulate(fit_bracl_p)
nam_mf <- names(model.frame(fit_bracl_p))
nam_simu <- names(simu_df)
expect_identical(nrow(simu_df),
                 nrow(stemcell) * nlevels(stemcell$research))
expect_identical(levels(simu_df$research),
                 levels(stemcell$research))
expect_identical(is.ordered(simu_df$research),
                 is.ordered(stemcell$research))
expect_identical(nam_mf[!(nam_mf %in% nam_simu)],
                 "(weights)")
expect_identical(nam_simu[!(nam_simu %in% nam_mf)],
                 as.character(fit_bracl_p$call$weights))
expect_identical(nam_simu[(nam_simu %in% nam_mf)],
                 nam_mf[(nam_mf %in% nam_simu)])

## simulate from bracl objects returns a data frame with expected characteristics
simu_df <- simulate(fit_bracl)
nam_mf <- names(model.frame(fit_bracl))
nam_simu <- names(simu_df)
expect_identical(nrow(simu_df),
                 nrow(stemcell) * nlevels(stemcell$research))
expect_identical(levels(simu_df$research),
                 levels(stemcell$research))
expect_identical(is.ordered(simu_df$research),
                 is.ordered(stemcell$research))
expect_identical(nam_mf[!(nam_mf %in% nam_simu)],
                 "(weights)")
expect_identical(nam_simu[!(nam_simu %in% nam_mf)],
                 as.character(fit_bracl$call$weights))
expect_identical(nam_simu[(nam_simu %in% nam_mf)],
                 nam_mf[(nam_mf %in% nam_simu)])
