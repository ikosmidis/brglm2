context("implementation of derivatives of mean with respect to eta (comparison with numerical derivatives)")

library("numDeriv")


links <- c("logit", "probit", "cauchit", "cloglog",
           "identity", "log", "sqrt",
           "1/mu^2", "inverse")

tol <- 0.00001 # sqrt(.Machine$double.eps)

for (link in links) {
    etas <- if (link == "1/mu^2") seq(0.1, 8, length = 100) else seq(-8, 8, length = 100)

    clink <- enrichLink(make.link(link))

    test_that(paste("mu.eta is correctly implemented for", link), {
        expect_equal(grad(clink$linkinv, etas), clink$mu.eta(etas), tolerance = tol)
    })

    test_that(paste("d2mu.deta is correctly implemented for", link), {
        expect_equal(grad(clink$mu.eta, etas), clink$d2mu.deta(etas), tolerance = tol)
    })

    test_that(paste("d3mu.deta is correctly implemented for", link), {
        expect_equal(grad(clink$d2mu.deta, etas), clink$d3mu.deta(etas), tolerance = tol)
    })

}

## Power family
etas <- seq(0.1, 8, length = 100)
tol <- 0.0001
for (lambda in c(0, 1, 0.2, 0.4, 0.553, 0.98) ) {
    clink <- enrichLink(power(lambda = lambda))

    test_that(paste("mu.eta is correctly implemented for power with lambda", lambda), {
        expect_equal(grad(clink$linkinv, etas), clink$mu.eta(etas), tolerance = tol)
    })

    test_that(paste("d2mu.deta is correctly implemented for power with lambda", lambda), {
        expect_equal(grad(clink$mu.eta, etas), clink$d2mu.deta(etas), tolerance = tol)
    })

    test_that(paste("d3mu.deta is correctly implemented for power with lambda", lambda), {
        expect_equal(grad(clink$d2mu.deta, etas), clink$d3mu.deta(etas), tolerance = tol)
    })
}
