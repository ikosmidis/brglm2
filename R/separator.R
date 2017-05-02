#  Port of the separator function from the safeBinaryRegression
#  (version 0.1-3) R package (see safeBinaryRegression/R/separator.q)
#
#  Copyright (C) 2017 Ioannis Kosmidis
#  Copyright (C) 2009-2013 Kjell Konis
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
separator <- function(x, y, linear_program = c("primal", "dual"), purpose = c("test", "find"),
                      beta_tolerance = 1e-03) {
    n <- dim(x)[1L]
    p <- dim(x)[2L]
    p_seq <- seq.int(p)
    zeros <- rep.int(0, n)
    dimnames(x) <- NULL ## FIXME: do we really need that?
    y.bar <- -sign(y - 0.5)
    x.bar <- y.bar * x
    ans <- list()
    linear_program <- match.arg(linear_program)
    purpose <- match.arg(purpose)
    if (linear_program == "primal" && purpose == "test") {
        lp <- lpSolveAPI::make.lp(n, p)
        for(j in p_seq) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[, j])
        }
        status <- lpSolveAPI::set.rhs(lp,  zeros)
        status <- lpSolveAPI::set.constr.type(lp, rep.int(1, n))
        status <- lpSolveAPI::set.objfn(lp, -colSums(x.bar))
        status <- lpSolveAPI::set.bounds(lp, lower = rep(-Inf, p), upper = rep(Inf, p))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "max",
                                          simplextype = c("primal", "primal"))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        ## end work
        if (status == 0) {
            ans$separation <- FALSE
        }
        else {
            if (status == 3) {
                ans$separation <- TRUE
            }
            else {
                stop("unexpected result from lpSolveAPI for primal test")
            }
        }
    }
    if (linear_program == "primal" && purpose == "find") {
        lp <- lpSolveAPI::make.lp(n, p)
        for (j in p_seq) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[, j])
        }
        status <- lpSolveAPI::set.rhs(lp, zeros)
        status <- lpSolveAPI::set.constr.type(lp, rep.int(1, n))
        status <- lpSolveAPI::set.objfn(lp, -colSums(x.bar))
        status <- lpSolveAPI::set.bounds(lp, lower = rep.int(-1, p), upper = rep.int(1, p))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "max",
                              simplextype = c("primal", "primal"))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        if (status != 0) {
            stop("unexpected result from lpSolveAPI for primal test")
        }
        beta <- lpSolveAPI::get.variables(lp)
        if (any(abs(beta) > beta_tolerance)) {
            ans$separation <- TRUE
        }
        else {
            ans$separation <- FALSE
        }
        ans$beta <- beta
    }
    if (linear_program == "dual" && purpose == "test") {
        lp <- lpSolveAPI::make.lp(p, n)
        for (j in 1:n) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[j, ])
        }
        status <- lpSolveAPI::set.rhs(lp, -colSums(x.bar))
        status <- lpSolveAPI::set.constr.type(lp, rep.int(3, p))
        status <- lpSolveAPI::set.objfn(lp, zeros)
        status <- lpSolveAPI::set.bounds(lp, lower = zeros, upper = rep(Inf, n))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "min",
                                          simplextype = c("primal", "primal"))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        if (status == 0) {
            ans$separation <- FALSE
        }
        else {
            if (status == 2) {
                ans$separation <- TRUE
            }
            else {
                stop("unexpected result from lpSolveAPI for dual test")}
        }
    }
    if (linear_program == "dual" && purpose == "find") {
        lp <- lpSolveAPI::make.lp(p, n + 2*p)
        for (j in 1:n) {
            status <- lpSolveAPI::set.column(lp, j, x.bar[j, ])
        }
        for (j in p_seq) {
            status <- lpSolveAPI::set.column(lp, n+j, -1.0, j)
        }
        ## IK, 12 April 2017: p_seq instead 1:n below;
        ## safeBinaryRegression:::separator (version 0.1-3) has 1:n
        for (j in p_seq) {
            status <- lpSolveAPI::set.column(lp, n+p+j, 1.0, j)
        }
        b <- -colSums(x.bar)
        status <- lpSolveAPI::set.rhs(lp, b)
        status <- lpSolveAPI::set.constr.type(lp, rep.int(3, p))
        status <- lpSolveAPI::set.objfn(lp, rep.int(c(0.0, 1.0), c(n, 2*p)))
        status <- lpSolveAPI::set.bounds(lp, lower = rep.int(0, n + 2*p), upper = rep(Inf, n + 2*p))
        control <- lpSolveAPI::lp.control(lp, pivoting = "firstindex", sense = "min",
                                          simplextype = c("primal", "primal"))
        basis <- p_seq
        basis[b >= 0.0] <- basis[b >= 0.0] + p
        status <- lpSolveAPI::set.basis(lp, -(n + p + basis))
        status <- lpSolveAPI::solve.lpExtPtr(lp)
        beta <- lpSolveAPI::get.dual.solution(lp)[2:(p+1)]
        if (all(abs(beta) > beta_tolerance)) {
            ans$separation <- TRUE
        }
        else {
            ans$separation <- FALSE
        }
        ans$beta <- beta
    }
    ans
}
