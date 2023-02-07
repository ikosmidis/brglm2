# Copyright (C) 2018- Ioannis Kosmidis

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

#' A [`"link-glm"`][make.link] object for misclassified responses in binomial regression models
#'
#' [mis()] is a [`"link-glm"`][make.link] object that specifies the link function in Neuhaus (1999, expression (8)) for handling misclassified responses in binomial regression models using maximum likelihood. A prior specification of the sensitivity and specificity is required.
#'
#' @param link the baseline link to be used.
#' @param sensitivity the probability of observing a success given that a success actually took place given any covariate values.
#' @param specificity the probability of observing a failure given that a failure actually took place given any covariate values.
#'
#' @details
#'
#' `sensitivity + specificity` should be greater or equal to 1,
#' otherwise it is implied that the procedure producing the responses
#' performs worse than chance in terms of misclassification.
#'
#' @references
#'
#' Neuhaus J M (1999). Bias and efficiency loss due to misclassified
#' responses in binary regression. Biometrika, **86**, 843-855.
#' \url{https://www.jstor.org/stable/2673589}.
#'
#' @seealso [glm()], [brglm_fit()]
#'
#' @examples
#'
#' ## Define a few links with some misclassification
#' logit_mis <- mis(link = "logit", sensitivity = 0.9, specificity = 0.9)
#'
#' lizards_f <- cbind(grahami, opalinus) ~ height + diameter + light + time
#'
#' lizardsML <- glm(lizards_f, family = binomial(logit), data = lizards)
#'
#' lizardsML_mis <- update(lizardsML, family = binomial(logit_mis),
#'                         start = coef(lizardsML))
#'
#' ## A notable change is coefficients is noted here compared to when
#' ## specificity and sensitity are 1
#' coef(lizardsML)
#' coef(lizardsML_mis)
#'
#' ## Bias reduction is also possible
#' update(lizardsML_mis, method = "brglmFit", type = "AS_mean",
#'        start = coef(lizardsML))
#'
#' update(lizardsML_mis, method = "brglmFit", type = "AS_median",
#'        start = coef(lizardsML))
#'
#' @export
mis <- function(link = "logit", sensitivity = 1, specificity = 1) {
    link <- make.link(link)
    linkfun <- function(mu) {
        link$linkfun((mu -1 + specificity) / (sensitivity + specificity - 1))
    }
    linkinv <- function(eta) {
        (sensitivity + specificity - 1) * link$linkinv(eta) + 1 - specificity
    }
    mu.eta <- function(eta) {
        (sensitivity + specificity - 1) * link$mu.eta(eta)
    }
    valideta <- function(eta) {
        TRUE
    }
    structure(list(linkfun = linkfun,
                   linkinv = linkinv,
                   mu.eta = mu.eta,
                   valideta = valideta,
                   name = "miss"),
              class = "link-glm")
}

