% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brmultinom.R, R/zzz_conventions.R
\name{residuals.brmultinom}
\alias{residuals.brmultinom}
\alias{residuals.bracl}
\title{Residuals for multinomial logistic regression and adjacent category logit models}
\usage{
\method{residuals}{brmultinom}(object, type = c("pearson", "response", "deviance", "working"), ...)

\method{residuals}{bracl}(object, type = c("pearson", "response", "deviance", "working"), ...)
}
\arguments{
\item{object}{the object coming out of \code{\link[=bracl]{bracl()}} and
\code{\link[=brmultinom]{brmultinom()}}.}

\item{type}{the type of residuals which should be returned.  The
options are: \code{"pearson"} (default), \code{"response"}, \code{"deviance"},
\code{"working"}. See Details.}

\item{...}{Currently not used.}
}
\description{
Residuals for multinomial logistic regression and adjacent category logit models
}
\details{
The residuals computed are the residuals from the equivalent
Poisson log-linear model fit, organized in a form that matches the
output of \code{fitted(object, type = "probs")}. As a result, the output
is residuals defined in terms of the object and expected
multinomial counts.
}
\seealso{
brmultinom bracl
}
