% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bracl.R
\name{predict.bracl}
\alias{predict.bracl}
\title{Predict method for \link{bracl} fits}
\usage{
\method{predict}{bracl}(object, newdata, type = c("class", "probs"), ...)
}
\arguments{
\item{object}{a fitted object of class inheriting from \code{\link[=bracl]{"bracl"}}.}

\item{newdata}{optionally, a data frame in which to look for
variables with which to predict.  If omitted, the fitted linear
predictors are used.}

\item{type}{the type of prediction required. The default is
\code{"class"}, which produces predictions of the response category
at the covariate values supplied in \code{"newdata"}, selecting the
category with the largest probability; the alternative
\code{"probs"} returns all category probabilities at the covariate
values supplied in \code{newdata}.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
If \code{type = "class"} a vector with the predicted response
categories; if \code{type = "probs"} a matrix of probabilities for all
response categories at \code{newdata}.
}
\description{
Obtain class and probability predictions from a fitted adjacent
category logits model.
}
\details{
If \code{newdata} is omitted the predictions are based on the data
used for the fit.
}
\examples{

data("stemcell", package = "brglm2")

# Adjacent category logit (non-proportional odds)
fit_bracl <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
                   data = stemcell, type = "ML")
# Adjacent category logit (proportional odds)
fit_bracl_p <- bracl(research ~ as.numeric(religion) + gender, weights = frequency,
                    data = stemcell, type = "ML", parallel = TRUE)

# New data
newdata <- expand.grid(gender = c("male", "female"),
                       religion = c("liberal", "moderate", "fundamentalist"))

# Predictions
sapply(c("class", "probs"), function(what) predict(fit_bracl, newdata, what))
sapply(c("class", "probs"), function(what) predict(fit_bracl_p, newdata, what))

}
