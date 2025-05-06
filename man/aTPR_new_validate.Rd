\name{aTPR_new_validate}
\alias{aTPR_new_validate}
\title{Adjusted True Positive Rate Estimation}
\description{
  Estimates the adjusted true positive rate (aTPR) with missing data present, across specified thresholds using risk recalibration and density ratio adjustment, optionally with bootstrapped confidence intervals.
}
\usage{
aTPR_new_validate(data, covariate, controlvar, validatevar, groupvar, ref, response,
                  risk, taus, calmethod = c("logit"), cal.args = list(2, TRUE),
                  drmethod = c("logit"), dr.args = list(2), cv = FALSE,
                  se.boot = FALSE, bootsize = 1000, alpha = 0.05, quietly = TRUE)
}
\arguments{
  \item{data}{A \code{data.frame} containing original and recalibrated risks, group labels, and control/validation indicators.}
  \item{covariate}{A \code{data.frame} containing covariates \eqn{Z} used for calibration and modeling.}
  \item{controlvar}{Binary variable indicating control set membership.}
  \item{validatevar}{Binary variable indicating validation set membership.}
  \item{groupvar}{The column indicating group assignment (e.g., demographic group).}
  \item{ref}{Reference group for evaluating relative performance.}
  \item{response}{Column name indicating the observed outcome variable.}
  \item{risk}{Column name of the original risk score to be evaluated.}
  \item{taus}{A numeric vector of thresholds at which to compute TPR.}
  \item{calmethod}{Calibration method to apply. Currently only \code{"logit"} is supported.}
  \item{cal.args}{Additional arguments for calibration; for \code{"logit"}, this includes polynomial degree and logit transformation flag.}
  \item{drmethod}{Method for estimating density ratios. Options are \code{"logit"} or \code{"none"}.}
  \item{dr.args}{Arguments for the density ratio estimation method.}
  \item{cv}{Logical; if \code{TRUE}, cross-validation is used during density ratio estimation.}
  \item{se.boot}{Logical; if \code{TRUE}, bootstrap is used to compute standard error and confidence intervals.}
  \item{bootsize}{Number of bootstrap replicates to sample. Default is 1000.}
  \item{alpha}{Confidence level for bootstrap intervals. Default is 0.05.}
  \item{quietly}{Logical; suppress messages if \code{TRUE}. Default is \code{TRUE}.}
}
\value{
  A list containing:
  \item{aTPR}{A \code{data.frame} with adjusted TPR estimates across groups and thresholds.}
  \item{boot}{(Optional) A \code{data.frame} with bootstrapped aTPR estimates, included if \code{se.boot = TRUE}.}
}
\details{
  The function performs the following steps:
  \enumerate{
    \item Recalibrates the original risk scores using the specified calibration method.
    \item Estimates density ratios to adjust for distributional shifts across groups.
    \item Computes adjusted true positive rates (aTPR) for each group at the specified thresholds.
    \item Optionally performs bootstrapping to provide standard errors and confidence intervals.
  }

  This is useful for fairness evaluation and adjustment in predictive modeling where the calibration and performance across subgroups are of interest.
}
\import dplyr
\importFrom stats glm predict poly sd quantile
\import FairRisk
\seealso{
  \code{\link{calibrateRisk_new_validate}}, \code{\link{estDensityRatioCV}}, \code{\link{get_aTPR}}
}
\examples{
# Example (requires supporting data and functions):
# result <- aTPR_new_validate(data, covariate, controlvar = "control",
#                              validatevar = "validate", groupvar = "s", ref = 1,
#                              response = "y", risk = "gX", taus = seq(0.1, 0.4, 0.9))
}
