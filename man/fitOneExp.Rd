\name{fitOneExp}
\alias{fitOneExp}
\title{Fit multiple dose response models for a single drug in a single cell line (One experiement unit)}
\usage{
  fitOneExp(dat, drug = NA, cellLine = NA, unit = NA,
    models = c("LL.3", "LL.3u", "sigEmax", "logistic", "linlog"),
    cols = NA, log.d = TRUE,
    percentile = seq(0.1, 0.9, by = 0.1), alpha = 0.01,
    fitCtr = FALSE, interpolation = TRUE, plot = FALSE,
    transparency = 1, ...)
}
\arguments{
  \item{dat}{a 2-column data frame with first column as
  dose and second column response. Controls are decided by
  dose=0}

  \item{drug}{drug for this analysis}

  \item{cellLine}{cell line for this analysis}

  \item{unit}{unit of drug concentration}

  \item{models}{models the user may specify}

  \item{cols}{colors of dose-response curves}

  \item{log.d}{in computed ICx values, whether to return
  the log10 dose or the dose without this transformation}

  \item{percentile}{IC percentile}

  \item{alpha}{outlier significance level}

  \item{fitCtr}{whether the model is fitted with control
  data}

  \item{plot}{whether to draw the dose response curve}

  \item{...}{additional parameters to plotOneExp()}
}
\value{
  a list containing elements with:

  fits,

  models,

  cols,

  ICmat, IC matrix from all specified models as well as RSE
  and model name

  ICx, IC values from the best model

  datWithOutlierStatus the input data with outlier status
  appended

  bestModel the best model by RSE

  RSEs
}
\description{
  the dose response mode is usually fitted with log10 dose
  as it is best fitted by the dose-response curves. The
  computed ICx values however, can be in either log10 scale
  or the original scale
}

