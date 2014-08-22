\name{fit_median_efect}
\alias{fit_median_efect}
\title{notice that the median-effect model is a special case of the sigma Emax model;
3 linear models are fitted; fixed ratio thus is required
mostly borrowed from chou8491.SSC; but also borrowed from CI_IIV2.SSC}
\usage{
  fit_median_efect(d1, d2, e, name1 = "Drug A",
    name2 = "Drug B", d2.d1, base = exp(1))
}
\arguments{
  \item{d1}{dose for drug 1}

  \item{d2}{dose for drug 2}

  \item{e}{corresponding response in the range [0, 1]}

  \item{name1}{drug 1 name}

  \item{name2}{drug 2 name}

  \item{d2.d1}{d2/d1, the fixed ratio}

  \item{base}{base of logarithm; disabled since the
  dose-response curve estimation also needs transformation}
}
\value{
  a list
}
\description{
  notice that the median-effect model is a special case of
  the sigma Emax model; 3 linear models are fitted; fixed
  ratio thus is required mostly borrowed from chou8491.SSC;
  but also borrowed from CI_IIV2.SSC
}

