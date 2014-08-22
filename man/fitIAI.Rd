\name{fitIAI}
\alias{fitIAI}
\title{Estimate the interaction index as well as its confidence interval using delta method for fixed ratio design}
\usage{
  fitIAI(d1, d2, e, E = seq(0.05, 0.95, 0.005),
    name1 = "Drug A", name2 = "Drug B", alpha = 0.05)
}
\arguments{
  \item{d1}{dose for drug 1}

  \item{d2}{dose for drug 2}

  \item{e}{corresponding response in the range [0, 1]}

  \item{E}{a vector of responses (between 0 and 1) where
  IAI and confidence interval are to be computed from.}

  \item{alpha}{significance level of confidence interval}
}
\value{
  a data frame with columns IAI, IAI.low, IAI.up, E, dx1
  (corresponding dose of drug 1), dx2 (corresponding dose
  of drug 1),
}
\description{
  this code is extracted from the source code distributed
  at https://biostatistics.mdanderson.org/SoftwareDownload/
}

