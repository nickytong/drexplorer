% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drexplorer_hillNIC60.R
\name{plot.nci60Fit}
\alias{plot.nci60Fit}
\title{plot method for nci60Fit class}
\usage{
\method{plot}{nci60Fit}(fit, xlab = "Dose", ylab = "Relative growth",
  main = "Fit of NCI60 method", xlim = NULL, ylim = NULL, cex.main = 1,
  cex.axis = 1, cex.lab = 1, axes = TRUE, pcol = "black", pch = 1,
  lcol = "black", lty = 1, ltyh = 1, lwd = 2, type = c("line",
  "points"), h = c(-0.5, 0, 0.5), ...)
}
\arguments{
\item{h}{horizontal line to add indicating e.g. GI (h=0.5), TGI (h=0), LD50 (h=-0.5)}

\item{...}{additional parameters, not implemented}

\item{x}{a nci60Fit object as returned by nci60Fit}
}
\description{
plot method for nci60Fit class
}

