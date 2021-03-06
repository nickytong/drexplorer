% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drexplorer.R
\docType{class}
\name{drFit-class}
\alias{drFit-class}
\title{Class \code{"drFit"}}
\description{
Class \code{"drFit"}
}
\section{Slots}{

\describe{
\item{\code{fit}}{A fit object by either drc or DoseFinding package.}

\item{\code{fitDat}}{The actual data matrix (after excluding outliers) used to fit a dose-response model.}

\item{\code{originalDat}}{The original dose response data, a matrix.}

\item{\code{alpha}}{A scalar for significance level. This specifies the significance level to identify outliers}

\item{\code{fitCtr}}{A logical value specifying whether to include the control points into the model fitting.}

\item{\code{tag}}{A string (either 'drc' or 'DoseFinding') tracking which package is used for model fitting.}

\item{\code{modelName}}{model name of the fitted model}

\item{\code{nPar}}{number of parameters for the fitted model}

\item{\code{info}}{A list that holds information related to the model, i.e. Residual Standard Error (rse).}
}}
\section{Creating Objects}{
 
Although objects can be created directly using \code{new("drFit", ...)}, the most common usage will be to pass dose-response data to the \code{drFit} function.
}
\seealso{
\code{\link{drOutlier}, \link{drModels}, \link{drFit}, \link{drFit-class}}
}

