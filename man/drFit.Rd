\name{drFit}
\alias{drFit}
\title{
Fit a dose-response model
}
\description{
This function can fit various dose-response models by specifying a model name and package source (either drc or DoseFinding).
}
\usage{
drFit(drMat, modelName = "sigEmax", alpha = 0.01, fitCtr = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drMat}{
dose-response matrix. the first column being dosage and second column being response. controls are included by specifying dose=0
}
  \item{modelName}{
a dose-response model. For available models to be specified, see drModels().
}
  \item{alpha}{
a scalar for significance level. This specifies the significance level to identify outliers which will be excluded from model fitting. To include
all data, set alpha=1. 
}
  \item{fitCtr}{
A logic vector specifying whether to include the control points into the model fitting.
}
}
\details{
This is a wrapper to fit dose response models from drc and DoseFinding packages. 

}
\value{
It returns a drFit S4 object.
}

\author{
Kevin R Coombes (kcoombes@mdanderson.org), Pan Tong (nickytong@gmail.com)
}
\seealso{
\code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit-class}}
}
\examples{
data(ryegrass) # use the ryegrass data from drc package
# fit a sigmaEmax model without outlier removal. the controls are excluded from model fitting
fit.LL.3 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3", alpha=0.01, fitCtr=FALSE)
fit.LL.3u <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3u", alpha=0.01, fitCtr=FALSE)
fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.01, fitCtr=FALSE)
}
