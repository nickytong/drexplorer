\name{drOutlier}
\alias{drOutlier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying outliers for dose-response data
}
\description{
A wrapper for the NewmanTest() function. It repeatedly call NewmanTest() for each dose level. 
}
\usage{
drOutlier(drMat, alpha = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drMat}{
dose-response matrix. the first column being dosage and second column being response. controls are included by specifying dose=0
}
  \item{alpha}{
a scalar for significance level. Only 0.01, 0.05 and 1 are allowed. alpha=1 is included for comparability issue. In this case,  no outliers will be identified. 
}
}
\value{
\item{indicator }{a vector with length nrow(drMat) specifying if a data point is outlier. Control points always have status FALSE}
}
\author{
Kevin R Coombes (kcoombes@mdanderson.org), Pan Tong (nickytong@gmail.com)
}
\seealso{
\code{\link{NewmanTest}, \link{drModels}, \link{drFit}, \link{drFit-class}}
}
\examples{
data(ryegrass) # use the ryegrass data from drc package
drOutlier(drMat=ryegrass[, c(2, 1)], alpha=0.05)

}
