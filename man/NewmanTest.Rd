\name{NewmanTest}
\alias{NewmanTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying outliers with the range to standard deviation ratio statistic.
}
\description{
This function implements the method described 
by D. Newman, Biometrika, 1939 to identify outliers.
}
\usage{
NewmanTest(ref, obs, alpha = 0.01, recursive = FALSE)
}
\arguments{
  \item{ref}{
a vector giving the reference values to estimate sigma
}
  \item{obs}{
a vector giving the observed values where potential outlier might exist
}
  \item{alpha}{
significance level. Only 0.01, 0.05 or 1 is allowed}
  \item{recursive}{
whether to recursively identify outliers. Currently not implemented
}
}
\details{
Given measurements from controls (no drug treated), we can compute the sample standard deviation (s). The range of responses from
treated samples (w) can be computed for a given dose level. Assuming the controls to have the same variation as the drug treated case,
the distribution of ratio statistic q=w/s can be derived and used to calculate if there is outliers in the treated responses as described 
by D. Newman, Biometrika, 1939.  

Note that this function works for a single dose level. When multiple dose levels exist, one need to repeatedly call this function to identify
outliers at each dose level or use the flagOutliers() function which is just a wrapper.
}
\value{
\item{indicator }{a logical vector specifying if the corresponding point is flagged as outlier}
}
\references{
Newman, D. (1939). The distribution of range in samples from a normal population, 
expressed in terms of an independent estimate of standard deviation. Biometrika, 31(1/2), 20-30.
}
\author{
Kevin R Coombes (kcoombes@mdanderson.org), Pan Tong (nickytong@gmail.com)
}

\seealso{
\code{\link{drOutlier}, \link{drModels}, \link{drFit}, \link{drFit-class}}
}
\examples{
set.seed(1)
x <- rnorm(10, 0, 1)
y <- c(rnorm(5, 0, 1), rnorm(1, 0, 1)+4)
# the last observation in y is an outlier
NewmanTest(x, y, alpha=0.05)
}
