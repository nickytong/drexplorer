\name{drModels}
\alias{drModels}
\title{
Show available dose-response models supported.
}
\description{
This function shows available models to be passed to drFit function as modelName.
}
\usage{
drModels(return=FALSE, verbose=TRUE)
}
\arguments{
  \item{return}{
whether to return the models in a list
}
  \item{verbose}{
whether to print out the models
}
}

\value{
Following is the print out:

Models implemented in drc package:
LL.2
LL.3
LL.3u
LL.4
LL.5
W1.2
W1.3
W1.4
W2.2
W2.3
W2.4
BC.4
BC.5
LL2.2
LL2.3
LL2.3u
LL2.4
LL2.5
AR.2
AR.3
MM.2
MM.3
Models implemented in DoseFinding package:
emax
emaxGrad
sigEmax
sigEmaxGrad
exponential
exponentialGrad
quadratic
quadraticGrad
betaMod
betaModGrad
linear
linearGrad
linlog
linlogGrad
logistic
logisticGrad
linInt
linIntGrad
More specific information can be found in the original packages.

}

\seealso{
\code{\link{NewmanTest}, \link{drOutlier}, \link{drFit}, \link{drFit-class}}

}
\examples{
drModels()
}