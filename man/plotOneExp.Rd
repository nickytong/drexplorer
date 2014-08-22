\name{plotOneExp}
\alias{plotOneExp}
\title{plot dose response curve from fitOneExp() result}
\usage{
  plotOneExp(fitRes, h = c(0.3, 0.5), main = NA, xlab = NA,
    ylab=NA, show='both', cexLegend=NA, showTopN=NA)
}
\arguments{
  \item{fitRes}{return value from fitOneExp()}

  \item{h}{horizontal line added to the figure, i.e.
  indicating IC50, IC70}

  \item{main}{main}

  \item{xlab}{xlab}

  \item{ylab}{ylab}
  \item{show}{either RSE, IC50 or both indicating if we need to show RSE and/or IC50 in addition to the model}
  \item{cexLegend}{legend cex}
  \item{showTopN}{if specified show best N model in figure to avoid busy plotting; otherwise show all successful models.}
}
\description{
  plot dose response curve from fitOneExp() result
}

