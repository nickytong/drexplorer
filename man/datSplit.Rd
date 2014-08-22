\name{datSplit}
\alias{datSplit}
\title{Split a drug-response data for the template from Suk Young}
\usage{
  datSplit(dat, colNames, by, minReplicate = 3,
    verbose = TRUE)
}
\arguments{
  \item{dat}{ input data with multiple cell lines for one drug}
  \item{colNames}{colNames to name to columns}

  \item{by}{string for the colName used for split}

  \item{minReplicate}{if a subset does not have nrow   exceeding this, throw an error}
  \item{verbose}{ whether to produce warning message }
}
\value{
  a list of the split data frames
}
\description{
  This split the data into a list so that each element is a
  drug-Cell Line combination for further analysis
}

