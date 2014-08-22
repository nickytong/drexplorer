\name{drFit-class}
\docType{class}
\alias{drFit-class}
\alias{lines,drFit-method}
\alias{plot,drFit,ANY-method}
\alias{predict,drFit-method}

\title{Class \code{"drFit"}}
\description{
  The \code{drFit} class is designed to fit various dose-response models and exploring the effects of outliers using the method developed by D. Newman.
  }

\section{Creating Objects}{
  Although objects can be created directly using \code{new("drFit", ...)}, the most
  common usage will be to pass dose-response data to the 
  \code{drFit} function. 
}

\section{Slots}{
  \describe{
    \item{\code{fit}:}{A fit object by either drc or DoseFinding package. }
    \item{\code{fitDat}:}{The actual data matrix (after excluding outliers) used to fit a dose-response model.}
    \item{\code{originalDat}:}{The original dose response data, a matrix.}
    \item{\code{alpha}:}{A scalar for significance level. This specifies the significance level to identify outliers}
    \item{\code{fitCtr}:}{A logical value specifying whether to include the control points into the model fitting.}
    \item{\code{tag}:}{A string (either 'drc' or 'DoseFinding') tracking which package is used for model fitting.}
    \item{\code{info}:}{A list that holds information related to the model, i.e. Residual Standard Error (rse).}
  }
}
\section{Methods}{
  \describe{
    \item{lines(x = "drFit", col=5, lwd=2, show_points=FALSE, pch=16}{For a given fitted object (drFit), 
	add a predicted dose-response curve to a existing figure constructed by \code{plot}. col, lwd and \dots are for 
	generic plotting options. Updates: show_points is added to indicate if original data will be added in the form of points. It will have
	the same color as the line. pch will all be the same. Notice current implementation does not add outlier information from color and shape.
	}
    \item{plot(x = "drFit", pchs=c(16, 17, 15), cols=c(1, 2, 3), 
			col=4, lwd=2, addLegend=TRUE, xlab="Log10(Dose)", ylab="Relative viability", ylim=NA, xlim=NA, main, style='Full', bty='n', h=c(0.5), cex.main=1)}{Plot the original dose-response data and indicate
	potential outliers. It also add the fitted dose-response curve from the drFit object. pchs specify symbols to show different data points. In particular,
	pchs[1] is used to show regular points. pchs[2] is used to show outliers at 0.05 significance level and pchs[3] is used 
	to show outliers at 0.01 significance level.  Similarly, cols[1] is used to color regular points, cols[2] to color outliers at 0.05 significance level and
	cols[3] to color outliers at 0.01 significance level. col and lwd are used to specify the appearance of fitted curve.
	We later add the style option. Default is 'full' which means the outlier status as well as its legend will be shown; if style=='simple', the points and legend
	for outlier status will be removed; if style=='points', only points will be shown; this is useful if to compare multiple curves from different drugs. }
    \item{predict(object = "drFit", newData)}{Predict response for a given vector of dose levels (newData). Default is to predict response at observed dose levels.}
	 }
}
\references{
Newman, D. (1939). The distribution of range in samples from a normal population, 
expressed in terms of an independent estimate of standard deviation. Biometrika, 31(1/2), 20-30.
}
\author{
Kevin R Coombes (kcoombes@mdanderson.org), Pan Tong (nickytong@gmail.com)
}

\seealso{
\code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit}, \link{computeIC}}
}
\examples{

###### needs to be commented out
data(ryegrass) # use the ryegrass data from drc package
# fit a sigmaEmax model without outlier removal. the controls are excluded from model fitting
fit.LL.3 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3", alpha=0.01, fitCtr=FALSE)
fit.LL.3u <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3u", alpha=0.01, fitCtr=FALSE)
fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.01, fitCtr=FALSE)
predict(fit.sigEmax)
predict(fit.LL.3u)

###
plot(fit.LL.3, main='', col=4, lwd=2)
lines(fit.LL.3u, col=5, lwd=2)
lines(fit.sigEmax, col=6, lwd=2)
legend("bottomleft", c('LL.3', 'LL.3u', 'sigEmax'), col=4:6, lwd=3)

}
\keyword{classes}
