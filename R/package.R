#' Dose-response Explorer
#'
#' There is a great variety of models developed for dose-response
#'    data, many of which have been implemented in the drc and DoseFinding
#'    packages. drexplorer combines both packages to aid the user to visually
#'    examine and compare how existing models perform on the data. Another
#'    important feature for drexplorer is to allow the user to identify outlier
#'    measurements and visually examine how these outliers affect the fitted
#'    model.
#' @docType package
#' @name drexplorer-package
#' @aliases drexplorer
#' @author Pan Tong, Kevin R Coombes
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit}}
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @importFrom epiR epi.ccc
#' @import stringr
#' @import fgui
#' @importFrom drc drm LL.2 LL.3 LL.3u LL.4 LL.5 W1.2 W1.3 W1.4 W2.2 W2.3 W2.4 BC.4 BC.5 LL2.2 LL2.3 LL2.3u LL2.4 LL2.5 AR.2 AR.3 
#' @importFrom DoseFinding fitMod emax sigEmax exponential quadratic linear linlog logistic 
NULL

## not working
# @importFrom drc drm LL.2 LL.3 LL.3u LL.4 LL.5 W1.2 W1.3 W1.4 W2.2 W2.3 W2.4 BC.4 BC.5 LL2.2 LL2.3 LL2.3u LL2.4 LL2.5 AR.2 AR.3 
# @importFrom DoseFinding fitMod emax sigEmax exponential quadratic linear linlog logistic 
