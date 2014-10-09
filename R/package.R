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
NULL

## not working
# @importFrom drc drm summary predict
# @importFrom DoseFinding fitMod
