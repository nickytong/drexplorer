#' Fit multiple dose response models for a single drug in a single cell line (One experiement unit)
#' 
#' the dose response mode is usually fitted with original dose (not log10 transformed). The computed ICx
#' values however, can be in either log10 scale or the original scale.  
#'
#' @param dat a 2-column data frame with first column as dose and second column response. Controls are decided by dose=0
#' @param drug drug for this analysis
#' @param cellLine cell line for this analysis
#' @param unit unit of drug concentration
#' @param models models the user may specify 
#' @param cols colors of dose-response curves
#' @param log.d in computed ICx values, whether to return the log10 dose or the dose without this transformation
#' @param percentile IC percentile
#' @param alpha outlier significance level
#' @param fitCtr whether the model is fitted with control data
#' @param standardize whether to standardize (scale) the data based on control points. This should be disabled when no control data is supplied
#' @param interpolation whether calculate ICx through interpolation or through root finding. This parameter is passed to computeIC() function. When interpolation=TRUE,
#' ICx value will be bounded by observed minimum and maximum dose; otherwise, ICx is calculated through root finding and thus can be outside this boundary (extrapolated). 
#' @param plot whether to draw the dose response curve
#' @param transparency a value between 0 and 1 specifying transparency through alpha blending; transparency=0 means totally transparent.
#' @param ... additional parameters to plotOneExp()
#' @return a list containing elements with:
#'
#' fits, 
#'
#' models, 
#'
#' cols,
#'
#' ICmat, IC matrix from all specified models as well as RSE and model name
#'
#' ICx, IC values from the best model  
#'
#' datWithOutlierStatus the input data with outlier status appended
#'
#' bestModel the best model by RSE 
#'
#' RSEs 
#'
#' @export
fitOneExp <- function(dat, ### data format specific to: i.e. ExportToR 2013 07 01 A549- HCC827 - 633.xlsx
					drug=NA,	
					cellLine=NA,
					unit=NA, # dose unit					
					 models=c("LL.3", "LL.3u", "sigEmax", "logistic", "linlog"), # models the user may specify
					 cols=NA, # colors of dose-response curves
					 log.d=TRUE,
					 percentile=seq(0.1, 0.9, by=0.1), # IC percentile
					 alpha=0.01, # outlier significance level
					 fitCtr=FALSE, # whether the model is fitted with control data
					 standardize=TRUE, # whether to standardize the data based on control
					 interpolation=TRUE, # interpolation for IC or not
					plot=FALSE, transparency=1, ...){
	require(drexplorer)
	if(is.na(cols[1])){ # use default colors when not specified or inappripriate
		if(length(models)<=9){
			cols <- alpha(brewer.pal(9, "Set1")[seq_along(models)], alpha=transparency)
		} else {
			cols <- alpha(rainbow(length(models)), alpha=transparency)
		}
	}
	if(length(cols)!=length(models) & !is.na(cols[1])) cols <- rep(cols[1], length(models))
	if(is.na(drug)){
		stop('Drug name not specified\n')
	}
	if(is.na(cellLine)){
		stop('Cell Line name not specified\n')
	}
	#browser()
	if(ncol(dat)>2) stop('Only accept data with 2 columns for the dat argument!\n')
	#browser()	
	dose <- dat[, 1]
	#resp <- dat[, 2]
	drMat <- dat
	#indCtrl <- which(round(dose, 5)==0) # dose = 0 at precision 1e-5 would be accepted as control
	indtrt <- which(round(dose, 15)!=0) # dose = 0 at precision 1e15 (not enough for 1e-5) would be accepted as control
	#browser()
	dmin <- min(dose[indtrt], na.rm=TRUE)
	dmax <- max(dose[indtrt], na.rm=TRUE)
	indicator <- drOutlier(drMat=drMat, alpha=alpha) 
	# this does not work: set list element to NULL actually removes it
	#fits <- vector('list')
	#for(i in 1:length(models)){
	#		tmfit <- try(drFit(drMat=drMat, modelName = models[i], alpha=alpha, fitCtr=fitCtr, standardize=standardize), silent=TRUE)
	#		## some model fails to model the data and numerically not computable by the original packages
	#		if(class(tmfit)!='try-error') {
	#			fits[[i]] <- tmfit
	#		} else {
	#			#fits[[i]] <- NULL: this makes length of fits != models
	#			fits[[i]] <- list(NULL)
	#			warning(sprintf('Drug: %s CellLine: %s Model: %s failed!\n', drug, cellLine, models[i]))
	#		}
	#}
	# no foreach 
	fits <- lapply(models, function(x){
		tmfit <- try(drFit(drMat=drMat, modelName = x, alpha=alpha, fitCtr=fitCtr, standardize=standardize), silent=TRUE)
			## some model fails to model the data and numerically not computable by the original packages
			if(class(tmfit)!='try-error') {
				res <- tmfit
			} else {
				#fits[[i]] <- NULL: this makes length of fits != models
				res <- NULL
				# bug here
				#warning(sprintf('Drug: %s CellLine: %s Model: %s failed!\n', drug, cellLine, models[i]))
				warning(sprintf('Drug: %s CellLine: %s Model: %s failed!\n', drug, cellLine, x))
			}
			res
	})
	#fits <- foreach(i=1:length(models)) %do% {
	#		tmfit <- try(drFit(drMat=drMat, modelName = models[i], alpha=alpha, fitCtr=fitCtr, standardize=standardize), silent=TRUE)
	#		## some model fails to model the data and numerically not computable by the original packages
	#		if(class(tmfit)!='try-error') {
	#			res <- tmfit
	#		} else {
	#			#fits[[i]] <- NULL: this makes length of fits != models
	#			res <- NULL
	#			warning(sprintf('Drug: %s CellLine: %s Model: %s failed!\n', drug, cellLine, models[i]))
	#		}
	#		res
	#}
	#browser()
	names(fits) <- models
	indSuccess <- which(!sapply(fits, is.null))
	###
	# on 2014/04/16: remove sigEmax when its coef has NA
	###
	# a bug in sigEmax model: the model fits correctly but Emax value might be NA and makes it impossible to calculate IC value!
	# now we detect this and remove it from successful model
	if('sigEmax' %in% names(indSuccess)) {
		indCheck <- indSuccess['sigEmax']
		Coef <- coef(fits[[indCheck]]@fit)
		if(any(is.na(Coef)))
			indSuccess <- indSuccess[indSuccess!=indCheck]
	}
	if(length(indSuccess)>0) {
			rMax <- max(sapply(fits[indSuccess], function(x) max(x@fitDat[, 2])))
			RSEs <- sapply(fits[indSuccess], function(x) x@info$RSE)
			indBest <- indSuccess[which.min(RSEs)] # find best model by RSE
			bestModel <- models[indBest]
	} else {
		stop(sprintf('None of the specified models can be fitted in cell line: %s drug: %s\n\tSpecified models are: %s', cellLine, drug, str_c(models, collapse=', ')))
		#message(sprintf('None of the specified models can be fitted in cell line: %s drug: %s\n\tSpecified models are: %s', cellLine, drug, str_c(models)))
		#res <- NA
		#attr(res, 'CellLine') <- cellLine
		#attr(res, 'drug') <- drug
		#attr(res, 'models') <- models
		#attr(res, 'call') <- as.list(match.call())
		#browser()
		#return(res)
	}
	#zz = computeICvariance(fits[indSuccess][[2]], percent=percentile)
	#zz = computeIC(fits[indSuccess][[2]], percent=percentile)
	# added variance estimation of IC values on 2015/07/21
	nPar <- sapply(fits[indSuccess], function(x) x@nPar)
	ICmat0_variance <- t(sapply(fits[indSuccess], computeICvariance, percent=percentile))
	ICmat0 <- t(sapply(fits[indSuccess], computeIC, percent=percentile, log.d=log.d, interpolation=interpolation))
	#AUCmat_untrnsf <- t(sapply(fits[indSuccess], computeAUC, dmin=dmin, dmax=dmax, islogd=F)) # AUC, AUC0, AUCs in original scale
	AUCmat_trnsf <- t(sapply(fits[indSuccess], computeAUC, dmin=log10(dmin), dmax=log10(dmax), islogd=T)) # AUC, AUC0, AUCs in log10 sacle
	#getAUC(fit_sigEmax_alpha_o5, dmin=-0.027, dmax=1.477, islogd=T)
	IC50 <- sapply(fits[indSuccess], computeIC, percent=0.5, log.d=log.d, interpolation=interpolation)
	names(IC50) <- models[indSuccess]
	# make sure to use: models[indSuccess] since RSEs is only for successful model
	ICmat <- data.frame(Drug=drug, CellLine=cellLine, Model=models[indSuccess], isBestModel=(RSEs==min(RSEs, na.rm=TRUE)), RSE=RSEs, ICmat0, AUCmat_trnsf)
	ICmat_variance <- data.frame(Drug=drug, CellLine=cellLine, Model=models[indSuccess], isBestModel=(RSEs==min(RSEs, na.rm=TRUE)), ICmat0_variance)
	datWithOutlierStatus <- data.frame(dat, isOutlier=indicator)
	# append min and max dose so that the user can use this to truncate the predicted value
	#ICx <- ICmat[indBest, ] # this has a mismatch! indBest is absolute index! ICmat removes the failures! Use name for tracking!!!
	ICx <- ICmat[names(indBest), ]
	ICx_variance <- ICmat_variance[names(indBest), ]
	ICx$minLog10Dose <- log10(dmin)
	ICx$maxLog10Dose <- log10(dmax)
	ICx$unit <- unit
	## add experiment setup info
	resPrepDR <- prepDRdat(drMat, alpha=alpha, fitCtr=fitCtr, standardize=standardize)
	info_experimentSetup <- prepDRdat2expInfo(resPrepDR)
	res <- list(fits=fits, 
		indSuccess=indSuccess, indBest=indBest, 
		models=models, cols=cols, unit=unit, 
		ICmat=ICmat, ICx=ICx, ICmat_variance=ICmat_variance, ICx_variance=ICx_variance, nParBest=nPar[names(indBest)], nPar=nPar, IC50=IC50,  
		datWithOutlierStatus=datWithOutlierStatus,
		bestModel=bestModel, RSEs=RSEs, drug=drug,
		cellLine=cellLine, info_experimentSetup=info_experimentSetup
		)
	#browser()
	if(plot){
		plotOneExp(res, ...)
	}	
	res
}

#' plot dose response curve from fitOneExp() result
#'
#' @param fitRes return value from fitOneExp()
#' @param ind2plot index for the models that will be plotted; default is NA which leads to all curves available; when specified as 'best', the best model is selected 
#' @param cols color for the lines. If cols is of length 1, then all lines (representing different models) have the same color as specified; if cols has a length larger than 1,
#'  the function further checks if this length equals the number of fitted models. If this is the case, cols specifies colors for all models. Otherwise, default color is used. 
#' @param pcols cols for points, similar to cols
#' @param type the plot type. The user can specify more flexibly by specifying a vector of entities to plot including:
#' either plot or line; when specified as line, it will only adds to an existing figure; When length(ind2plot)>1, type will be reset to plot
#' 	which means the first curve will be made with plot() and additional ones with lines()
#' 	plot: to start a new figure
#' 	line: dose response curve
#' 	points: observed data points
#' 	control: control data points
#' 	legend: legend for outlier status
#' @param h horizontal line added to the figure, i.e. indicating IC50, IC70
#' @param tag tag before main
#' @param main main
#' @param cex.main cex.main to adjust main title size
#' @param ylim ylim
#' @param xlim xlim
#' @param xlab xlab
#' @param ylab ylab
#' @param style deprecated. use type instead
#' @param show whether to show RSE ('RSE'), IC50 ('IC50'), both ('both') or nothing ('None') in addition to the best model as legend
#' @param cexLegend legend cex
#' @param posLegend position of legend for model information 
#' @param showTopN if specified show best N model in figure to avoid busy plotting; otherwise show all successful models. 
#' @param lwd line width for the curves
#' @param lwd line width for the curves
#' @param lty line type for the curves
#' @param cex.axis cex for axis annotation
#' @param cex.lab cex for axis label
#' @param addOutlierLegend whether to add legend for outlier status
#' @param axes whether to add axes. only effective if type='plot'
#' @export
plotOneExp <- function(fitRes, ind2plot=NA, cols=NA, pcols=NA, type='plot', style='full', h=NULL, tag=NA, main=NA, 
	cex.main=1, xlab=NA, ylab=NA, ylim=NA, xlim=NA, show='both', cexLegend=NA, posLegend='bottomleft', showTopN=NA, lwd=2, lty=1, cex.axis=1, cex.lab=1, 
	addOutlierLegend=TRUE, axes=TRUE){
	#browser()
	# calculate the actual color to be used
	# notice cols here is different from fitRes$models: it ensures col_use have equal length as length(fitRes$models)
	# therefore, cols should always have length of length(fitRes$models)
	if(is.na(cols[1]) | length(cols)!=length(fitRes$models)) {
		# when col is not specified, use the col in the model
		warning('cols do not have the same length as number of models; using the cols attached to the fitted object!\n')
		col_use <- fitRes$cols
	} else {
		col_use <- cols # use the col specified when length matches
	}
	plotBest <- ifelse(!is.na(ind2plot[1]) & ind2plot[1]=='best', TRUE, FALSE)
	if(plotBest & length(cols)==1 & !is.na(cols)){
		# the only case to specify incompatible cols is when only plot the best model: populate all cols as the only col specified to get it done
		col_use <- rep(cols, length(fitRes$models))
	}
	if(is.na(pcols[1])) pcols <- col_use
	#attach(fitRes) # this is a place for bug!
	with(fitRes, {
	# draw dots
	## when the user does not specify which index of the models to plot, plot all available ones; of course, the user can specify the one for the best model
	RSEs <- sapply(fits[indSuccess], function(x) x@info$RSE)
	RSEs_aug <- sapply(fits, function(x) ifelse(is.null(x), NA, x@info$RSE)) # length matched with fits
	indBest <- indSuccess[which.min(RSEs)] # find best model by RSE
	bestModel <- models[indBest]
	if(is.na(xlab)) xlab <- 'Log10(Dose)'
	#xlab <- sprintf('Log10(Dose) of %s', drug)
	if(is.na(ylab)) ylab <- 'Relative viability'
	#browser()
	if(is.na(tag)) tag <- sprintf('Drug: %s\nCell line:%s', drug, cellLine)
	if(is.na(main)) main <- sprintf('%s\nBest model:%s', tag, bestModel)
	if(is.na(ind2plot[1])){
		if(length(indSuccess)==0) {
			ind2plot <- NULL # no ind to plot
		} else {
			if(!is.na(showTopN)){
				if(showTopN<=length(indSuccess)){
					# just extract top models from successful ones
					ind2plot <- indSuccess[order(RSEs, decreasing=FALSE)[1:showTopN]]
				} else {
					ind2plot <- indSuccess[order(RSEs, decreasing=FALSE)] # in case topN > successful models, show all successful one
				}
			} else {
				ind2plot <- indSuccess[order(RSEs, decreasing=FALSE)]
			}
		}	
	}
	if(plotBest){ # just plot the best model
		ind2plot <- which(models==bestModel)
	}
	#browser()
	if(length(ind2plot)>1)  {# when there are multiple curves to plot, type must be plot and additional lines 
		type=c(type, 'plot') # so that to open a new figure
	}
	if(plotBest){
		type=c(type, 'line')
	}
	if(length(pcols)==1){
		pcols <- rep(pcols, length(col_use)) # in case only one color is specified
	}
	#browser()
	if('plot' %in% type){ 
		rMax <- getMaxResponse(fits[indSuccess]) # max(sapply(fits[indSuccess], function(x) max(x@fitDat[, 2]))) # maxum response value to control ylim
		if(is.na(ylim[1]))
			ylim <- c(0, max(1.2,rMax))
		#browser() # pcol=pcols[ind2plot[1]], # unable to do at this moment
		# cols: 3 elements, first one is for regular points
		#browser()
		plot(fits[[ind2plot[1]]], col=col_use[ind2plot[1]], cols=c(pcols[ind2plot[1]], 2, 3), lwd=lwd, lty=lty, main=main, type=type, 
				xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, bty='n', h=h, cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab, axes=axes)		
	}
	#browser()
	if ('line' %in% type){
		#browser()
		#show_points <- ifelse('points' %in% type, TRUE, FALSE)
		lines(fits[[ind2plot[1]]], col=col_use[ind2plot[1]], pcol=pcols[ind2plot[1]], lwd=lwd,type='line', lty=lty)
	}
	if ('points' %in% type){
		#browser()
		#show_points <- ifelse('points' %in% type, TRUE, FALSE)
		lines(fits[[ind2plot[1]]], col=col_use[ind2plot[1]], pcol=pcols[ind2plot[1]], lwd=lwd,type='points', lty=lty)
	}
	if ('sem' %in% type){
		#show_points <- ifelse('points' %in% type, TRUE, FALSE)
		lines(fits[[ind2plot[1]]], col=col_use[ind2plot[1]], pcol=pcols[ind2plot[1]], lwd=lwd,type='sem', lty=lty)
	}
	#browser()
	# draw all remaining models
	if(length(ind2plot)>1) {
			#browser()
			for(i in ind2plot[-1]) {
				lines(fits[[i]], col=col_use[i], pcol=pcols[i], lwd=lwd, lty=lty)
			}
	}
	#models_ <- models
	#if(any(sapply(fits, is.null))) models_[sapply(fits, is.null)] <- paste(models[sapply(fits, is.null)], ': failed', sep='')
	#legend("bottomleft", models_, col=cols, lwd=3, bty='n')
	models_more <- models
	if(show=='IC50'){
		if(is.na(cexLegend)) cexLegend <- 1
		models_more <- str_c(models, signif(IC50, 3), sep=', IC50=')
	}
	if(show=='RSE'){
		if(is.na(cexLegend)) cexLegend <- 1
		models_more <- str_c(models, signif(RSEs_aug, 4), sep=', RSE=')
	}
	if(show=='both'){
		if(is.na(cexLegend)) cexLegend <- 0.8
		models_more <- str_c(models, ', IC50=', signif(IC50, 3), ', RSE=', signif(RSEs_aug, 4), sep='')
	}
	if(any(sapply(fits, is.null))) models_more[sapply(fits, is.null)] <- paste(models[sapply(fits, is.null)], ': failed', sep='')
	models_show <- models_more[ind2plot]
	#browser()
	if(show!='None') { # None disable the legend
		if(is.na(ind2plot[1])) { # by default, show legend of all models when ind2plot is not specified
			legend(posLegend, models_show, col=col_use, lwd=lwd*2, bty='n', cex=cexLegend)
		} else { # when ind2plot=='best', need to match with customized color
			legend(posLegend, models_show, col=col_use[ind2plot], lwd=lwd*2, bty='n', cex=cexLegend)
		}	
	}
	#browser()
	#detach("fitRes")
	})
}

getMaxResponse <- function(fits){
	max(sapply(fits, function(x) max(x@fitDat[, 2]))) # maxum response value to control ylim
}
