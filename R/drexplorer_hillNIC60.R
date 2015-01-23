
# Hill equation to return y=f(Dose). HS>0 is assumed 
#' the Hill equation function
#'
#' @param Dose the concentration of drug (original scale)
#' @param Einf the bottom asymptote of response (e.g. at drug concentration of Inf)
#' @param E0 the top asymptote of response (e.g. at drug concentration of 0)
#' @param logEC50 natural logarithm of EC50
#' @param HS Hill slope
#' @return the response at dose=Dose (between 0 and 1)
#' @export
fHill <- function(Dose, Einf, E0, logEC50, HS) {
	res <- Einf + (E0-Einf) /(1+(Dose/exp(logEC50))^HS)
	res
}

#' fit Hill equation on dose response data
#'
#' This function implements S3-style OOP. Methods such as predict, plot and lines are available.
#' Notice that control data (dose=0) is usually not fitted (supplied) in Hill equation.
#'
#' @param d dose original dose
#' @param y response relative viability 
#' @return it returns a hillFit object consisting a vector of estimations for 'E0', 'Einf', 'EC50', 'HS', 'Emax', 
#' 'MSE' (Mean Squared Error), 'Rsq' (R squared), 'RSE' (Residual Standard Error)
#' @export
hillFit <- function(d, y){
	if(length(d)!=length(y)){
		stop('Please specify equal length in d (dose) and y (response)!\n')
	}
	dat <- data.frame(dose=d, response=y)
	# this might fail due to convergence
	# need to add boundary: in the order of start
	# EC50 must between dmin and dmax
	lb <- c(0, 0, log(min(dat$dose)), 1e-5)
	ub <- c(2, 1, log(max(dat$dose)), 1e5)
	# warnOnly=TRUE: always give an result although it might be non-convergent ---> this will always return a result but Rsq can be really poor
	fit <- try(nls(y ~ Einf + (E0-Einf) /(1+(dose/exp(logEC50))^HS), start=list(E0=1, Einf=0, logEC50=log(median(dat$dose)), HS=1), 
		algorithm="port", trace=FALSE, lower=lb, upper=ub, control=list(maxiter=50000, warnOnly=TRUE), data=dat), silent=TRUE)
	#browser()
	res <- rep(NA, 8)
	names(res) <- c('E0', 'Einf', 'EC50', 'HS', 'Emax', 'MSE', 'Rsq', 'RSE')
	if(!inherits(fit, 'try-error')){
		cc <- coef(fit)
		rss <- sum(residuals(fit)^2) # residual sum of squared
		tss <- sum((y-mean(y, na.rm=T))^2, na.rm=T)
		MSE <- rss/length(y)
		Rsq <- 1-rss/tss
		df <- sum(!is.na(y))-4 # four parameters in for Hill model
		rse <- sqrt(rss/df)# residual standard error: sqrt(RSS/df)
		res['EC50'] <- exp(cc['logEC50'])
		res['E0'] <- cc['E0']
		res['Einf'] <- cc['Einf']
		res['HS'] <- cc['HS']
		res['MSE'] <- MSE
		res['Rsq'] <- Rsq
		res['RSE'] <- rse
		res['Emax'] <- predict(fit, list(dose=max(d, na.rm=T))) # response at maximum dose using nls predict method
	}	
	attr(res, 'class') <- c('numeric', 'hillFit')
	attr(res, 'fitDat') <- dat
	attr(res, 'fit') <- fit
	#browser()
	res
}
#fit = hillFit(d=mydat$Dose, y=mydat$y)

#' predict method for hillFit class
#' @method predict hillFit
#' @param x a hillFit object as returned by hillFit
#' @param newData dose in original scale to be used for prediction of response
#' @return predicted response at specified dose
#' @export
predict.hillFit <- function(fit, newData=NULL){
	fitDat <- Getter(fit, 'fitDat')
	if(is.null(newData)) newData <- fitDat$dose 
	dose <- newData
	y <- fHill(Dose=dose, Einf=fit['Einf'], E0=fit['E0'], logEC50=log(fit['EC50']), HS=fit['HS'])
	#y <- fit['Einf'] + (fit['E0']-fit['Einf']) /(1+(dose/fit['EC50'])^fit['HS'])
	y
}

getPlotDat_hill <- function(fit){
	fitDat <- Getter(fit, 'fitDat')
	gg <- format_grid(fitDat$dose)
	top <- gg$top
	bot <- gg$bot
	xGrid <- gg$xGrid
	y <- predict(fit, newData=xGrid)
	## too many points and leads to a large pdf: use a subset of the points
	ind1 <- which(diff(log10(xGrid))>1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
	ind2 <- floor(seq(max(ind1)+1, length(xGrid), length.out=1000))
	indSel <- c(ind1, ind2)
	list(fitDat=fitDat, y=y, indSel=indSel, xGrid=xGrid)
}
#' plot method for hillFit class
#' @method plot hillFit
#' @param x a hillFit object as returned by hillFit
#' @param ... additional parameters, not implemented
#' @export
plot.hillFit <- function(fit, xlab="Log10(Dose)", ylab="Relative viability", main='Fit of Hill equation', 
	xlim=NULL, ylim=NULL, cex.main=1, cex.axis=1, pcol='black', lcol='black', lwd=2, ...){
	plL <- getPlotDat_hill(fit)
	fitDat <- plL$fitDat
	xGrid <- plL$xGrid
	indSel <- plL$indSel
	y <- plL$y
	if(is.null(ylim)) ylim <- range(pretty(fitDat$response))
	if(is.null(xlim)) xlim <- range(pretty(log10(fitDat$dose)))		
	with(fitDat, plot(log10(dose), response, 
		xlab=xlab, ylab=ylab, main=main, xlim=xlim, ylim=ylim, cex.main=cex.main, cex.axis=cex.axis))
	#browser()
	lines(log10(xGrid)[indSel], y[indSel], col=lcol, lwd=lwd)
	#browser()
} 
#plot.hillFit(fit)

#' lines method for hillFit class
#' @method lines hillFit
#' @param x a hillFit object as returned by hillFit
#' @param col line color	   
#' @param lwd line width  
#' @param show_points whether to add points for dose-response paires (dose not 0) 
#' @param pcol color for points; only relevant if show_points is TRUE
#' @param pch pch for points; only relevant if show_points is TRUE	   
#' @param ... additional parametrs passed to generic lines() function
#' @export
lines.hillFit <- function(fit, col=1, lwd=2, show_points=FALSE, pcol='black', pch=16, ...){
	plL <- getPlotDat_hill(fit)
	fitDat <- plL$fitDat
	xGrid <- plL$xGrid
	indSel <- plL$indSel
	y <- plL$y
	lines(log10(xGrid)[indSel], y[indSel], col=col, lwd=lwd)
	if(show_points){
		with(fitDat, points(log10(dose), response, col=pcol, pch=pch))
	}
	#browser()
} 
#lines(fit)

#' implementation of NCI60 method (GI50, TGI, LC50) for dose response data
#'
#' This function implements S3-style OOP. Methods such as predict, plot and lines are available.
#' Notice that here the response is assumed to be calculated as in: http://dtp.nci.nih.gov/branches/btb/ivclsp.html. The response has range
#' -1 to 1. 
#' The original NCI60 method applies linear interpolation to estimate GI50, TGI and LC50. Tong, F. P. improved this by using 4-parameter
#' logistic model, which is the LL.4 model in drexplorer (from drc package). For resistant cell lines, even the LL.4 cannot be fitted since
#' the curve is essentially flat. Here we fit multiple models throught the model selection framework in drexplorer and use the best model to estimate GI50, TGI and LC50.
#' This will almost give an estimate for all extreme data without giving NA.
#'
#' @param d dose original dose
#' @param y response percentage of growth inhibition as defined by NCI60 method
#' @param interpolation whether to use interpolatuon to estimat GI50, TGI and LC50.
#' @param log.d whether to return log10 transformed value for GI50, TGI and LC50 (the dose concentration)
#' @return a nci50Fit object consisting different attributes. It includes a vector with GI50, TGI and LC50. 
#' @export
#' @references Shoemaker, R. H. The NCI60 Human Tumour Cell line Anticancer Drug Screen. Nature Reviews, 6: 813-823, 2006.
#' @references Tong, F. P. (2010). Statistical Methods for Dose-Response Assays.
nci60Fit <- function(d, y, interpolation=TRUE, log.d=FALSE){
	#browser()
	dat <- data.frame(dose=d, response=y)
	#fit <- try(drFit(drMat=dat, modelName='LL.4', standardize=F), silent=TRUE)
	fitL <- fitOneExp(dat=dat, standardize=F, drug='', cellLine='', unit='', models=drModels(return=TRUE, verbose=FALSE)[['recommendedModels']], alpha=1, interpolation=TRUE)
	fitBest <- fitL$fits[[fitL$indBest]]
	fit <- fitBest
	if(!inherits(fit, 'try-error')){
		#m1 <- drm(y~d, fct = LL.4())
		res <- findDoseGivenResponse(fit, response=c(0.5, 0, -0.5), interpolation=interpolation, log.d=log.d)
	} else {
		res <- rep(NA, 3)
	}
	base <- ifelse(log.d, 10, 1)
	names(res) <- c('GI50', 'TGI', 'LC50')
	attr(res, 'class') <- c('numeric', 'nci60Fit')
	attr(res, 'fitDat') <- dat
	attr(res, 'fit') <- fit
	attr(res, 'base') <- base
	#browser()
	res
}
#fit2 <- nci60Fit(d=mydat$Dose, y=mydat$Growth/100)

#fit <- drFit(drMat=cbind(mydat$Dose, mydat$Growth/100), modelName='LL.4', standardize=F)
#findDoseGivenResponse(fit, response=0.5, interpolation=F)


#' predict method for nci60Fit class
#' @method predict nci60Fit
#' @param x a nci60Fit object as returned by nci60Fit
#' @param newData dose in original scale to be used for prediction of response
#' @return predicted response at specified dose
#' @export
predict.nci60Fit <- function(fit, newData=NULL){
	fitDat <- Getter(fit, 'fitDat')
	fitDR <- Getter(fit, 'fit') #drFit object
	if(is.null(newData)) newData <- fitDat$dose
	#browser()
	if(!inherits(fitDR, 'try-error')){ 
		y <- predict(fitDR, newData=newData)
	} else {
		y <- rep(NA, length(newData))
	}
	#browser()
	y
}
# predict(fit2)

getPlotDat_nci60 <- function(fit){
	fitDat <- Getter(fit, 'fitDat')
	gg <- format_grid(fitDat$dose)
	top <- gg$top
	bot <- gg$bot
	xGrid <- gg$xGrid
	#browser()
	y <- predict(fit, newData=xGrid)
	## too many points and leads to a large pdf: use a subset of the points
	ind1 <- which(diff(log10(xGrid))>1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
	ind2 <- floor(seq(max(ind1)+1, length(xGrid), length.out=1000))
	indSel <- c(ind1, ind2)
	list(fitDat=fitDat, y=y, indSel=indSel, xGrid=xGrid)
}

#' plot method for nci60Fit class
#' @method plot nci60Fit
#' @param x a nci60Fit object as returned by nci60Fit
#' @param h horizontal line to add indicating e.g. GI (h=0.5), TGI (h=0), LD50 (h=-0.5)
#' @param ... additional parameters, not implemented
#' @export
plot.nci60Fit <- function(fit, xlab="Log10(Dose)", ylab="Relative growth", main='Fit of NCI60 method', 
	xlim=NULL, ylim=NULL, cex.main=1, cex.axis=1, pcol='black', lcol='black', h=c(-0.5, 0, 0.5), lwd=2, ...){
	plL <- getPlotDat_nci60(fit)
	fitDat <- plL$fitDat
	xGrid <- plL$xGrid
	indSel <- plL$indSel
	y <- plL$y
	#if(is.null(ylim)) ylim <- range(pretty(fitDat$response))
	if(is.null(ylim)) ylim <- c(-1, 1)
	if(any(fitDat$response>1)) ylim[2] <- max(fitDat$response)
	if(is.null(xlim)) xlim <- range(pretty(log10(fitDat$dose)))		
	with(fitDat, plot(log10(dose), response, 
		xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main, cex.main=cex.main, cex.axis=cex.axis))
	#browser()
	#xlim1 <- xlim[1]
	#segments(xlim1, 0.5, fit['GI50'], 0.5, lty=2, col='gray')
	#segments(xlim1, 0.0, fit['TGI'], 0.0, lty=2, col='gray')
	#segments(xlim1, -0.5, fit['LC50'], -0.5, lty=2, col='gray')
	lines(log10(xGrid)[indSel], y[indSel], col=lcol, lwd=lwd)

	#browser()
	abline(h=h, col='grey')
} 
#plot.nci60Fit(fit2)



###########################
# optim for Hill
###########################

# objective function: RSS
RSS_Hill <- function(par, d, y, isLogEC50=TRUE){
	Einf <- par[1]; E0 <- par[2]
	logEC50 <- par[3]; HS <- par[4]
	res <- sum((y-fHill(Dose=d, Einf, E0, logEC50, HS))^2)
	if(is.infinite(res) | is.na(res)){
		#browser()
		warning(sprintf('irregular RSS at par: %s\n', str_c(par, collapse=', ')))
		res <- 1e100 # maybe not necessary
	}
	res
}
# param inits a vector of initial values for (Einf, E0, EC50, HS)
hillFit_optim <- function(d, y, inits=NULL){
	dat <- data.frame(dose=d, response=y)
	res <- rep(NA, 8)	# res=c(Einf, E0, EC50, HS)
	names(res) <- c('E0', 'Einf', 'EC50', 'HS',  'Emax', 'MSE', 'Rsq', 'RSE')
	initials <- rep(NA, 4)
	if(is.null(inits)){
		initials[1] <- min(y)
		initials[2] <- max(y)
		initials[3] <- log(median(d))
		initials[4] <- 1.5
	} else {
		if(length(inits)!=4) stop('4 initial values need to be specified!')
		initials <- inits
	}
	constrL <- c(0, 0, -100, 1) # box constraint. HS can be both positive and negative
	constrU <- c(2, 3, 100, 10)
	#optimRes <- try(optim(par=initials, RSS_Hill, y=y, d=d, lower=constrL, 
	#		upper=constrU, control=list(fnscale=1, pgtol=1e-16, factr=1e3, maxit=3000), method="L-BFGS-B"), silent=TRUE) # L-BFGS-B
	#optimRes <- try(optim(par=initials, RSS_Hill, y=y, d=d, control=list(fnscale=1, pgtol=1e-16, factr=1e3, maxit=3000), method="CG"), silent=F) 
	optimRes <- optim(par=initials, RSS_Hill, y=y, d=d, control=list(fnscale=1, pgtol=1e-16, factr=1e3, maxit=3000),lower=constrL,upper=constrU, method="L-BFGS-B")
	# this RSS might be worse than nls; need fine tuning
	# returns in the order (Einf, E0, logEC50, HS)
	#browser()
	if(class(optimRes)!='try-error'){
		#browser()
		temp <- optimRes$par
		res['Einf'] <- temp[1]
		res['E0'] <- temp[2]
		res['EC50'] <- exp(temp[3])
		res['HS'] <- temp[4]
		rss <- optimRes$value # residual sum of squares
		tss <- sum((y-mean(y, na.rm=T))^2, na.rm=T)
		MSE <- rss/length(y)
		Rsq <- 1-rss/tss
		df <- sum(!is.na(y))-4 # four parameters in for Hill model
		rse <- sqrt(rss/df)
		res['MSE'] <- MSE
		res['Rsq'] <- Rsq
		res['RSE'] <- rse
		#res['Emax'] <- predict(fit, list(dose=max(d, na.rm=T))) # response at maximum dose using nls predict method
		res['Emax'] <- fHill(Dose=max(d, na.rm=T), Einf=res['Einf'], E0=res['E0'], logEC50=log(res['EC50']), HS=res['HS'])
	}	
	attr(res, 'class') <- c('numeric', 'hillFit')
	attr(res, 'fitDat') <- dat
	attr(res, 'fit') <- optimRes
	res	
}	
#hillFit_optim(d=mydat$Dose, y=mydat$y, inits=NULL)		
# source('nls_Hill.R')
