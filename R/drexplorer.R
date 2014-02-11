### 0.9.4: add RSE in the model fitting
### 0.9.3: computeIC as Kevin's way. The original approach predicted from the model is still preserved. But default is to use interpolation. 
### name changed in S4 version: flagOutliers4SingleDose ---> NewmanTest
### a function to flag outliers based on table III and IV in Newman's paper. We later change the name to NewmanTest which essentially implements the idea 
### but instead of returning p values, it returns logic values.
### Currently only removes the first outlier. It is possible to sequentially remove the outliers which means n is reduced by 1 in each step as suggested in the paper
### significance level: 0.01 or 0.05
### Input: 
#		  ref: a vector giving the reference values to estimate sigma
#		  obs: a vector giving the observed values where potential outlier might exist
#		  alpha: significance level. Only 0.01 or 0.05 is allowed
#		  recursive: whether to recursively identify outliers. Currently not implemented
### Output: 
#		  indicator: a logic vector specifying if the corresponding point is flagged as outlier
NewmanTest <- function(ref, obs, alpha=0.01, recursive=FALSE){
	check.alpha(alpha)
	f <- sum(!is.na(ref))-1 # f defined in the paper 
	n <- sum(!is.na(obs)) # n defined in the paper
	if(!(f %in% c(5:20, 24, 30, 40, 60)) | !(n %in% c(3:12, 20))) 
		stop("The scenario (sample size combination) is not available for lookup!")
	s <- sd(ref, na.rm=TRUE) # the standard deviation estimated from control group. It is assumed to be the same as the treatment group 
	w <- diff(range(obs, na.rm=TRUE)) # the range statistic from the treatment group
	q <- w/s # the statistic used in D. Newman, Biometrika 1939
	tab <- onePercentTab # default lookup table, table IV in the paper
	if(alpha==0.05) tab <- fivePercentTab
	cutoff <- tab[as.character(f), as.character(n)]
#browser()
	indicator <- rep(FALSE, length(obs))
	if(alpha[1]==1) 
		return(indicator) # at significance level 1, no outlier at all
	if(q>cutoff) {
		temp <- obs-median(obs, na.rm=TRUE)
		if(abs(min(temp, na.rm=TRUE))>abs(max(temp, na.rm=TRUE)))
			ind <- which.min(obs) # minimum is the extreme and excluded. note that NA does not affect the index due to the design of which.min()
		else
			ind <- which.max(obs)
		indicator[ind] <- TRUE
	}
	indicator	
}
#dose <- ryegrass[, 2]
#response <- ryegrass[, 1]
#NewmanTest(ref=response[dose==0], obs=response[dose==3.75], alpha=0.05)

### name changed in S4 version: flagOutliers ---> drOutlier
### flag outliers when multiple dosages exist. This means repeatedly identify outlier at each dose. 
### the usual format for dose-response data: control has dosage=0
### Input: 
#		  drMat: dose-response matrix. the first column being dosage and second column being response. controls are included by specifying dose=0
### Output: 
#		  indicator: a vector with length nrow(drMat) specifying if a data point is outlier. Control points always have status FALSE
drOutlier <- function(drMat, alpha=0.01) {
	check.drMat(drMat)
	dose <- drMat[, 1]
	response <- drMat[, 2]
	ctr <- response[dose==0] # ctr measurements have dose=0
	indicator <- rep(FALSE, length(dose))
	for(d in unique(dose[dose!=0])){
		sel <- which(dose==d)
		trt <- response[sel]
		indicator[sel] <- NewmanTest(ref=ctr, obs=trt, alpha=alpha)
	}
	indicator
}
model.drc <- c('LL.2', 'LL.3', 'LL.3u', 'LL.4', 'LL.5', 'W1.2', 'W1.3', 'W1.4', 'W2.2', 'W2.3', 'W2.4', 'BC.4', 'BC.5', 
			   'LL2.2', 'LL2.3', 'LL2.3u', 'LL2.4', 'LL2.5', 'AR.2', 'AR.3', 'MM.2', 'MM.3')
model.DoseFinding <- c('emax', 'emaxGrad', 'sigEmax', 'sigEmaxGrad', 'exponential', 'exponentialGrad', 'quadratic', 'quadraticGrad', 
			   'betaMod', 'betaModGrad', 'linear', 'linearGrad', 'linlog', 'linlogGrad', 'logistic', 'logisticGrad', 'linInt', 'linIntGrad')
drModels <- function(){
	cat("Models implemented in drc package:\n")
	cat(model.drc, sep='\n')
	cat("Models implemented in DoseFinding package:\n")
	cat(model.DoseFinding, sep='\n')
	cat('More specific information can be found in the original packages.\n')
}	
getPackageName <- function(modelName){
	if(modelName %in% model.drc) {
		return('drc' )
	} else if(modelName %in% model.DoseFinding) {
		return('DoseFinding')
	} else {
	stop('The specified model is not supported. Check supported models by typing drModels()!\n')
	}	
}
### define a class to hold fits from both packages
setOldClass("drc")
setOldClass("DRMod")
setClassUnion("drFit0", list("drc", "DRMod"))
# info: a list element holding the model information, i.e. RSE
setClass('drFit', representation(fit='drFit0', fitDat='matrix', originalDat='matrix', alpha='numeric', fitCtr='logical', tag='character', info='list'))

## added residual standard error (RSE) on 2014/01/06 so that model selection by RSE can be implemented.
drFit <- function(drMat, modelName = "sigEmax", alpha=0.01, fitCtr=FALSE){
	package <- getPackageName(modelName) # which package source is the model implemented
	indicator <- drOutlier(drMat=drMat, alpha=alpha) # for outlier removal
	dose <- drMat[, 1]
	response <- drMat[, 2]
	ctr <- response[dose==0] # ctr measurements have dose=0
	trt <- response[dose!=0] 
	# scale the control and treatment effects with control mean
	ctrScaled <- ctr/mean(ctr)
	trtScaled <- trt/mean(ctr)
#browser()
	dose1 <- dose[dose!=0] # nonzero dosage
	## prepare input matrix: the response is scaled by controls
	if(fitCtr) { 
		dat <- data.frame(dose=dose[!indicator], response=(response/mean(ctr))[!indicator])
	} else {
		tempSel <- (!indicator)[dose!=0]
		dat <- data.frame(dose=dose1[tempSel], response=trtScaled[tempSel])
	}
	rse <- NA
	if(package=='drc'){
		f <- get(modelName) ## i.e. the LL.3 or LL.3u function from drc package
        model <- try(drm(response ~ dose, data=dat, fct=f())) ## pass f to drm in the drc package to fit dose-response curve
        if(inherits(model, "try-error")) {
			stop(cat("Model", modelName, "failed\n"))
        } else {
			rse <- summary(model)$rseMat[1]
		}
	} else {
		model <- try(fitMod('dose', 'response', data=dat, model=modelName))
        if(inherits(model, "try-error")) {
			stop(cat("Model", modelName, "failed\n"))
        } else {
			rse <- sqrt(model$RSS/model$df)
		}
	}
#browser()
	res <- new('drFit', fit=model, fitDat=data.matrix(dat), originalDat=data.matrix(drMat), alpha=alpha, fitCtr=fitCtr, tag=package, info=list(RSE=rse))
	res
}	
	 
setMethod('predict', signature(object='drFit'),
          function(object, newData) {
	if(missing(newData)) newData <- object@originalDat[, 1] 	
	fitObj <- object@fit
	tag <- object@tag
	y <- rep(NA, length(newData))
	if(tag=="drc"){ # using model fitted from drc
		y <- predict(fitObj, newdata=data.frame(dose=newData), od=TRUE)
	}
	if(tag=="DoseFinding"){ # using model fitted from DoseFinding
		modelName <- attributes(fitObj)$model
		f <- get(modelName, pos=which(search() == "package:DoseFinding")) # function name, i.e. sigEmax
		coy <- as.list(fitObj$coefs)
       coy$dose <- newData
		y <- do.call(f, coy)
	}
	y
})


### find IC values by univariate rootfinding from fitted curve. options in uniroot can be specified
RootFindingIC <- function(drFit, percent=0.5, log.d=TRUE, lower, upper, ...) {
	tag <- drFit@tag
	res <- NA
	dose <- drFit@originalDat[, 1]
	if(missing(lower)) lower <- min(c(0, min(dose)))
	if(missing(upper)) upper <- max(dose)*1e60
	if(tag=="drc"){
		objFct <- drFit@fit$fct
		pm <- drFit@fit$parmMat # par mat
		objFct$"fct"(dose, t(pm))
		# define root finding function
		f.drc <- function(d, parm, IC) objFct$"fct"(d, parm)-IC
		res <- uniroot(f.drc, parm=t(pm), IC=percent, lower=lower, upper=upper, ...)$root
	}
	if(tag=="DoseFinding"){
		modelName <- attributes(drFit@fit)$model
		fname <- get(modelName, pos=which(search() == "package:DoseFinding")) # function name, i.e. sigEmax
		coy <- as.list(drFit@fit$coefs) # coefs, a list for do.call
		# define root finding function
		f.DoseFinding <- function(d, coef, IC) {
			coef$dose <- d
			do.call(fname, coef)-IC
		}
		#myf(3.75, coef=coy, IC=0.5)
		res <- uniroot(f.DoseFinding, coef=coy, IC=percent, lower=lower, upper=upper, ...)$root
	}
	if(log.d) res <- log10(res)
	names(res) <- paste('IC', percent*100, sep='')
	res
}

## IC value is approximated by the nearest point in the grid
### it first fits a dose response model from drc or DoseFinding package and then searches for the nearest point from the predictions.
### a model is specified by a modelName and package name.
### models from DoseFinding: modelName="linlog", "linear", "quadratic", "emax", "exponential", "sigEmax", "betaMod" and "logistic"
### models from drc: a lot of models. see getMeanFunctions()
### percent: the inhibition ratio to be searched against. A value between 0 and 1. When percent is beyond the predicted
###			 range, it will be truncated at at either lowest dose or highest dose observed
### log.d: whether to return log10(dose) or the raw dose. Default is set to TRUE.
### nBin:  number of bins to approximate IC values.
### interpolation: whether to use interpolation to estimate IC50. 
### niter: number of equally spaced intervals during interpolation. Only used when interpolation=TRUE.

computeIC <- function(drFit, percent=0.50, log.d=TRUE, interpolation=TRUE, niter=1001, lower, upper, ...) {
	if(interpolation==FALSE){
	if(length(percent)==1) return(RootFindingIC(drFit, percent, log.d, lower=lower, upper=upper, ...))
	if(length(percent)>1) {
		res <- rep(NA, length(percent))
		for(i in 1:length(percent)) {
			tm <- try(RootFindingIC(drFit, percent[i], log.d, lower=lower, upper=upper, ...), silent=TRUE)
			if(class(tm)!='try-error')
			res[i] <- tm
		}
		names(res) <- paste('IC', percent*100, sep='')
		return(res)
	}
	} else {
	drMat <- drFit@originalDat
	dose <- drMat[, 1]
	response <- drMat[, 2]
	ctr <- response[dose==0] # ctr measurements have dose=0
	trt <- response[dose!=0] 
	# scale the control and treatment effects with control mean
	ctrScaled <- ctr/mean(ctr)
	trtScaled <- trt/mean(ctr)
	#browser()
	dose1 <- dose[dose!=0] # nonzero dosage
	# define xlim based on responses treated with drug
	top <- 10^ceiling(log10(max(dose1))) ######### modified: sometimes only 6 doses are observed. This modification guarantees all 7 doses are present
    bot <- 10^floor(log10(min(dose1)))
    xGrid <- seq(bot, top, length=niter) ## set the grid for x-axis, at log10 scale. A lucky choice of niter would give xGrid with only 1 or 2 digits.
#browser()
	yv <- predict(drFit, newData=xGrid) ## predicted values. dose at the original scale
	if(length(percent)==1) {
		res <- icByInterpolation(percent, xv=xGrid, yv, max(dose1), min(dose1))
	} else {
		res <- rep(NA, length(percent))
		for(i in 1:length(percent)) {
			tm <- try(icByInterpolation(percent[i], xv=xGrid, yv, max(dose1), min(dose1)), silent=TRUE)
			if(class(tm)!='try-error')
			res[i] <- tm
		}
		names(res) <- paste('IC', percent*100, sep='')
	}
	if(log.d) {
	return (log10(res))
	} else {
	return (res)
	}
	}
}
# computeIC(fit_sigEmax_alpha_o5, percent=0.5, log.d=TRUE, niter=500)

### IC by interpolation of predicted responses
icByInterpolation <- function(perc, xv, yv, ubound=max(xv), lbound=min(xv)){
	av <- abs(yv-perc)
	wv <- which(av==min(av))[1]
	ic <- xv[wv]
	if (min(yv) > perc) ic <- ubound ## saturated case
	if (ic < lbound) ic <- lbound
	ic
}

setMethod('plot', signature(x='drFit'),
          function(x, pchs=c(16, 17, 15), cols=c(1, 2, 3), col=4, lwd=2, addLegend=TRUE, xlab="Log10(Dose)", ylab="Scaled Response", ylim, xlim, main, ...) {
  	if(missing(main)) main <- attributes(x@fit)$model
	if(missing(ylim)) ylim <- c(0,1.2)
	drMat <- x@originalDat
	## plot the original data points
	# (1) outlier at both levels: for graphical purpose. the actual outlier identification is embedded in model fitting. 
	indicator1 <- drOutlier(drMat=drMat, alpha=0.05)  
	indicator2 <- drOutlier(drMat=drMat, alpha=0.01)
	pCols <- rep(cols[1], nrow(drMat)) # cols[1] for regular points
	pCols[indicator1] <- cols[2] # cols[2] for outliers at 5% significance level
	pCols[indicator2] <- cols[3] # cols[3] for outliers at 1% significance level
	pPchs <- rep(pchs[1], nrow(drMat))
	pPchs[indicator1] <- pchs[2]
	pPchs[indicator2] <- pchs[3]
	# (2) data points and coordinates
	dose <- drMat[, 1]
	response <- drMat[, 2]
	ctr <- response[dose==0] # ctr measurements have dose=0
	trt <- response[dose!=0] 
	# scale the control and treatment effects with control mean
	ctrScaled <- ctr/mean(ctr)
	trtScaled <- trt/mean(ctr)
#browser()
	dose1 <- dose[dose!=0] # nonzero dosage
	# define xlim based on responses treated with drug
	top <- 10^ceiling(log10(max(dose1))) ######### modified: sometimes only 6 doses are observed. This modification guarantees all 7 doses are present
    bot <- 10^floor(log10(min(dose1)))
    xGrid <- seq(bot, top, length=500) ## set the grid for x-axis, at log10 scale
    # initialize the plot: the scaled values for ctr and trt
    par(bg="white")
	## the actual data points
	if(missing(xlim)) xlim <- log10(c(bot, top))
	plot(log10(dose1), trtScaled, ylim=ylim, xlim=xlim,
           xlab=xlab, ylab=ylab, main=main, col=pCols[dose!=0], pch=pPchs[dose!=0])
	points(rep(log10(bot), length(ctrScaled)), ctrScaled, pch=8)	
	abline(h=c(0.5, 0.75), col='grey')
	### add curve from the fitted model
	y <- predict(x, newData=xGrid)
#browser()
	lines(log10(xGrid), y, col=col, lwd=lwd)
	## add legend
	#if(addLegend) legend("topright", c("okay", "95%", "99%"), col=cols, pch=pchs) ## use significance level to remove confusion (modified on 09/03/2013)
	if(addLegend) legend("topright", c("okay", "5%", "1%"), col=cols, pch=pchs)
})
	   
setMethod('lines', signature(x='drFit'),
    function(x, col=5, lwd=2) {
	dose <- x@originalDat[, 1]
	dose1 <- dose[dose!=0] # nonzero dosage
	top <- 10^ceiling(log10(max(dose1))) ######### modified: sometimes only 6 doses are observed. This modification guarantees all 7 doses are present
    bot <- 10^floor(log10(min(dose1)))
    xGrid <- seq(bot, top, length=500) ## set the grid for x-axis, at log10 scale
	y <- predict(x, newData=xGrid)
	lines(log10(xGrid), y, col=col, lwd=lwd)	  
})		  
######### check utility
check.drMat <- function(drMat){
	#browser()
	if(missing(drMat)) stop("Input drMat missing!") 
	if(!is.matrix(drMat) & !is.data.frame(drMat)) stop('Input drMat must be a matrix or a data frame!')
	if(ncol(drMat)>2) warning('More than 2 columns are present in drMat. Only the first 2 columns will be used!')
	if(nrow(drMat)<=3) stop('At least 3 rows are needed for drMat!')
	if(!all(apply(drMat, 2, is.numeric))) stop('All elements in drMat must be numeric!')
}
check.alpha <- function(alpha) {
	if(length(alpha)>1) stop('More than 1 alpha levels specified. Only 1 value can be specified!')
	if(!(alpha %in% c(0.01, 0.05, 1))) stop('Alpha can only be 0.01, 0.05 or 1!')
}
