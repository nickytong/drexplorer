#'  Table III and IV published by Newman
#'
#' This dataset contains two 20-by-12 matrices corresponding to table III and table IV published by Newman, D. In
#' particular, onePercentTab gives the quantiles at significance level of 0.01 and fivePercentTab gives the quantiles at significance level of 0.01.
#' The degree of freedom f and sample size n are used as row and column names for each table.
#'
#'
#' @format Each of the two objects (onePercentTab, fivePercentTab)
#' is a matrix with 20 rows (degree of freedom f) and 12 columns (sample size n).
#' @references Newman, D. (1939). The distribution of range in samples from a normal population, 
#' expressed in terms of an independent estimate of standard deviation. Biometrika, 31(1/2), 20-30.
#' @name NewmanTables
#' @aliases onePercentTab,fivePercentTab
#' @seealso \code{\link{NewmanTest}}
NULL

#'  ryegrass dataset from drc package
#'
#' This dataset is from the ryegrass data in drc package. It is a 2 column data frame where the conc column is the concentration of ferulic acid in mM
#' and rootl column is the root length of ryegrass. 
#'
#'
#' @references Inderjit and J. C. Streibig, and M. Olofsdotter (2002) Joint action of phenolic acid mixtures and its
#' significance in allelopathy research, Physiologia Plantarum, 114, 422–428, 2002
#' @name ryegrass
NULL



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

#' Identifying outliers with the range to standard deviation ratio statistic.
#' 
#' This function implements the method described 
#'  by D. Newman, Biometrika, 1939 to identify outliers.
#' 
#' Given measurements from controls (no drug treated), we can compute the sample standard deviation (s). The range of responses from
#' treated samples (w) can be computed for a given dose level. Assuming the controls to have the same variation as the drug treated case,
#' the distribution of ratio statistic q=w/s can be derived and used to calculate if there is outliers in the treated responses as described 
#' by D. Newman, Biometrika, 1939.  
#' 
#' Note that this function works for a single dose level. When multiple dose levels exist, one need to repeatedly call this function to identify
#' outliers at each dose level or use the flagOutliers() function which is just a wrapper.
#'
#' @param ref a vector giving the reference values to estimate sigma
#' @param obs a vector giving the observed values where potential outlier might exist
#' @param alpha significance level. Only 0.01, 0.05 or 1 is allowed
#' @param recursive whether to recursively identify outliers. Currently not implemented
#' @return indicator a logical vector specifying if the corresponding point is flagged as outlier
#' @export
#' @references Newman, D. (1939). The distribution of range in samples from a normal population, 
#' expressed in terms of an independent estimate of standard deviation. Biometrika, 31(1/2), 20-30.\url{http://www.jstor.org/stable/2334973}
#' @seealso \code{\link{drOutlier}, \link{drModels}, \link{drFit}, \link{drFit-class}}
#' @examples
#' set.seed(1)
#' x <- rnorm(10, 0, 1)
#' y <- c(rnorm(5, 0, 1), rnorm(1, 0, 1)+4)
#' # the last observation in y is an outlier
#' NewmanTest(x, y, alpha=0.05)
#'
NewmanTest <- function(ref, obs, alpha=0.01, recursive=FALSE){
	check.alpha(alpha)
	indicator <- rep(FALSE, length(obs))
	if(alpha[1]==1) 
		return(indicator) # at significance level 1, no outlier at all
	f <- sum(!is.na(ref))-1 # f defined in the paper 
	n <- sum(!is.na(obs)) # n defined in the paper
	#if(!(f %in% c(5:20, 24, 30, 40, 60)) | !(n %in% c(3:12, 20))) 
	#	stop("The scenario (sample size combination) is not available for lookup!")
	# modified on 2014/02/21: when there is no test result, just return all FALSE and give an warning
	if(!(f %in% c(5:20, 24, 30, 40, 60)) | !(n %in% c(3:12, 20))) {
		warning("The scenario (sample size combination) is not available for lookup!\n Treat all measurements as non-outlier")
		return(indicator)
	}
	s <- sd(ref, na.rm=TRUE) # the standard deviation estimated from control group. It is assumed to be the same as the treatment group 
	w <- diff(range(obs, na.rm=TRUE)) # the range statistic from the treatment group
	q <- w/s # the statistic used in D. Newman, Biometrika 1939
	tab <- onePercentTab # default lookup table, table IV in the paper
	if(alpha==0.05) tab <- fivePercentTab
	cutoff <- tab[as.character(f), as.character(n)]
#browser()
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
#' Identifying outliers for dose-response data
#'
#' A wrapper for the NewmanTest() function. It repeatedly call NewmanTest() for each dose level. 
#'
#' @param drMat dose-response matrix. the first column being dosage and second column being response. 
#'   controls are included by specifying dose=0
#' @param alpha a scalar for significance level. Only 0.01, 0.05 and 1 are allowed. 
#' 	alpha=1 is included for comparability issue. In this case,  no outliers will be identified. 
#' @return indicator a vector with length nrow(drMat) specifying if a data point is outlier. 
#'  Control points always have status FALSE
#' @export
#' @author Kevin R Coombes (\email{kcoombes@@mdanderson.org}), Pan Tong (\email{nickytong@@gmail.com})
#' @seealso \code{\link{NewmanTest}, \link{drModels}, \link{drFit}, \link{drFit-class}}
#' @examples
#' data(ryegrass, package='drc') # use the ryegrass data from drc package
#' drOutlier(drMat=ryegrass[, c(2, 1)], alpha=0.05)
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
### these models are before 0.9.7 where not all are feasible
model.drc_old <- c('LL.2', 'LL.3', 'LL.3u', 'LL.4', 'LL.5', 'W1.2', 'W1.3', 'W1.4', 'W2.2', 'W2.3', 'W2.4', 'BC.4', 'BC.5', 
			   'LL2.2', 'LL2.3', 'LL2.3u', 'LL2.4', 'LL2.5', 'AR.2', 'AR.3', 'MM.2', 'MM.3')
model.DoseFinding_old <- c('emax', 'emaxGrad', 'sigEmax', 'sigEmaxGrad', 'exponential', 'exponentialGrad', 'quadratic', 'quadraticGrad', 
			   'betaMod', 'betaModGrad', 'linear', 'linearGrad', 'linlog', 'linlogGrad', 'logistic', 'logisticGrad', 'linInt', 'linIntGrad')
### after 0.9.7, only below models are allowed			   
model.drc <- c('LL.2', 'LL.3', 'LL.3u', 'LL.4', 'LL.5', 'W1.2', 'W1.3', 'W1.4', 'W2.2', 'W2.3', 'W2.4', 'BC.4', 'BC.5', 
			   'LL2.2', 'LL2.3', 'LL2.3u', 'LL2.4', 'LL2.5', 'AR.2', 'AR.3')
model.DoseFinding <- c('emax', 'sigEmax', 'exponential', 'quadratic',
			   'linear', 'linlog', 'logistic')
### frequently used models for drug screening

#' frequently used models
#'
#' this function returns frequently used dose-response models for anti-cancer drug screening.
#'
#' @export
recommendedModels <- function(){
	c('sigEmax', 'LL.4', 'LL.5', 'linear', 'logistic') # LL.3 might be better than linlog; linlog removed since it is not appropriate, especially predict always -1.37 due to the default off parameter in the original DoseFinding package
}

#' Show available dose-response models with direct support. 
#' 
#' This function shows available models to be passed to drFit function as modelName.
#'
#' @param return whether to return the models in a list
#' @param verbose whether to print out the models
#' @return Following is the print out:
#'   Models implemented in drc package:
#'   LL.2
#'   LL.3
#'   LL.3u
#'   LL.4
#'   LL.5
#'   W1.2
#'   W1.3
#'   W1.4
#'   W2.2
#'   W2.3
#'   W2.4
#'   BC.4
#'   BC.5
#'   LL2.2
#'   LL2.3
#'   LL2.3u
#'   LL2.4
#'   LL2.5
#'   AR.2
#'   AR.3
#'   MM.2
#'   MM.3
#'   Models implemented in DoseFinding package:
#'   emax
#'   emaxGrad
#'   sigEmax
#'   sigEmaxGrad
#'   exponential
#'   exponentialGrad
#'   quadratic
#'   quadraticGrad
#'   betaMod
#'   betaModGrad
#'   linear
#'   linearGrad
#'   linlog
#'   linlogGrad
#'   logistic
#'   logisticGrad
#'   linInt
#'   linIntGrad
#'   More specific information can be found in the original packages.
#' 
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drFit}, \link{drFit-class}}
#' @examples
#' drModels()
#' @export
drModels <- function(return=FALSE, verbose=TRUE){
	if(verbose) {
	cat("Models implemented in drc package:\n")
	cat(model.drc, sep='\n')
	cat("Models implemented in DoseFinding package:\n")
	cat(model.DoseFinding, sep='\n')
	cat('More specific information can be found in the original packages.\n')
	cat("Models frequently used in drug screening experiments:\n")
	cat(recommendedModels(), sep='\n')
	}
	if(return){
		return(list(drc=model.drc, doseFinding=model.DoseFinding, recommendedModels=recommendedModels()))
	}
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

## this is wrong
## plot method for drFit class.
## @name plot
## @rdname drFit-class
setGeneric('plot', package='graphics' )


## lines method for drFit class.
## @name lines
## @rdname drFit-class
setGeneric('lines', package='graphics' )

## predict method for drFit class.
## @name predict
## @rdname drFit-class
setGeneric('predict', package='stats' )

# info: a list element holding the model information, i.e. RSE
#' Class \code{"drFit"}
#'
#' @name drFit-class
#' @rdname drFit-class
#' @exportClass drFit
#' @section Creating Objects: 
#' Although objects can be created directly using \code{new("drFit", ...)}, the most common usage will be to pass dose-response data to the \code{drFit} function.
#'
#' @slot fit A fit object by either drc or DoseFinding package. 
#' @slot fitDat The actual data matrix (after excluding outliers) used to fit a dose-response model.
#' @slot originalDat The original dose response data, a matrix.
#' @slot alpha A scalar for significance level. This specifies the significance level to identify outliers
#' @slot fitCtr A logical value specifying whether to include the control points into the model fitting. 
#' @slot tag A string (either 'drc' or 'DoseFinding') tracking which package is used for model fitting. 
#' @slot info A list that holds information related to the model, i.e. Residual Standard Error (rse).  
#' @seealso \code{\link{drOutlier}, \link{drModels}, \link{drFit}, \link{drFit-class}}
setClass('drFit', representation(fit='drFit0', fitDat='matrix', originalDat='matrix', alpha='numeric', fitCtr='logical', tag='character', info='list'))

noNA <- function (dat) 
{
    sel <- complete.cases(dat)
    if (is.null(dim(dat))) {
        res <- dat[sel]
    }
    else {
        res <- dat[sel, ]
    }
    res
}
#' prepare dose-reponse data
#'
#' this function scaled the response by mean reponse in control when necessary
#'
#' the standardization means response (e.g. count) is to be scaled by mean control response so that
#' the standardized response is relative viability, usually between 0 and 1.
#' this function first detects outlier data points (using supplied data, either scaled or not scaled). It
#' then scale the data and split the data into two data frames: one for dose not 0 and one for dose at 0 
#' Ideally, outlier detection is better at original count (or signal intensity) due to asymptotic assumption
#' But for data already scaled, detection outlier at this level is what we can do
#'
#' @param drMat dose-response matrix as for drFit. the first column being dosage and second column being response. controls are included by specifying dose=0
#' @param alpha a scalar for significance level. This specifies the significance level to identify outliers which will be excluded from model fitting. To include
#'  all data, set alpha=1. 
#' @param fitCtr A logic vector specifying whether to include the control points into the model fitting.
#' @param standardize whether to standardize (scale) the data based on control points. This should be disabled when no control data is supplied
#' @return a list
prepDRdat <- function(drMat, alpha=0.01, fitCtr=FALSE, standardize=TRUE){
	# if drMat itself is scaled data, then outlier remover is on scaled data
	indicator <- drOutlier(drMat=drMat, alpha=alpha) # for outlier removal
	indicator1 <- suppressWarnings(drOutlier(drMat=drMat, alpha=0.05)) # for plot purpose indicating different levels of outlier status
	indicator2 <- suppressWarnings(drOutlier(drMat=drMat, alpha=0.01)) # for plot purpose indicating different levels of outlier status
	drMat <- noNA(drMat)
	hasCtrl <- any(drMat[, 1]==0) # logical indicating if there is control data by testing if any dose is 0
	dose <- drMat[, 1]
	response <- drMat[, 2]
	#indCtrl <- round(dose, 5)==0 # dose might be 1e-8 and thus it is not control!
	indCtrl <- dose==0
	indTrt <- !indCtrl
	ctr <- response[indCtrl] # ctr measurements have dose=0
	trt <- response[indTrt] 		
	if(standardize){
		# scale the control and treatment effects with control mean
		## datAll is scalled data including control and outlier
		datAll <- data.frame(dose=dose, response=(response/mean(ctr)))
	} else {
		# no standardization
		datAll <- data.frame(dose=dose, response=response)
	}
	# dat is the data used for fitting models
	if(fitCtr) { 
		dat <- datAll[!indicator, ] # just remove outlier and include control
	} else {
		tempSel <- !indicator & indTrt
		dat <- datAll[tempSel, ]
	}
	isCtrl <- indCtrl # this include outliers; 
	isTrt <- indTrt # this include outliers; 
	#browser()
	list(dat=dat, datAll=datAll, isCtrl=isCtrl, isTrt=isTrt, isOutlier=indicator, indicator1=indicator1, indicator2=indicator2, 
		standardize=standardize, fitCtr=fitCtr, alpha=alpha, hasCtrl=hasCtrl)
}
## added residual standard error (RSE) on 2014/01/06 so that model selection by RSE can be implemented.
## notice: ******** drFit internally scales the response and fit the models on the scaled values. It is prohibited to scale the data before feeding to drexplorer!!!! *********
#' Fit a dose-response model
#' 
#' This function can fit various dose-response models by specifying a model name and package source (either drc or DoseFinding).
#'
#' This is a wrapper to fit dose response models from drc and DoseFinding packages. 
#' When fit the model, the response is internally scaled with mean response in control. Therefore, the fitted response is a ratio.
#' In visualization, i.e. dose-response curve plotting, the dose is log10 transformed so that we can see the sigmoid curve.
#' However, in computing IC50 and doing prediction, the dose is in original scale since the model is trained in the original scale of dose. Correspondingly, predicted
#' response is also scaled response.
#' When the user requires a log10 transformed dose, it is done after estimating dose in original scale through the computeIC() function.
#'
#' @param drMat dose-response matrix. the first column being dosage and second column being response. controls are included by specifying dose=0
#' @param modelName a dose-response model. For available models to be specified, see drModels().
#' @param alpha a scalar for significance level. This specifies the significance level to identify outliers which will be excluded from model fitting. To include
#'  all data, set alpha=1. 
#' @param fitCtr A logic vector specifying whether to include the control points into the model fitting.
#' @param standardize whether to standardize (scale) the data based on control points. This should be disabled when no control data is supplied
#' @return the function returns a drFit S4 object.
#' @export
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit-class}}
#' @examples
#' data(ryegrass, package='drc') # use the ryegrass data from drc package
#' # fit a sigmaEmax model without outlier removal. the controls are excluded from model fitting
#' fit.LL.3 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3", alpha=0.01, fitCtr=FALSE)
#' fit.LL.3u <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3u", alpha=0.01, fitCtr=FALSE)
#' fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.01, fitCtr=FALSE)
#' @author Kevin R Coombes (\email{kcoombes@@mdanderson.org}), Pan Tong (\email{nickytong@@gmail.com})
drFit <- function(drMat, modelName = "sigEmax", alpha=0.01, fitCtr=FALSE, standardize=TRUE){
	package <- getPackageName(modelName) # which package source is the model implemented
	#browser()
	resPrepDR <- prepDRdat(drMat, alpha=alpha, fitCtr=fitCtr, standardize=standardize)
	dat <- resPrepDR$dat
	rse <- NA
	if(package=='drc'){
		f <- get(modelName) ## i.e. the LL.3 or LL.3u function from drc package
        model <- try(drc::drm(response ~ dose, data=dat, fct=f())) ## pass f to drm in the drc package to fit dose-response curve
        if(inherits(model, "try-error")) {
			stop(cat("Model", modelName, "failed\n"))
        } else {
			rse <- summary(model)$rseMat[1]
		}
	} else {
		model <- try(DoseFinding::fitMod('dose', 'response', data=dat, model=modelName))
        if(inherits(model, "try-error")) {
			stop(cat("Model", modelName, "failed\n"))
        } else {
			rse <- sqrt(model$RSS/model$df) # sqrt(RSS/df)
		}
	}
#browser()
	res <- new('drFit', fit=model, fitDat=data.matrix(dat), originalDat=data.matrix(drMat), alpha=alpha, fitCtr=fitCtr, tag=package, 
			info=list(RSE=rse, standardize=standardize, resPrepDR=resPrepDR))
	res
}
#fit <- drFit(drMat=mydat[, c('Dose', 'y')], standardize=F, modelName = "sigEmax", alpha=1, fitCtr=FALSE)	
###
### when predicting, what dose is used, the original or the scaled dose?????
###	there is no scaling on dose, only scalin on value, only transformed. So this is a wrong question; but the answer is, the original dose is used!!!


#' Calculate AUC from a fitted object
#'
#' AUC is calculated through the integrate() function based on dose-response curve.
#'
#' @param fit usually a drFit object. However, any object with a predict method would work.
#' @param dmin minimum dose. The integral range is [dmin, dmax].
#' @param dmax maximum dose. The integral range is [dmin, dmax].
#' @param islogd whether the supplied dose dmin/dmax is in log10 scale. The user should be responsible for the consistency between the actual value
#' of dmin, dmax and islogd. If log10 transformed dmin/dmax is supplied with islogd=TRUE, the AUC is calculated based on dose-response curve with x-axis being log10(dose); On the
#' other hand, if dmin, dmax is original scale and islogd=FALSE, the AUC is calculated based on a dose-response curve with x-axis being dose
#' @return a vector of AUC, AUC0 and AUCs. AUC is the area under dose response curve; AUC0 is the area under the line response=1; AUCS is AUC/AUC0
#' @export
computeAUC <- function(fit, dmin, dmax, islogd=TRUE) {
	# assumption: the predict function is trained on original dose
	# the rationale: predict response at each dose; apply integration;
	# AUC can be calculated based on log10 dr curve or original dr curve (in this case, absolute AUC is smaller, scaled AUC is also smaller usually); Notice that this only changes the integration range/scale (x) without affecting response (y)
	# However, the user should supply dmin/dmax in log10 scale if islogd==TRUE (for log10 AUC) and dmin/dmax in original scale for islogd=FALSE (for untransformed AUC)
	# according to the Nature paper, AUC is defined on the log10 scale; this is also more intuitive since it matches dr curve which is in log10 scale
	# islogd: this makes it possible to calculate AUC either on log10dose or original dose. However, the user
	#      should make sure dmin and dmax is on the same scale (take log10 correspondingly)
	f_response <- function(Dose, islogd=islogd) {
		if(islogd) {
			dd <- 10^(Dose)
		} else {
			dd <- Dose
		}
		## currently only implements for drFit
		#if(class(fit)=='drFit'){
		#	res <- predict(fit, newData=dd)
		#}
		res <- predict(fit, newData=dd)
		res  
	}
	# reference response, e.g. the line response=1; this can be used to scale AUC
	f_responseRef <- function(Dose, reference=1) {
		res <- rep(reference, length(Dose))
		res 
	}
	#browser()
	AUC <- try(integrate(f_response, lower=dmin, upper=dmax, islogd=islogd)$value, silent=TRUE)  
	AUC0 <- try(integrate(f_responseRef, lower=dmin, upper=dmax)$value, silent=TRUE)  
	if(inherits(AUC, 'try-error')) {
		cat(AUC)
		AUC=NA
	}
	if(inherits(AUC0, 'try-error')) {
		cat(AUC0)
		AUC0=NA
	}
	res <- c(AUC=AUC, AUC0=AUC0, AUCs=AUC/AUC0)
	res
}
#computeAUC(fit_sigEmax_alpha_o5, dmin=-0.027, dmax=1.477, islogd=T)
#computeAUC(fit, dmin=-8, dmax=-5.5, islogd=T)

#browser()
#### @rdname drFit-class
#### @aliases predict,drFit-method
#' predict method for drFit object
#'
#' @param object a drFit object
#' @param newData a vector of dose levels to predict at. Default is to predict response at observed dose levels.
#' @return a vector of predicted values
#' @aliases predict,drFit-method
#' @export 
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
		#browser()
		# linlog always return y=-1.37; problematic for this model!
		modelName <- attributes(fitObj)$model
		#f <- get(modelName, pos=which(search() == "package:DoseFinding")) # function name, i.e. sigEmax
		func <- get(modelName) # function name, i.e. sigEmax
		coy <- as.list(fitObj$coefs)
        coy$dose <- newData
		y <- do.call(func, coy)
	}
	y
})


### find IC values by univariate rootfinding from fitted curve. options in uniroot can be specified
# use dmin to indicate non-achievable IC due to to small IC (drug too effective based on the curve): this is original scale, so minimum dose is 0
# use dmax to indicate non-achievable IC due to to high IC (drug too ineffective based on the curve)
# 2014/11/28: to make the function more general, percent means desired response from fitted model, not IC. Thus, percent
#   essentially can be negative. computeIC function needs to convert between IC and percent (response) by: 1-x(IC)/100=percent
# percent is actually desired response here
#
# Sometimes, the curve is increasing (e.g. due to experimental noise or unknown reason). We need rule to override the model based estimate for response out of the range
# Notice this rule only deals with cases where required response is out of range. The emirical rule solved flat curve when cell line is extremely sensitive (e.g. y=0.1). 
# If the curve is increasing from 0.1 to 0.9, such rule might be wrong
# This rule also applies to interpolation case
# Empirical rule:
#	(1) myIC is out of range, then assign dmin or dmax based on if myIC is larger than largest or smaller than smallest
#	(2) myIC is within range, then root finding
RootFindingIC <- function(drFit, percent=0.5, log.d=TRUE, lower, upper, dmin=0, dmax=Inf, ...) {
	#### the predicted response is bounded with theoretical lower and upper bound. When the required response specified by IC is out of the theoretical range, we need to specify the returned dose
	#### originally we specify NA, in which case we may get IC50=NA, which can be due to too sensitive or too negative. Thus, it is better to give a value.
	#### we assign dmin for the scaled response larger than theoretical response; we assign dmax for response smaller than theretical minimum response.
	#### notice this only affects the IC by prediction result
	#dmin <- 1e-30 # rediculously low dose
	#dmax <- 1e30 # rediculously high dose
	#dmin <- -Inf # rediculously low dose
	#dmax <- Inf
	tag <- drFit@tag
	res <- NA
	dose <- drFit@originalDat[, 1]
	#if(missing(lower)) lower <- min(c(1e-30, min(dose)))
	if(missing(lower)) lower <- 1e-60 # cannot allow log(lower)=log(0) as uniroot lower bound
	#if(missing(upper)) upper <- max(dose)*1e60
	# modified on 02/22/2014: upper too large makes f(d) Nan in sigEmax. sigEmax is between e0+eMax ~ e0
	if(missing(upper)) upper <- 1e60 # assume maximum dose will not exceed 1e10
	myIC <- percent #
	if(tag=="drc"){
		objFct <- drFit@fit$fct
		pm <- drFit@fit$parmMat # par mat
		objFct$"fct"(dose, t(pm))
		# define root finding function
		#f.drc <- function(d, parm, IC) objFct$"fct"(d, parm)-IC
		fl.drc <- function(ld, parm, IC) objFct$"fct"(exp(ld), parm)-IC
		#res <- uniroot(f.drc, parm=t(pm), IC=percent, lower=lower, upper=upper, ...)$root
		#browser()
		# max: f.drc(d=lower, parm=t(pm), IC=0)
		# min: f.drc(d=upper, parm=t(pm), IC=0)
		# myIC is the scaled response
		#myIC <- 1-percent
		#ymin <- f.drc(d=upper, parm=t(pm), IC=0)
		#ymax <- f.drc(d=lower, parm=t(pm), IC=0)
		ymin <- fl.drc(ld=log(upper), parm=t(pm), IC=0)
		ymax <- fl.drc(ld=log(lower), parm=t(pm), IC=0)
		# in case the curse is increasing, ymax < ymin
		# we swap them
		if(ymax < ymin) {
			tt <- ymax
			ymax <- ymin
			ymin <- tt
		}
		#browser()
		if(myIC>ymax) {
			res <- dmin # required response > theoretical maximum response, throw NA
		} else if(myIC<ymin){
			res <- dmax # required response < theoretical minimum response, throw NA
		} else {	
			#res <- uniroot(f.drc, parm=t(pm), IC=myIC, lower=lower, upper=upper, ...)$root # to feed with biological IC
			res <- exp(uniroot(fl.drc, parm=t(pm), IC=myIC, lower=log(lower), upper=log(upper), ...)$root) 
		}
	}
	if(tag=="DoseFinding"){
		modelName <- attributes(drFit@fit)$model
		#fname <- get(modelName, pos=which(search() == "package:DoseFinding")) # function name, i.e. sigEmax
		fname <- get(modelName) # function name, i.e. sigEmax
		coy <- as.list(drFit@fit$coefs) # coefs, a list for do.call
		# define root finding function: a function of d. coef specifies 4 essential parameter: e0, emax, ed50, h
		fl.DoseFinding <- function(ld, coef, IC) {
			coef$dose <- exp(ld)
			do.call(fname, coef)-IC
		}
		ymax <- fl.DoseFinding(ld=log(lower), coef=coy, IC=0)
		ymin <- fl.DoseFinding(ld=log(upper), coef=coy, IC=0)
		#browser()
		# in case the curse is increasing, ymax < ymin
		# we swap them
		if(ymax < ymin) {
			tt <- ymax
			ymax <- ymin
			ymin <- tt
		}
		# y at dose=d:  f.DoseFinding(d=0, coef=coy, IC=0)
		#myf(3.75, coef=coy, IC=0.5)
		#res <- uniroot(f.DoseFinding, coef=coy, IC=percent, lower=lower, upper=upper, ...)$root
		#browser()
		if(myIC>ymax) {
			res <- dmin # required response > theoretical maximum response, throw NA
		} else if(myIC<ymin){
			res <- dmax # required response < theoretical minimum response, throw NA
		} else {	
			# now myIC is between ymax and ymin
			res <- exp(uniroot(fl.DoseFinding, coef=coy, IC=myIC, lower=log(lower), upper=log(upper), ...)$root)
		}	
	}
	if(log.d) res <- log10(res)
	names(res) <- paste('response_', percent, sep='')
	#browser()
	res
}

findDoseGivenResponse <- function(drFit, response=0.50, log.d=TRUE, interpolation=TRUE, stepLen=NA, lower, upper, ...) {
	#percent <- 1-percent ### convert to biological percent; updated on 2014/02/13---> not ok: the order is reversed!
	if(interpolation==FALSE){
		if(length(response)==1) return(RootFindingIC(drFit, response, log.d, lower=lower, upper=upper, ...))
		if(length(response)>1) {
			res <- rep(NA, length(response))
			for(i in 1:length(response)) {
				tm <- try(RootFindingIC(drFit, response[i], log.d, lower=lower, upper=upper, ...), silent=TRUE)
				if(class(tm)!='try-error')
				res[i] <- tm
			}
			names(res) <- paste('response_', response, sep='')
			#return(res)
		}
	} else {
		resPrepDR <- drFit@info$resPrepDR # result of prep DR data 
		datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
		isTrt <- resPrepDR$isTrt # global indicator if data is treatment
		dose <- datAll$dose # may include 0 dose
		dose1 <- dose[isTrt]
		#
		gg <- format_grid(dose1=dose1, stepLen=stepLen, resolution=100)
		top <- gg$top
		bot <- gg$bot
		xGrid <- gg$xGrid
		yv <- predict(drFit, newData=xGrid) ## predicted values. dose at the original scale
		if(length(response)==1) {
			res <- icByInterpolation(response, xv=xGrid, yv, max(dose1), min(dose1))
		} else {
			res <- rep(NA, length(response))
			for(i in 1:length(response)) {
				tm <- try(icByInterpolation(response[i], xv=xGrid, yv, max(dose1), min(dose1)), silent=TRUE)
				if(class(tm)!='try-error')
				res[i] <- tm
			}
			names(res) <- paste('response_', response, sep='')
		}
		# final presentation
		if(log.d) {
		res <- log10(res)
		}
	}
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
### stepLen: number of equally spaced intervals during interpolation. Only used when interpolation=TRUE.
### updated on 2014/02/13: percent is actually 1-percent, which means IC70 is actually IC30 in biology; We can solve this simply
###				by replace percent with 1-percent in each calculation
### Interpretation of IC30: this is the concentration needed to inhibit 30% of the cells. That is, this concentration makes the cell
### to have 70% of the numbers when compared to no drug. Easy to show that IC30<IC50<IC70. 
#' Compute the IC values at specified percentiles
#' 
#' This function uses rootfinding with fitted curve to compute IC values.
#'
#' @param drFit A drFit object as returned by drFit() function.
#' @param percent the inhibition ratio to be searched against. A vector between 0 and 1. Corresponding
#'   response is 1-percent
#' @param log.d whether to return log10(dose) or the raw dose. Default is set to TRUE.
#' @param interpolation whether to use interpolation to estimate IC values. In this case, the computed IC values will be bound by the observed dosages.
#' @param stepLen step length to construct equally spaced intervals during interpolation. Only used when interpolation=TRUE.
#' @param lower lower bound in root search. Default is min(c(0, min(dose))) where dose is the observed dose levels from drFit@@originalDat.
#' @param upper upper bound in root search. Default is max(dose)*1e6 where dose is the observed dose levels from drFit@@originalDat.
#' @param ... other optimization parameters to be passed to uniroot().
#' @return A named vector giving the IC values at specified percentiles.
#' @seealso \code{\link{NewmanTest}, \link{drOutlier}, \link{drModels}, \link{drFit}}
#' @author Kevin R Coombes (\email{kcoombes@@mdanderson.org}), Pan Tong (\email{nickytong@@gmail.com})
#' @examples
#' data(ryegrass, package='drc') # use the ryegrass data from drc package
#' fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.01, fitCtr=FALSE)
#' computeIC(fit.sigEmax, percent=seq(0, 1, by=0.1), log.d=FALSE)	
#' @export
computeIC <- function(drFit, percent=0.50, log.d=TRUE, interpolation=TRUE, stepLen=NA, lower, upper, ...) {
	res <- findDoseGivenResponse(drFit, response=1-percent, log.d=log.d, interpolation=interpolation, stepLen=stepLen, lower=lower, upper=upper, ...)
	names(res) <- paste('IC', percent*100, sep='')
	res
}

# computeIC(fit_sigEmax_alpha_o5, percent=0.5, log.d=TRUE, niter=500)

### IC by interpolation of predicted responses
# Sometimes, the curve is increasing (e.g. due to experimental noise or unknown reason). We need rule to override the model based estimate for response out of the range
# Notice this rule only deals with cases where required response is out of range
# Empirical rule:
#	(1) myIC is out of range, then assign dmin or dmax based on if myIC is larger than largest or smaller than smallest
#	(2) myIC is within range, then root finding
icByInterpolation <- function(perc, xv, yv, ubound=max(xv), lbound=min(xv)){
	#browser()
	#perc <- 1- perc # revert to biological IC starting from 0.9.5
	# 2014/11/28: perc is the desired response now
	av <- abs(yv-perc)
	wv <- which(av==min(av))[1]
	yr <- range(yv, na.rm=TRUE)
	if(perc>yr[2]){ # response too large
		ic <- lbound 
	} else if(perc<yr[1]) { # response too small
		ic <- ubound
	} else { # response in range
		ic <- xv[wv]
	}
	#browser()
	## modified on 02/22/2014: this fixed that interpolated IC might be larger than maximum dose
	#if (min(yv) > perc) ic <- ubound ## saturated case
	# truncation
	if (ic>ubound) ic <- ubound ## saturated case
	if (ic < lbound) ic <- lbound
	#browser()
	ic
}
#computeIC(fit_sigEmax_alpha_o5, percent=0.6, log.d=T, interpolation=TRUE)

## a function that computes grid (xGrid) and xlim (top, bot) for plotting
# resolution: mid_d/resolution is the step length
# stepLen: this overrides the default way of computing step length
format_grid <- function(dose1, resolution=50, stepLen=NA){
	top <- 10^ceiling(log10(max(dose1))) ######### modified: sometimes only 6 doses are observed. This modification guarantees all 7 doses are present
	bot <- 10^floor(log10(min(dose1)))
    dif <- diff(sort(dose1, decreasing=FALSE)); 
	min_d <- min(dif[dif!=0]) # add sort so that dif will always be positive
	if(is.na(stepLen)){
		stepLen <- min_d/resolution
	}
	xGrid <- seq(bot, top, by=stepLen) ## set the grid for x-axis, at log10 scale
	# modified on 02/22/2014: the grid might be too coarse that misses some observed dose; we need to force the observed dose in
	xGrid <- sort(unique(c(xGrid, dose1)), decreasing=FALSE)
	list(top=top, bot=bot, xGrid=xGrid, stepLen=stepLen)
}


#' plot method for drFit object
#' 
#' @param x a drFit object
#' @param pchs pchs (a vector) specify symbols to show different data points. In particular,
#'	pchs[1] is used to show regular points. pchs[2] is used to show outliers at 0.05 significance level and pchs[3] is used 
#'	to show outliers at 0.01 significance level.  
#' @param cols Similarly, cols[1] is used to color regular points, cols[2] to color outliers at 0.05 significance level and
#'	cols[3] to color outliers at 0.01 significance level. 
#' @param col color used to specify the appearance of fitted curve
#' @param lwd line width used to specify the appearance of fitted curve
#' @param addLegend whether to add legend indicating outlier status
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ylim y limit for display
#' @param xlim x limit for display
#' @param main main title 
#' @param style when style='full', observed dose-response pairs are plotted including controls, fitted curve is superimposed and legend for outliers indicated;
#'  when style='simple', only fitted curve is plotted (without dose-response points and of course, not outlier status). 
#' @param bty bty passed to legend
#' @param h horizontal line to add indicating e.g. IC50 (h=0.5)
#' @param cex.main cex for main title
#' @param cex.axis cex for axis annotation
#' @param cex.lab cex for axis label
#' @aliases plot,drFit-method
#' @export 
setMethod('plot', signature(x='drFit'),
          function(x, pchs=c(16, 17, 15), cols=c(1, 2, 3), col=4, lwd=2, addLegend=TRUE, xlab="Log10(Dose)", 
		  ylab="Relative viability", ylim=NA, xlim=NA, main, style='full', bty='n', h=c(0.5), cex.main=1, cex.axis=1, cex.lab=1) {
  	if(missing(main)) main <- attributes(x@fit)$model
	#browser()
	#tt <- prepDRdat(drMat, alpha=1, fitCtr=FALSE, standardize=x@info$standardize)
	resPrepDR <- x@info$resPrepDR # result of prep DR data 
	datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
	isCtrl <- resPrepDR$isCtrl # global indicator if data is control
	isTrt <- resPrepDR$isTrt # global indicator if data is treatment
	hasCtrl <- resPrepDR$hasCtrl
	indicator1 <- resPrepDR$indicator1
	indicator2 <- resPrepDR$indicator2
	dose <- datAll$dose # may include 0 dose
	dose1 <- dose[isTrt]
	trtScaled <- datAll$response[isTrt] # response at nonzero dosage
	ctrScaled <- datAll$response[isCtrl] # response at nonzero dosage
	#browser()
	if(is.na(ylim[1])) ylim <- range(pretty(datAll$response))
	## plot the original data points
	# (1) outlier at both levels: for graphical purpose. the actual outlier identification is embedded in model fitting. 
	#indicator1 <- drOutlier(drMat=drMat, alpha=0.05)  
	#indicator2 <- drOutlier(drMat=drMat, alpha=0.01)
	pCols <- rep(cols[1], nrow(datAll)) # cols[1] for regular points
	pCols[indicator1] <- cols[2] # cols[2] for outliers at 5% significance level
	pCols[indicator2] <- cols[3] # cols[3] for outliers at 1% significance level
	pPchs <- rep(pchs[1], nrow(datAll))
	pPchs[indicator1] <- pchs[2]
	pPchs[indicator2] <- pchs[3]
	# setup grid
	gg <- format_grid(dose1)
	top <- gg$top
	bot <- gg$bot
	xGrid <- gg$xGrid
	#browser()
	# initialize the plot: the scaled values for ctr and trt
	par(bg="white")
	## the actual data points
	if(is.na(xlim[1])) xlim <- log10(c(bot, top))
	#browser()
	if(style!='simple') { # full or points
		# points for dose>0
		#browser()
		plot(log10(dose1), trtScaled, ylim=ylim, xlim=xlim,
           xlab=xlab, ylab=ylab, main=main, col=pCols[isTrt], pch=pPchs[isTrt], cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab)
		## more intuitive axis 
	} else {
		#browser()#simple no points
		plot(log10(dose1), trtScaled, ylim=ylim, xlim=xlim,
           xlab=xlab, ylab=ylab, main=main, cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab, col=pCols[isTrt], pch=pPchs[isTrt], type='n')
	}	
	# only add the points if style is  full
	if(style!='simple'){ # full or points
		# points for dose=0
		if(hasCtrl)	points(rep(log10(bot), length(ctrScaled)), ctrScaled, pch=8)	
	}
	abline(h=h, col='grey')
	### add curve from the fitted model
	y <- predict(x, newData=xGrid)
	#browser()
	# modified on 02/23/2014: use a subset of points so that the pdf figure will not be too large
	## too many points and leads to a large pdf: use a subset of the points
	ind1 <- which(diff(log10(xGrid))>1e-3) # at log10 dose scale, a step length=1e-3 should be small enough to produce smooth curves
	ind2 <- floor(seq(max(ind1)+1, length(xGrid), length.out=1000))
	indSel <- c(ind1, ind2)
	#lines(log10(xGrid), y, col=col, lwd=lwd)
	#browser()
	lines(log10(xGrid)[indSel], y[indSel], col=col, lwd=lwd)
	## add legend
	#if(addLegend) legend("topright", c("okay", "95%", "99%"), col=cols, pch=pchs) ## use significance level to remove confusion (modified on 09/03/2013)
	if(style!='full')
		addLegend <- FALSE # simple style removes the outlier status
	# if alpha=1, of course no need for legend
	if(x@alpha==1) addLegend <- FALSE
	if(addLegend) legend("topright", c("okay", "5%", "1%"), col=cols, pch=pchs, bty=bty)
	#browser()
})
## @rdname drFit-class
#' lines method for drFit object
#' @param x a drFit object	   
#' @param col line color
#' @param lwd line width  
#' @param show_points whether to add points for dose-response paires (dose not 0) 
#' @param pcol color for points; only relevant if show_points is TRUE
#' @param pch pch for points; only relevant if show_points is TRUE	   
#' @param ... additional parametrs passed to generic lines() function
#' @aliases lines,drFit-method
#' @export 
setMethod('lines', signature(x='drFit'),
    function(x, col=5, lwd=2, show_points=FALSE, pcol='black', pch=16, ...) {
	resPrepDR <- x@info$resPrepDR # result of prep DR data 
	datAll <- resPrepDR$datAll # scaled data including control (if available) and outliers
	isTrt <- resPrepDR$isTrt # global indicator if data is treatment
	dose <- datAll$dose # may include 0 dose
	dose1 <- dose[isTrt]
	trtScaled <- datAll$response[isTrt] # response at nonzero dosage
	# setup grid	
	gg <- format_grid(dose1)
	top <- gg$top
	bot <- gg$bot
	xGrid <- gg$xGrid
	y <- predict(x, newData=xGrid)
	# modified on 02/23/2014: use a subset of points so that the pdf figure will not be too large
	ind1 <- which(diff(log10(xGrid))>1e-3)
	ind2 <- floor(seq(max(ind1)+1, length(xGrid), length.out=1000))
	indSel <- c(ind1, ind2)
	lines(log10(xGrid)[indSel], y[indSel], col=col, lwd=lwd, ...)
	# sometimes the added lines can be accompanied with points to show original data
	if(show_points){
		points(log10(dose1), trtScaled, col=pcol, pch=pch)
	}
	#lines(log10(xGrid), y, col=col, lwd=lwd)	  
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
