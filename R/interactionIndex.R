#'  The drug combination data between SCH66336 and 4-HPR. This is an all-possible combination design.
#'
#' This dataset is published in Lee et al, 2007, Journal of Biopharmaceutical Statistics, Table 1
#'
#' \itemize{
#'   \item schd SCH66336 dose
#'   \item hpr 4-HPR dose
#'   \item y1 response
#' }
#'
#' @format A data frame with 30 rows and 3 columns
#' @source \url{https://biostatistics.mdanderson.org/SoftwareDownload/}
#' @references Lee, J. Jack, et al. Interaction index and different 
#'    methods for determining drug interaction in combination therapy. Journal of biopharmaceutical statistics 17.3 (2007): 461-480.
#' @name nl22B2
NULL


#' The drug combination data of squamous cell carcinoma cells (UMSCC22B) surviving after 72 hours of treatment by single and
#' combination dose levels of SCH66336 and 4-HPR. This is a fixed-ratio design.
#'
#' This dataset is published in Lee et al, 2009, Statistics in biopharmaceutical research, Table 2
#'
#' \itemize{
#'   \item dose1 SCH66336 dose
#'   \item dose2 4-HPR dose
#'   \item e response
#' }
#'
#' @format A data frame with 13 rows and 3 columns
#' @source \url{https://biostatistics.mdanderson.org/SoftwareDownload/}
#' @references Lee, J. J., & Kong, M. (2009). Confidence intervals 
#'   of interaction index for assessing multiple drug interaction. Statistics in biopharmaceutical research, 1(1), 4-17.
#' @name UMSCC22B
NULL


#' detect if a experiment is fixed-ray design and prepare input for CI computation
#' @param d1 d1
#' @param d2 d2
#' @param e e
#' @param tol tolerance in declaring fixed ratio. 
#' @param d2.d1.force the ratio d2/d1 to be forced in fitting the linear model; this is only effective if the data is not a fixed ratio design
#'  usually, the all-possible combination design; specifying d2.d1.force will use the specified ratio (from available ratios) to fit the model and generate dose-response curve.
#' @return a list
#' @export
detect_ray_design <- function(d1, d2, e, tol=0.001, d2.d1.force=NA){
	# Here d2.d1.force use ray design method to fix grid data. What is the rationale and how to justify it?
	ind_eAbnorm <- e<0 | e>1
	if(sum(ind_eAbnorm)>0) {
		warning(sprintf('%d responses not in the range [0, 1]; truncated automatically', sum(ind_eAbnorm)))
		e[e<0] <- 0
		e[e>1] <- 1
	}
	ind1 <- d1!=0 & d2==0
	ind2 <- d1==0 & d2!=0
	#ind12_ <- d1!=0 & d2!=0 # where the dose differs
	ind12 <- d1!=0 & d2!=0 # where the dose differs
	ratio_d <- d2[ind12]/d1[ind12] # pointwise ratio
	#browser()
	isRayDesign <- all.equal(ratio_d, rep(mean(ratio_d), length(ratio_d)), tolerance=tol) # TRUE or a string
	if(isRayDesign!=TRUE) isRayDesign <- FALSE # make string as false
	# d1, e1, d2, e2, d12, e12 added for CI.delta() estimation; d12 is just the sum; d2.d1 is the mean of ratio when not a ray design
	if(isRayDesign!=TRUE){
		# not a ray design
		if(is.na(d2.d1.force)){
			d2.d1.force <- 1
		}
		if(!d2.d1.force %in% ratio_d){
			stop(sprintf("This is not a fixed ratio design; the specified d2.d1 ratio is not supported in the data!"))
		}
		# force d2.d1 as the specified ratio or the default
		d2.d1 <- d2.d1.force
	} else {
		d2.d1 <- mean(ratio_d) # this is not necessary for fixed ratio design
	}
	list(drMat=data.frame(d1=d1, d2=d2, e=e), isRayDesign=isRayDesign, ind_eAbnorm=ind_eAbnorm,  
		d1=d1[ind1], e1=e[ind1], d2=d2[ind2], e2=e[ind2], d12=d1[ind12]+d2[ind12], e12=e[ind12], ind12=ind12, d2.d1=d2.d1, d2.d1_avail=ratio_d, d2.d1.force=d2.d1.force)
}
#res_design <- detect_ray_design(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1)
#detect_ray_design(d1=dose1, d2=dose2, e=fa)

#### todo: deal with all-possible combination case -- it is not solved completely in Lee 2009 paper. 
#### mostly borrowed from chou8491.SSC; but also borrowed from CI_IIV2.SSC
#' fit median-effect model for 2-drug combination with a fixed ratio desgin
#'
#' notice that the median-effect model is a special case of the sigma Emax model;
#' 3 linear models are fitted; fixed ratio thus is required.
#' all 3 linear models are fitted using logit(E) ~ log(d1+d2) using a subset of the data depending on if the data is for durg 1, drug 2 and the mixture
#'
#' @param d1 dose for drug 1 
#' @param d2 dose for drug 2
#' @param e corresponding response in the range [0, 1]
#' @param name1 drug 1 name
#' @param name2 drug 2 name
#' @param d2.d1 d2/d1, the fixed ratio
#' @param base base of logarithm; disabled since the dose-response curve estimation also needs transformation
#' @return a list
#' @export
fit_median_efect <- function(d1, d2, e, name1='Drug A', name2='Drug B', d2.d1, base=exp(1)){
	dose1 <- d1; dose2 <- d2; fa <- e
	#browser()
	fu <- 1-fa
	##     define an indicator for identifying doses satisfying dose2/dose1=d2.d1
	ind_mix <- !(dose1+dose2==0 | fa<=0 | fa>=1) # exclude data outside range of the model
	#totdose <- dose1[ind] + dose2[ind]
	#totdose <- dose1[ind] + dose2[ind]
	totdose <- rep(NA, length(dose1)) # match with input length,i.e. ind.ratio
	totdose[ind_mix] <- (dose1+dose2)[ind_mix] # obtain a total dose to fit 3 linear models
	####
	#resp <- rep(NA, length(fa[ind]))
	#BUG in the original code: length not sum(ind): resp[ind] <- log(fa[ind]/fu[ind]) # logit(E)
	#resp <- log(fa/fu)[ind] # logit(E)
	####
	resp <- rep(NA, length(fa)) # match with input length as in plot(medianEffect$logd,medianEffect$resp
	ind_logit <- fa<=0 | fa>=1 # logit transform for good data; others are NA
	resp[!ind_logit] <- log(fa[!ind_logit]/fu[!ind_logit]) 	
	logd <- rep(NA, length(totdose))
	#logd <- changeBase(log(totdose), oldBase=exp(1), newBase=base) # change base of logarithm done here
	logd <- log(totdose) 
	#log10d <- log10(totdose) # log10(dose)	
	##     define an indicator for identifying doses satisfying dose2/dose1=d2.d1
	# this is very sensitive to mis-ratio;
	## modified 2017/02/24: let the user to make sure the data is fixed ratio; using all data in combo experiment
	#ind.ratio  <- dose1!=0 & abs(dose2-d2.d1*dose1)<0.00001
	ind.ratio  <- dose1!=0 & dose2!=0
	ind2  <- abs(dose2-d2.d1*dose1)<0.0001 # to also include 0 dose for dose-response curve without taking log
	##     Estimate the parameters using median-effect plot for two single drugs and 
	##     their combination at the fixed ratio (dose of drug 2)/(dose of drug 1)=d2.d1.
	#browser()
	lm1 <- lm(resp[dose2==0 & dose1!=0]~logd[dose2==0 & dose1!=0])
	dm1 <- exp(-summary(lm1)$coef[1,1]/summary(lm1)$coef[2,1])
	lm2 <- lm(resp[dose1==0 & dose2!=0]~logd[dose1==0 & dose2!=0])
	dm2 <- exp(-summary(lm2)$coef[1,1]/summary(lm2)$coef[2,1])
	lmcomb <- lm(resp[ind.ratio]~logd[ind.ratio])
	dmcomb <- exp(-summary(lmcomb)$coef[1,1]/summary(lmcomb)$coef[2,1]) 
	#browser()
	# some non-essential data are attached to facilitate visualization
	list(lm1=lm1, dm1=dm1, lm2=lm2, dm2=dm2, lmcomb=lmcomb, dmcomb=dmcomb, ind.ratio=ind.ratio, ind2=ind2, base=base, 
		logd=logd, dose1=dose1, dose2=dose2, resp=resp, totdose=totdose, fa=fa, name1=name1, name2=name2, d2.d1=d2.d1, d1=d1, d2=de, e=e)
}



#' generate plot for the fit from median-effect model
#' 
#' @param medianEffect the fitted result by fit_median_efect()
#' @param type type of plot, options are from c('medianEffect', 'doseResponseCurve', 'contour')
#' @param contour.level levels (between 0 and 1 representing different response values) to draw contour plot
#' @param logd whether to plot logd in Dose-response curve; only effective if type=='doseResponseCurve'
#' @param align whether to align dose range. Imagine that 2 drugs are mixed in fixed ratio, the doses in one drug will have low doses that the other drug and mixture do not have 
#'  which some biologists deams as a waste of graphic region. An easy fix is to truncate xlim. this is only only effective if type=='doseResponseCurve'
#' @param cex.legend cex for the legend text'
#' @param legend.position legend position passed to legend() function
plot_median_effect <- function(medianEffect, type=c('medianEffect', 'doseResponseCurve', 'contour'), 
	contour.level=(1:9)/10, logd=FALSE, cex.legend=1, align=TRUE, ylim=c(0, 1), legend.position='topright'){
	if(type=='medianEffect'){
			plot(medianEffect$logd,medianEffect$resp,type='n',xlab=sprintf('Log(Dose)'),ylab='Log(E/(1-E))', 
				main=sprintf('Median Effect Plot \n-- d2.d1=%.1g', medianEffect$d2.d1))
            #abline(lm1,lty=4)
            abline(medianEffect$lm1,lty=3)
            with(medianEffect, points(logd[dose2==0 & dose1!=0],resp[dose2==0 & dose1!=0],pch=1))
            abline(medianEffect$lm2,lty=2)
            with(medianEffect, points(logd[dose1==0 & dose2!=0],resp[dose1==0 & dose2!=0],pch=2))
            abline(medianEffect$lmcomb,lty=1)
            with(medianEffect, points(logd[ind.ratio],resp[ind.ratio], pch=16))
            #title("Median Effect Plots" )
			 with(medianEffect, legend('topright', c(name1, name2, 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n'))
			 #browser()
	} else if(type=='doseResponseCurve'){
			#browser()
			if(logd){
				#dose <- with(medianEffect, seq(0,max(totdose[ind2], na.rm=TRUE),max(totdose[ind2], na.rm=TRUE)/500000)) # need more points here for finer grid
				#browser()
				dose <- noNA(drexplorer:::format_grid(dose1=with(medianEffect, c(totdose, dose1, dose2)), n=300)$xGrid)
				if(align){
					mindose <- 	min(log(dose[dose!=0]), na.rm=T)
				} else {
					dd <- c(dose, with(medianEffect, c(dose1, dose2)))
					mindose <- min(log(dd[dd!=0]), na.rm=T)
            	}
            	#browser()
            	# notice: when plot log dose, at low doses, the predicted response might not be present in medianEffect due to grid selection. 
				with(medianEffect, plot(log(totdose[ind2]),fa[ind2],type='n',xlab='Log(Dose)',ylab='Relative viability', ylim=ylim, xlim=c(mindose, max(log(dose)))))
            	#y1<-with(medianEffect, 1.0/(1.0+(dm1/dose)^(summary(lm1)$coef[2,1])))
            	with(medianEffect, lines(log(dose), with(medianEffect, 1.0/(1.0+(dm1/dose)^(summary(lm1)$coef[2,1]))),type="l",lty=3))
            	with(medianEffect, points(log(dose1[dose2==0]),fa[dose2==0],pch=1))
            	with(medianEffect, lines(log(dose), 1.0/(1.0+(dm2/dose)^(summary(lm2)$coef[2,1])),type="l",lty=2))
            	with(medianEffect, points(log(dose2[dose1==0]),fa[dose1==0],pch=2))
            	with(medianEffect, lines(log(dose), 1.0/(1.0+(dmcomb/dose)^(summary(lmcomb)$coef[2,1])),type="l",lty=1))
            	with(medianEffect, points(log(totdose[ind.ratio]),fa[ind.ratio], pch=16))
            	title( "Dose-Response Curves" )
            	# legend('topright', c(var.name[1:2], 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n')		
		    	with(medianEffect, legend(legend.position, c(name1, name2, 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n', cex=cex.legend))	
			} else {
				dose <- with(medianEffect, seq(0,max(totdose[ind2], na.rm=TRUE),max(totdose[ind2], na.rm=TRUE)/100))
            	with(medianEffect, plot(totdose[ind2],fa[ind2],type='n',xlab='Dose',ylab='Relative viability', ylim=ylim, xlim=c(0, max(dose))))
            	#y1<-with(medianEffect, 1.0/(1.0+(dm1/dose)^(summary(lm1)$coef[2,1])))
            	with(medianEffect, lines(dose, with(medianEffect, 1.0/(1.0+(dm1/dose)^(summary(lm1)$coef[2,1]))),type="l",lty=3))
            	with(medianEffect, points(dose1[dose2==0],fa[dose2==0],pch=1))
            	with(medianEffect, lines(dose, 1.0/(1.0+(dm2/dose)^(summary(lm2)$coef[2,1])),type="l",lty=2))
            	with(medianEffect, points(dose2[dose1==0],fa[dose1==0],pch=2))
            	with(medianEffect, lines(dose, 1.0/(1.0+(dmcomb/dose)^(summary(lmcomb)$coef[2,1])),type="l",lty=1))
            	with(medianEffect, points(totdose[ind.ratio],fa[ind.ratio], pch=16))
            	title( "Dose-Response Curves" )
            	# legend('topright', c(var.name[1:2], 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n')		
		    	with(medianEffect, legend(legend.position, c(name1, name2, 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n', cex=cex.legend))
			}
	} else if(type=='contour'){
			##        The contour plot of raw data.
            #browser()
			temp1 <- sort(unique(medianEffect$dose1))
		    temp2 <- sort(unique(medianEffect$dose2))
			isAllPossible <- length(temp1)*length(temp2)==length(medianEffect$fa) # all possible combination design
			if(isAllPossible){
				obsfa <- matrix(medianEffect$fa, length(temp1),length(temp2), byrow=T)
				contour(temp1, temp2, obsfa, levels=contour.level, xlab=sprintf('%s dose', medianEffect$name1), ylab=sprintf('%s dose', medianEffect$name2)) 
				title("Contour Plot")
			} else {
				#browser()#
				ww <- sprintf('%s has %d doses; %s has %d doses; only %d total responses; This is not a all-possible design and thus no contour plot can be constructed!', 
					medianEffect$name1, length(temp1), 
					medianEffect$name2, length(temp2), length(medianEffect$fa))
				warning(ww)
			}
            
	}
}
#plot_median_effect(fit_allPoss$medianEffect, type='medianEffect')
#plotIAI(fit_allPoss, type='medianEffect', mode='both')



#' Estimate the interaction index as well as its confidence interval using delta method for fixed ratio design
#'
#' this code is extracted from the source code distributed at https://biostatistics.mdanderson.org/SoftwareDownload/
#'
#' Two papers have been published by Lee et al, one in 2007 (Lee2007) and on in 2009 (Lee2009). The Lee2007 paper 
#' described five methods to assess interaction: (1) Lowewe additivity model using interaction index (IAI) (2) Model of Greco et al 1990.
#' This approach uses \deqn{\alpha}{alpha} as the metric and it can be related to IAI (3) Model of Machado and Robinson which uses a metric denoted
#' as \deqn{\eta}{eta} (4) Model of Plummer and Short which can also be linked to IAI through the parameter \deqn{\beta_4}{beta_4} (5) Model of
#' Carter et al that can be linked to IAI through the parameter \deqn{\beta_{12}}{beta_12}. For more details of these models, please refer to Lee2007.
#' 
#' The Lee2009 paper provided generalization of IAI to multiple drugs using Lowewe additivity model and assumption of Chou and Talalay's median effect
#' equation. The Chou and Talalay's median effect equation can be expressed as: \deqn{log(E/(1-E))=m(log d - log Dm)}{log(E/(1-E))=m(log d - log Dm)} 
#' where E is the effect at dose d for a
#' compound whose median effective dose. 
#'
#' Some notes about experiment design. Usually the data is either fixed ratio design (ray design) or grid design which means all-possible combination of drug concentrations
#' between two drugs are available. The Lee2007 paper provided an example of grid design. However, specific fixed ratio is used to fit the median effect model which is the basis to estimate IAI. The Lee2009 paper
#' considered with fixed ratio design.  
#'
#' @param d1 dose for drug 1 
#' @param d2 dose for drug 2
#' @param e corresponding response in the range [0, 1]
#' @param E a vector of responses (between 0 and 1) where IAI and confidence interval are to be computed from.                                           
#' @param name1 name of drug 1                               
#' @param name2 name of drug 2                               
#' @param alpha significance level of confidence interval             
#' @param d2.d1.force a ratio passed to detect_ray_design() function so as to specify a fixed ratio for grid design   
#' @return a data frame with columns IAI, IAI.low, IAI.up, E, dx1 (corresponding dose of drug 1), dx2 (corresponding dose of drug 1), 
#' @export
#' @references Lee, J. J., & Kong, M. (2009). Confidence intervals 
#'   of interaction index for assessing multiple drug interaction. Statistics in biopharmaceutical research, 1(1), 4-17.
#' @references Lee, J. Jack, et al (2007). Interaction index and different 
#'    methods for determining drug interaction in combination therapy. Journal of biopharmaceutical statistics 17.3 461-480.
fitIAI <- function(d1, d2, e, E=seq(0.05, 0.95, 0.005), name1='Drug A', name2='Drug B', alpha=0.05, d2.d1.force=NA, tol=0.01){
	res_design <- detect_ray_design(d1=d1, d2=d2, e=e, d2.d1.force=d2.d1.force, tol=tol)
	#browser()
	# fit median-effect model
	medianEffect <- fit_median_efect(d1=d1, d2=d2, e=e, d2.d1=res_design$d2.d1, name1=name1, name2=name2)
	if(res_design$isRayDesign!=TRUE) { # not fixed ratio design, fitCI is NULL
		#warning("This is not a ray design; computation aborted!")
		warning("This is not a ray design; computation is done with method derived from fixed ratio!")
		#fitCI <- NULL
	} else {
		warning(sprintf("Fixed ratio design detected. Fixed ratio is %.3f\n", res_design$d2.d1))
	}
	## fit CI
	fitCI <- with(res_design, 
		CI.delta(d1=d1, e1=e1, d2=d2, e2=e2, d12=d12, 
				e12=e12, d2.d1, E=E, alpha=alpha)
				)
	
	#}
	meta <- list(d1=d1, d2=d2, e=e, E=E, alpha=alpha, d2.d1=res_design$d2.d1, isRayDesign=res_design$isRayDesign)
	list(CI=fitCI, meta=meta, medianEffect=medianEffect)
	#browser()
}
#fitIAI(d1=dose1, d2=dose2, e=fa)





#' Visualize interaction index 
#'
#' @param fit the fitted result from fitIAI()
#' @param type type of plot from c('IAI', 'medianEffect', 'doseResponseCurve', 'contour')
#' @param contour.level contour level. only effective if type='contour'
#' @param ylim y axis limit
#' @param mode specify if to plot against response, dose or both. only effective if type=='IAI'. can be either 'response', 'dose', or 'both'
#' @param logd logd passed to plot_median_effect
#' @param align align passed to plot_median_effect
#' @param cex.legend cex.legend passed to plot_median_effect
#' @param legend.position legend.position passed to plot_median_effect
#' @export
plotIAI <- function(fit, type=c('IAI', 'medianEffect', 'doseResponseCurve', 'contour'), contour.level=(1:9)/10, ylim=NULL, mode='both', 
	logd=FALSE, align=TRUE, cex.legend=1, legend.position='topright'){
	if(type=='IAI' & !is.null(fit$CI)){
		resCI <- fit$CI
		#browser()
		#if(is.null(ylim)) ylim <- range(resCI$IAI)
		if(is.null(ylim)) ylim <- c(0.01, 10)
		op <- par(mar=c(8, 5, 4, 4) + 0.1)
		lty_E <- 4
		# CI
		ind.CId.l <- resCI$IAI.low >=min(resCI$IAI)
		ind.CId.u <- resCI$IAI.up <= max(resCI$IAI)+1
		xlim_log10d <- range(log10(resCI$dx12))
		#col_d <- '#38383A'
		col_d <- 'blue2'
		lty_d <- 5
		if(mode=='both'){
		# IAI ~ E
		fa <- fit$meta$e
		E <- resCI$E
		#browser()
		with(resCI, plot(E, IAI, log="y", type="l", xlab="", xaxt='n', ylab="Interaction Index", 
			xlim=c(min(E, fa),max(E, fa)), ylim=ylim,cex=0.6, main='Interaction plot')) 
		axis(1, at=pretty(range(E)), col="black",lwd=1)
		mtext("Relative viability",side=1,col="black",line=2)
		abline(h=1, lty=4)
		lines(E[ind.CId.l], resCI$IAI.low[ind.CId.l], lty=lty_E)  
		lines(E[ind.CId.u], resCI$IAI.up[ind.CId.u], lty=lty_E)
		### add dose info	
		par(new=TRUE)
		## IAI ~ dose
		plot(resCI$dx12, resCI$IAI, log="xy", type="l", ylim=ylim, axes=F, xlab='', ylab='', lty=lty_d, col=col_d)
		axis(1, xlim=xlim_log10d,lwd=1,line=3.2, col=col_d, col.axis=col_d, lty=lty_d)
		mtext("Dose",side=1,col=col_d,line=5.2)
		lines(resCI$dx12[ind.CId.l], resCI$IAI.low[ind.CId.l], lty=lty_d, col=col_d)  
		lines(resCI$dx12[ind.CId.u], resCI$IAI.up[ind.CId.u], lty=lty_d, col=col_d)
		} else if(mode=='response'){
		# IAI ~ E
		fa <- fit$meta$e
		E <- resCI$E
		plot(resCI$E, resCI$IAI, log="y", type="l", xlab="", xaxt='n', ylab="Interaction Index", xlim=c(min(E, fa),max(E, fa)), ylim=ylim,cex=0.6, main='Interaction plot') 
		abline(h=1, lty=4)
		axis(1, at=pretty(range(resCI$E)), col="black",lwd=1)
		mtext("Relative viability",side=1,col="black",line=2)
		abline(h=1, lty=4)
		lines(E[ind.CId.l], resCI$IAI.low[ind.CId.l], lty=lty_E)  
		lines(E[ind.CId.u], resCI$IAI.up[ind.CId.u], lty=lty_E)
		} else if(mode=='dose'){
		## IAI ~ dose
		#plot(resCI$dx12, resCI$IAI, log="xy", type="l", ylim=ylim, xaxt='n', xlab='', ylab='Interaction Index', lty=lty_d, col=col_d)
		#browser()
		plot(resCI$dx12, resCI$IAI, log="xy", type="l", xlab="Dose", ylab="Interaction Index", ylim=ylim, cex=0.6, main='Interaction plot') 
		abline(h=1, lty=4)
		#axis(1, at=pretty(range(resCI$dx12)),col="black",lwd=1)
		#mtext("Dose (axis in log10 scale)",side=1,col="black",line=2)
		#mtext("Dose (axis in log10 scale)",side=1,col=col_d,line=5.2)
		lines(resCI$dx12[ind.CId.l], resCI$IAI.low[ind.CId.l], lty=lty_d)  
		lines(resCI$dx12[ind.CId.u], resCI$IAI.up[ind.CId.u], lty=lty_d)
		}
		par(op)
	} else {
	#browser()
		medianEffect <- fit$medianEffect
		if(type=='contour'){
			plot_median_effect(medianEffect, type=type, contour.level=contour.level) 
		} else {
			plot_median_effect(medianEffect, type=type, logd=logd, align=align, ylim=ylim, legend.position=legend.position, cex.legend=cex.legend)   
		}	
	}
}

### downloaded from: https://biostatistics.mdanderson.org/SoftwareDownload/
##################################################################################################################
###                                                                                                            ###
###   Confidence Bound of Interaction Indices vs Effects at a fixed ray based on Lee and Kong (2006)           ###
###           Section 3 in this paper                                                                          ###
###                                                                                                            ###
###  INPUT:                                                                                                    ###
###  d1 and e1:     observed doses and effects for drug 1                                                      ###
###  d2 and e2:     observed doses and effects for drug 2                                                      ###
###  d12 and e12:   observed doses and effects for mixture at the fixed ratio d2/d1=d2.d1                      ###
###  E:             fixed effects, their corresponding interaction indices and confidence intervals are estimated.
###  alpha:         1-alpha is the size of the confidence intervals, alpha has the default value of 0.05.      ###
###                                                                                                            ###
###  OUTPUT:                                                                                                   ###
###  ii:            the estimated interaction indices corresponding to the input effects E                     ###
###  ii.low:        the estimated lower confidence intervals for ii                                            ###                        
###  ii.up:         the estimated upper confidence intervals for ii                                            ###                        
###                                                                                                            ###
##################################################################################################################  
#
# originally called: CI.delta from CI_IIV2.SSC                          
#
#
# truncate effect so that it will never be 0 or 1

truncate_effect <- function(e, min=1e-4, max=1-1e-4){
	e[e<min] <- min
	e[e>max] <- max
	e
}

 
#' Estimate the interaction index as well as its confidence interval using delta method for fixed ratio design
#'
#' This code is extracted from the source code distributed at https://biostatistics.mdanderson.org/SoftwareDownload/
#' The paper is  
#' Two versions are included in the source code: CI_IIV2 2008.SSC and CI_IIV2.SSC
#' At first we implement CI_IIV2.SSC; this gives wider CI band; so we decide to use  CI_IIV2 2008.SSC. Further, CI_IIV2 2008.SSC has
#' a more recent date and thus most updated.
#'
#' @param d1 dose for drug 1 where drug 2 has dose equal to 0                                         
#' @param e1 corresponding scaled response for drug 1, must between 0 and 1                                          
#' @param d2 dose for drug 2 where drug 1 has dose equal to 0                                               
#' @param e2 corresponding scaled response for drug 2, must between 0 and 1                                
#' @param d12 the sum of doses from drug 1 and drug 2 when both drugs are administered                                          
#' @param e12 corresponding scaled response when both drugs are administered                                          
#' @param d2.d1 the ratio of the ray design. This is the fixed ratio of dose 2 divided by dose 1
#' @param E a vector of responses (between 0 and 1) where IAI and confidence interval are to be computed from.                                           
#' @param alpha significance level of confidence interval             
#' @param min response values close to 0 and 1 will lead to extreme values when logit transformation is applied. min is used to truncate the response value so that values less than min are
#' truncated at min          
#' @param max response values close to 0 and 1 will lead to extreme values when logit transformation is applied. max is used to truncate the response value so that values larger than max are
#' truncated at max          
#' @return a data frame with columns IAI, IAI.low, IAI.up, E, dx1 (corresponding dose of drug 1), dx2 (corresponding dose of drug 1), 
#'  dx12 (corresponding dose of combined drug, same as definition of d12)
#' @references Lee, J. J., & Kong, M. (2009). Confidence intervals 
#'   of interaction index for assessing multiple drug interaction. Statistics in biopharmaceutical research, 1(1), 4-17.
CI.delta <- function(d1, e1, d2, e2, d12, e12, d2.d1, E, alpha=0.05, min=0.02, max=0.98)
{
     # min and max are used to truncate the effect so that logit transform would not be obsurd; the original code is not very robust!
	 #browser()
	 e1 <- truncate_effect(e1, min=min, max=max)
	 e2 <- truncate_effect(e2, min=min, max=max)
	 e12 <- truncate_effect(e12, min=min, max=max)
	 # if e1, e2 or e12 has value of 0 or 1, the logit would be Inf and leads to error for lm.
	 # thus use safe log
	 lm1 <- lm(log(e1/(1-e1))~log(d1))
     dm1 <- exp(-summary(lm1)$coef[1,1]/summary(lm1)$coef[2,1])
     lm2 <- lm(log(e2/(1-e2))~log(d2))
     dm2 <- exp(-summary(lm2)$coef[1,1]/summary(lm2)$coef[2,1])
     lmcomb <- lm(log(e12/(1-e12))~log (d12))
     dm12 <- exp(-summary(lmcomb)$coef[1,1]/summary(lmcomb)$coef[2,1]) 
     Dx1 <- dm1*(E/(1-E))^(1/summary(lm1)$coef[2,1])
     Dx2 <- dm2*(E/(1-E))^(1/summary(lm2)$coef[2,1])
     dx12 <- dm12*(E/(1-E))^(1/summary(lmcomb)$coef[2,1])
     iix <- (dx12/(1+d2.d1))/Dx1+(dx12*d2.d1/(1+d2.d1))/Dx2
     lm1.s <-summary(lm1)
     lm2.s <-summary(lm2)
     lm12.s <-summary(lmcomb)
     c1 <- 1.0/lm1.s$coef[2,1]^2*lm1.s$coef[1,2]^2
     temp <- - mean(log(d1))*lm1.s$coef[2,2]^2
    ### temp <- lm1.s$coef[1,2]*lm1.s$coef[2,2]*lm1.s$cor[1,2]   ### covariance of b0 and b1
     c1 <- c1+2.0*(log(E/(1-E))-lm1.s$coef[1,1])/lm1.s$coef[2,1]^3*temp
     c1 <- c1+(log(E/(1-E))-lm1.s$coef[1,1])^2/lm1.s$coef[2,1]^4*lm1.s$coef[2,2]^2
     c2 <- 1.0/lm2.s$coef[2,1]^2*lm2.s$coef[1,2]^2
     temp <- - mean(log(d2))*lm2.s$coef[2,2]^2
    ### temp <- lm2.s$coef[1,2]*lm2.s$coef[2,2]*lm2.s$cor[1,2]   ### covariance of b0 and b1
     c2 <- c2+2.0*(log(E/(1-E))-lm2.s$coef[1,1])/lm2.s$coef[2,1]^3*temp
     c2 <- c2+(log(E/(1-E))-lm2.s$coef[1,1])^2/lm2.s$coef[2,1]^4*lm2.s$coef[2,2]^2
     c12 <- 1.0/lm12.s$coef[2,1]^2*lm12.s$coef[1,2]^2
     temp <- - mean(log(d12))*lm12.s$coef[2,2]^2
   ### temp <- lm12.s$coef[1,2]*lm12.s$coef[2,2]*lm12.s$cor[1,2]   ### covariance of b0 and b1
     c12 <- c12+2.0*(log(E/(1-E))-lm12.s$coef[1,1])/lm12.s$coef[2,1]^3*temp
     c12 <- c12+(log(E/(1-E))-lm12.s$coef[1,1])^2/lm12.s$coef[2,1]^4*lm12.s$coef[2,2]^2
     var.ii <-((dx12/Dx1)^2*c1+(dx12*d2.d1/Dx2)^2*c2+(1.0/Dx1+d2.d1/Dx2)^2*dx12^2*c12)/(1+d2.d1)^2 
     t975 <- qt(1-alpha/2,length(d1)+length(d2)+length(d12)-6)
     iix.low1 <- iix*exp(-t975*var.ii^0.5/iix)
     iix.up1 <- iix*exp(t975*var.ii^0.5/iix)
	 #browser()
     #return(list(ii=iix, ii.low=iix.low1, ii.up=iix.up1))
     # add the corresponding dose for visualization
	 IAI <- data.frame(E=E, IAI=iix, IAI.low=iix.low1, IAI.up=iix.up1, dx1=Dx1, dx2=Dx2, dx12=dx12)
	 #browser()
	 return(IAI)
}

