if(FALSE){
setwd('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/')
library(devtools)
build('drexplorer')
install('drexplorer')

detach("package:drexplorer", unload=TRUE)
library(drexplorer)

load_all('../drexplorer')


##
detach("package:drexplorerExtra", unload=TRUE)
library(drexplorerExtra)


source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/drexplorer.R'))
source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/drexplorerAdded.R'))
source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/aux_from_Extra.R'))
source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/GUI_1_source_v2.R'))
source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/interactionIndex.R'))
source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/GUI_2_source.R'))


}

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
		d2.d1 <- mean(ratio_d)
	}
	list(drMat=data.frame(d1=d1, d2=d2, e=e), isRayDesign=isRayDesign, ind_eAbnorm=ind_eAbnorm,  
		d1=d1[ind1], e1=e[ind1], d2=d2[ind2], e2=e[ind2], d12=d1[ind12]+d2[ind12], e12=e[ind12], ind12=ind12, d2.d1=d2.d1, d2.d1_avail=ratio_d, d2.d1.force=d2.d1.force)
}
#res_design <- detect_ray_design(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1)
#detect_ray_design(d1=dose1, d2=dose2, e=fa)


# fit median-effect model for 2-drug combination with a fixed ratio desgin
#'
#' notice that the median-effect model is a special case of the sigma Emax model;
#' 3 linear models are fitted; fixed ratio thus is required
#' mostly borrowed from chou8491.SSC; but also borrowed from CI_IIV2.SSC
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
#### todo: deal with all-possible combination case
fit_median_efect <- function(d1, d2, e, name1='Drug A', name2='Drug B', d2.d1, base=exp(1)){
	dose1 <- d1; dose2 <- d2; fa <- e
	#browser()
	fu <- 1-fa
	##     define an indicator for identifying doses satisfying dose2/dose1=d2.d1
	ind_mix <- !(dose1+dose2==0 | fa<=0 | fa>=1) # exclude data outside range of the model
	#totdose <- dose1[ind] + dose2[ind]
	#totdose <- dose1[ind] + dose2[ind]
	totdose <- rep(NA, length(dose1)) # match with input length,i.e. ind.ratio
	totdose[ind_mix] <- (dose1+ + dose2)[ind_mix]
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
	ind.ratio  <- dose1!=0 & abs(dose2-d2.d1*dose1)<0.00001
	ind2  <- abs(dose2-d2.d1*dose1)<0.0001 # to also include 0 dose for dose-response curve without taking log
	##     Estimate the parameters using median-effect plot for two single drugs and 
	##     their combination at the fixed ratio (dose of drug 2)/(dose of drug 1)=d2.d1.
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
plot_median_effect <- function(medianEffect, type=c('medianEffect', 'doseResponseCurve', 'contour'), contour.level=(1:9)/10){
	if(type=='medianEffect'){
	#browser()
		plot(medianEffect$logd,medianEffect$resp,type='n',xlab=sprintf('Log(Dose)'),ylab='Log(E/(1-E))', 
				main=sprintf('Median Effect Plot \n-- d2.d1=%.1f', medianEffect$d2.d1))
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
			dose <- with(medianEffect, seq(0,max(totdose[ind2], na.rm=TRUE),max(totdose[ind2], na.rm=TRUE)/100))
            with(medianEffect, plot(totdose[ind2],fa[ind2],type='n',xlab='Dose',ylab='Relative viability', ylim=c(0, 1), xlim=c(0, max(dose))))
            #y1<-with(medianEffect, 1.0/(1.0+(dm1/dose)^(summary(lm1)$coef[2,1])))
            with(medianEffect, lines(dose, with(medianEffect, 1.0/(1.0+(dm1/dose)^(summary(lm1)$coef[2,1]))),type="l",lty=3))
            with(medianEffect, points(dose1[dose2==0],fa[dose2==0],pch=1))
            with(medianEffect, lines(dose, 1.0/(1.0+(dm2/dose)^(summary(lm2)$coef[2,1])),type="l",lty=2))
            with(medianEffect, points(dose2[dose1==0],fa[dose1==0],pch=2))
            with(medianEffect, lines(dose, 1.0/(1.0+(dmcomb/dose)^(summary(lmcomb)$coef[2,1])),type="l",lty=1))
            with(medianEffect, points(totdose[ind.ratio],fa[ind.ratio], pch=16))
            title( "Dose-Response Curves" )
           # legend('topright', c(var.name[1:2], 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n')		
		   with(medianEffect, legend('topright', c(name1, name2, 'Mixture'), pch=c(1, 2, 16), lty=c(3, 2, 1), bty='n'))
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
#' @param d1 dose for drug 1 
#' @param d2 dose for drug 2
#' @param e corresponding response in the range [0, 1]
#' @param E a vector of responses (between 0 and 1) where IAI and confidence interval are to be computed from.                                           
#' @param alpha significance level of confidence interval             
#' @return a data frame with columns IAI, IAI.low, IAI.up, E, dx1 (corresponding dose of drug 1), dx2 (corresponding dose of drug 1), 
#' @export
fitIAI <- function(d1, d2, e, E=seq(0.05, 0.95, 0.005), name1='Drug A', name2='Drug B', alpha=0.05){
	res_design <- detect_ray_design(d1=d1, d2=d2, e=e)
	# fit median-effect model
	medianEffect <- fit_median_efect(d1=d1, d2=d2, e=e, d2.d1=res_design$d2.d1, name1=name1, name2=name2)
	if(res_design$isRayDesign!=TRUE) { # not fixed ratio design, fitCI is NULL
		#warning("This is not a ray design; computation aborted!")
		warning("This is not a ray design; computation is done with method derived from fixed ratio!")
		#fitCI <- NULL
	} #else {
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
#' @param mode specify if to plot against response, dose or both. only effective if type=='IAI'. can be either 'response', 'dose', or 'both'
#' @export
plotIAI <- function(fit, type=c('IAI', 'medianEffect', 'doseResponseCurve', 'contour'), contour.level=(1:9)/10, mode='both'){
	if(type=='IAI' & !is.null(fit$CI)){
		resCI <- fit$CI
		#browser()
		ylim <- range(resCI$IAI)
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
		plot(resCI$dx12, resCI$IAI, log="xy", type="l", xlab="Dose", ylab="Interaction Index", ylim=ylim,cex=0.6, main='Interaction plot') 
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
			plot_median_effect(medianEffect, type=type)   
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
 
#' Estimate the interaction index as well as its confidence interval using delta method for fixed ratio design
#'
#' this code is extracted from the source code distributed at https://biostatistics.mdanderson.org/SoftwareDownload/
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
#' @return a data frame with columns IAI, IAI.low, IAI.up, E, dx1 (corresponding dose of drug 1), dx2 (corresponding dose of drug 1), 
#'  dx12 (corresponding dose of combined drug, same as definition of d12)
CI.delta <- function(d1, e1, d2, e2, d12, e12, d2.d1, E, alpha=0.05)
{
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
	 IAI <- data.frame(IAI=iix, IAI.low=iix.low1, IAI.up=iix.up1, E=E, dx1=Dx1, dx2=Dx2, dx12=dx12)
	 #browser()
	 return(IAI)
}



if(FALSE){
	# example
	#par(mfrow=c(1,2),mai=c(0.4,0.4,0.4,0.0), mgp=c(1.0, 0.2, 0))
# data published in Table 2 of: Confidence Intervals of Interaction Index for Assessing Multiple Drug Interaction
dose1 <-c(0.1,0.5,1,2,4,0,0,0,0, 0.1, 0.5, 1, 2)
dose2 <-c(0,0,0,0,0,0.1,0.5,1,2, 0.1, 0.5, 1, 2)
fa <-c(0.6701,0.6289,0.5577, 0.455,0.3755,0.7666, 0.5833,0.5706,0.4934,0.6539,0.4919,0.3551,0.2341)
UMSCC22B <- data.frame(dose1=dose1, dose2=dose2, e=fa)
colnames(UMSCC22B) <- c('SCH66336', '4HPR', 'Resp')
#save(UMSCC22B, file='UMSCC22B.RData')
d2.d1 <-1
#####  median effect plots (Figure 3, Panel A )               #########
median.effect.plots(dose1, dose2, fa, name1="SCH66336", name2="4HPR", d2.d1=1/1)
####  Construct confidence intervals and confidence bounds  (Figure 3, Panel B )  ###
ind1 <- dose1!=0 & dose2==0
ind2 <- dose1==0 & dose2!=0
ind12 <- dose1!=0 & abs(dose2-d2.d1*dose1)<0.00001
E<- seq(0.05, 0.95, 0.005)
Er <-c(E,fa[ind12])
#CId.out <- fitIAI_ray(d1=dose1[ind1], e1=fa[ind1], d2=dose2[ind2], e2=fa[ind2], d12=dose1[ind12]+dose2[ind12], e12=fa[ind12], d2.d1, E, alpha=0.05)

    osbaseZdrive <- ifelse(.Platform$OS.type=="windows", "//mymdafiles/usersdqs1", "/home") 
	source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/base.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/moonshot.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/drexplorerExtra/R/interactionIndex.R'))

#library(drexplorerExtra)
library(drexplorer)
data(nl22B2)	
data(UMSCC22B)	
fit_fixedRay <- fitIAI(d1=UMSCC22B[, 1], d2=UMSCC22B[, 2], e=UMSCC22B[, 3], name1='SCH66336', name2='4HPR')
plotIAI(fit_fixedRay, type='IAI', mode='both') # this reproduce the example but not the paper at first; then we shift to CI_IIV2 2008.SSC and now it works!
plotIAI(fit_fixedRay, type='medianEffect', mode='both') # this reproduce the CI paper
plotIAI(fit_fixedRay, type='doseResponseCurve', mode='both')

plotIAI(fit_fixedRay, type='IAI', mode='response') 

res_design <- detect_ray_design(d1=UMSCC22B[, 1], d2=UMSCC22B[, 2], e=UMSCC22B[, 3])


#### this reproduces the example figures exactly in the source code: SYNERGY_V3_original
fit_allPoss <- fitIAI(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1, name1='SCH66336', name2='4HPR')
#medianEffect <- fit_median_efect(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1, d2.d1=res_design$d2.d1)
res_design <- detect_ray_design(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1)
plotIAI(fit_allPoss, type='IAI', mode='both') # this use the CI method; approximately
plotIAI(fit_allPoss, type='contour', mode='both') # this is slightly different from the paper, due to the author's inconsistency
plotIAI(fit_allPoss, type='medianEffect', mode='both')
plotIAI(fit_allPoss, type='doseResponseCurve', mode='both')

fit <- fitIAI(d1=dose1, d2=dose2, e=fa)
plotIAI(fit, type='IAI', mode='both')
x11()
plotIAI(fit, type='IAI', mode='dose')
x11()
plotIAI(fit, type='IAI', mode='response')
plotIAI(fit, type='medianEffect')
plotIAI(fit, type='doseResponseCurve')
plotIAI(fit, type='contour')
}