
## ----, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(fig.align='center', fig.show='asis', 
		dev='png', dpi=72*2, fig.width=5, fig.height=5)


## ----, message = FALSE---------------------------------------------------
library(drexplorer)


## ------------------------------------------------------------------------
#library(drc)
#data(ryegrass, package='drc')
dose <- ryegrass[, 2]
response <- ryegrass[, 1]


## ------------------------------------------------------------------------
## potential outlier at dose 3.75
NewmanTest(ref=response[dose==0], obs=response[dose==3.75], alpha=0.05)


## ------------------------------------------------------------------------
drOutlier(drMat=ryegrass[, c(2, 1)], alpha=0.05)


## ------------------------------------------------------------------------
set.seed(100)
r1 <- runif(28)
r2 <- r1+rnorm(28, 0, 0.1)
ccc <- plotCCC(r1, r2, xlab='Simulated response, replicate 1', ylab='Simulated response, replicate 2')
ccc


## ------------------------------------------------------------------------
fit_sigEmax_alpha1 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", 
	alpha=1, fitCtr=FALSE)


## ------------------------------------------------------------------------
fit_sigEmax_alpha_o5 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", 
	alpha=0.05, fitCtr=FALSE)
fit_sigEmax_alpha1@fit
fit_sigEmax_alpha_o5@fit


## ------------------------------------------------------------------------
y <- predict(fit_sigEmax_alpha_o5)
y


## ------------------------------------------------------------------------
computeIC(fit_sigEmax_alpha_o5, percent=seq(0, 1, by=0.1), log.d=FALSE, interpolation=TRUE)
computeIC(fit_sigEmax_alpha_o5, percent=seq(0, 1, by=0.1), log.d=FALSE, interpolation=FALSE)


## ------------------------------------------------------------------------
fit.LL.3 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3", alpha=0.05, fitCtr=FALSE)
fit.LL.3u <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "LL.3u", alpha=0.05, fitCtr=FALSE)
fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.05, fitCtr=FALSE)
###
plot(fit.LL.3, main='', col=4, lwd=2)
lines(fit.LL.3u, col=5, lwd=2)
lines(fit.sigEmax, col=6, lwd=2)
legend("bottomleft", c('LL.3', 'LL.3u', 'sigEmax'), col=4:6, lwd=3)


## ------------------------------------------------------------------------
sapply(list(fit.LL.3, fit.LL.3u, fit.sigEmax), function(x) x@info$RSE)


## ------------------------------------------------------------------------
# no outlier excluded
fit.sigEmax0 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=1, fitCtr=FALSE)
###
plot(fit.sigEmax0, main='sigEmax model', col=7, lwd=2)
lines(fit.sigEmax, col=6, lwd=2)
legend("bottomleft", c('alpha=0.05', 'ignored'), col=c(6, 7), lwd=3)


## ----, eval=TRUE---------------------------------------------------------
fitExp <- fitOneExp(ryegrass[, c(2, 1)], drug='', cellLine='', unit='', models=c('sigEmax', 'LL.4', 'LL.5', 'LL.3', 'LL.3u', 'logistic'), alpha=0.05, interpolation=TRUE)


## ----, fig.width=7, fig.height=7, eval=TRUE------------------------------
plotOneExp(fitExp, main='')



## ------------------------------------------------------------------------
sdat <- prepDRdat(drMat=ryegrass[, c(2, 1)], alpha=0.01, fitCtr=FALSE, standardize=TRUE)$dat
fit_hill <- hillFit(d=sdat$dose, y=sdat$response)
plot(fit_hill)
fit_hill[1:length(fit_hill)]


## ------------------------------------------------------------------------
fit_nci60 <- nci60Fit(d=datNCI60$Dose, y=datNCI60$Growth/100)
fit_nci60[1:length(fit_nci60)]
plot(fit_nci60)


## ------------------------------------------------------------------------
#data(UMSCC22B)	
fit_fixedRay <- fitIAI(d1=UMSCC22B[, 1], d2=UMSCC22B[, 2], 
			e=UMSCC22B[, 3], name1='SCH66336', name2='4HPR')


## ------------------------------------------------------------------------
# IAI vs response
plotIAI(fit_fixedRay, type='IAI', mode='response') 


## ------------------------------------------------------------------------
# IAI versus dose
plotIAI(fit_fixedRay, type='IAI', mode='dose') 


## ------------------------------------------------------------------------
# median effect
plotIAI(fit_fixedRay, type='medianEffect') 


## ------------------------------------------------------------------------
#data(nl22B2)	
fit_allPoss_1 <- fitIAI(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1, name1='SCH66336', name2='4HPR',d2.d1.force=1)
fit_allPoss_2 <- fitIAI(d1=nl22B2$schd, d2=nl22B2$hpr, e=nl22B2$y1, name1='SCH66336', name2='4HPR',d2.d1.force=2)


## ------------------------------------------------------------------------
# median effect
plotIAI(fit_allPoss_1, type='medianEffect') 


## ------------------------------------------------------------------------
plotIAI(fit_allPoss_2, type='medianEffect') 


## ------------------------------------------------------------------------
plotCCC(fit_allPoss_1$CI$IAI, fit_allPoss_2$CI$IAI)


## ----sessInfo, results='asis', echo=FALSE--------------------------------
sessionInfo()


