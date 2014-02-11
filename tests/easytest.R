library(drexplorer)
data(ryegrass)
#(1) identify outliers
dose <- ryegrass[, 2]
response <- ryegrass[, 1]
NewmanTest(ref=response[dose==0], obs=response[dose==3.75], alpha=0.05)
#
drOutlier(drMat=ryegrass[, c(2, 1)], alpha=0.05)

#(2) fit 
fit.sigEmax <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=0.05, fitCtr=FALSE)

#(3) predict
predict(fit.sigEmax)


#(4) plot functions
fit.sigEmax0 <- drFit(drMat=ryegrass[, c(2, 1)], modelName = "sigEmax", alpha=1, fitCtr=FALSE)
###
plot(fit.sigEmax0, main='sigEmax model', col=7, lwd=2)
lines(fit.sigEmax, col=6, lwd=2)
legend("bottomleft", c('alpha=0.05', 'ignored'), col=c(6, 7), lwd=3)
