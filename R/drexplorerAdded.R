#' Calculates Concordance Correlation Coefficient (CCC) to access reproducibility
#' 
#' This function calculates and plots concordance for paired data.
#'
#'
#' infinite values in x and y are masked as NA so as to compute CCC and corr (otherwise, nan).
#' @param x x vector
#' @param y y vector 
#' @param xlab xlab
#' @param ylab ylab
#' @param maskBeyondLim whether mask values beyond xlim (for x) or ylim (for y) when xlim and ylim is specified. Default is FALSE so that even xlim ylim specified, CCC will not be affected
#' @param tag add a tag
#' @param col_legend color for legend text
#' @param pch pch of the points
#' @param main main title 
#' @param cex for the dots 
#' @param plot whether draw a figure
#' @param ylim y axis limit 
#' @param xlim x axis limit 
#' @param plotOutlier logical whether add an outlier plot comparing Diff vs Mean plot and related statistics
#' @param sampleName symbols of the text to be shown in outlier plot
#' @param alpha_outlier alpha to call outliers based on the observed differences assuming from Normal(mean_diff, sd_diff)
#' @return a vector of c('ccc', 's_shift', 'l_shift', 'ccc_lo', 'ccc_hi', 'Cb', 'corr'). s_shift is scale shift; l_shift is location shift;
#'			ccc_lo and ccc_hi represent 95% confidence interval for CCC. Cb for the bias correction term satisfying CCC=corr*Cb where corr is the
#' 			Pearson correlation
#' @export
plotCCC <- function(x, y, xlab=NA, ylab=NA, tag='', col_legend='red', pch=1, main=NA, cex=2, plot=TRUE, ylim=NA, xlim=NA, maskBeyondLim=FALSE, 
	plotOutlier=FALSE, sampleName=NA, alpha_outlier=0.01){
	#require(epiR)
	#require(calibrate) # text around points in base graphics
	if(is.na(xlab)) xlab <- deparse(substitute(x)) # should happen at very beginning
	if(is.na(ylab)) ylab <- deparse(substitute(y))
	#browser()
	if(!is.vector(x)) x <- as.vector(x)
	if(!is.vector(y)) y <- as.vector(y)
	if(length(x)!=length(y)) stop('Length of x and y is unequal!')
	if(!is.na(ylim[1]) & maskBeyondLim) y <- maskByLimit(y, Lower=ylim[1], Upper=ylim[2])
	if(!is.na(xlim[1]) & maskBeyondLim) x <- maskByLimit(x, Lower=xlim[1], Upper=xlim[2])
	## inf value will lead to NaN in CCC and cor: for simplicity, mask as NA
	isInf_x <- is.infinite(x)
	isInf_y <- is.infinite(y)
	if(sum(isInf_x)>0) warning(sprintf('%d infinite values observed in x, masked as NA in order to compute CCC and Cor!', sum(isInf_x)))
	if(sum(isInf_y)>0) warning(sprintf('%d infinite values observed in y, masked as NA in order to compute CCC and Cor!', sum(isInf_y)))
	x[isInf_x] <- NA
	y[isInf_y] <- NA
	## remove NA: so that cor and CCC will not be NA
	tdf <- noNA(data.frame(x=x, y=y))
	x <- tdf$x
	y <- tdf$y
	#browser()
	## the default epi.ccc use sd(y)/sd(x) and y-x as the scale and location shift!
	ccc_fit <- epi.ccc(x, y, ci = "z-transform", conf.level = 0.95) 
	############### outlier calculation
	dat_outlier <- data.frame(avg=(x+y)*0.5, dif=y-x)
	### so inconsistent!: this is x-y!; use my own calc!
	#diff_mean <- mean(ccc_fit$blalt$delta)
	#diff_sd <- sqrt(var(ccc_fit$blalt$delta))
	diff_mean <- mean(dat_outlier$dif, na.rm=TRUE)
	diff_sd <- sd(dat_outlier$dif, na.rm=TRUE)
	### basic statistics that may be used for filtering
	mean_x <- mean(x, na.rm=TRUE)
	mean_y <- mean(y, na.rm=TRUE)
	sd_x <- sd(x, na.rm=TRUE)
	sd_y <- sd(y, na.rm=TRUE)
	#browser()
	extremes <- qnorm(c(alpha_outlier/2, 1-alpha_outlier/2), diff_mean, diff_sd) # lower and upper bound for outliers
	dat_outlier$isOutlier <- dat_outlier[, 'dif'] < extremes[1] | dat_outlier[, 'dif'] > extremes[2] 
	if(is.na(sampleName[1])) sampleName <- 1:nrow(dat_outlier)
	dat_outlier$sampleName <- sampleName
	#diff_mean-2*diff_sd
	#diff_mean+2*diff_sd
	############### CCC calculation
	tt <- unlist(ccc_fit$rho.c)
	s_shift <- ccc_fit$s.shift # scale shift , that is sigma1/sigma2
	l_shift <- ccc_fit$l.shift # location shift, that is (mu1-mu2)/sqrt(sigma1*sigma2)
	ccc <- tt['est'] # CCC estimate
	ccc_lo <- tt['lower'] # 95% CI lower
	ccc_hi <- tt['upper'] # 95% CI upper
	corr <- cor(x, y, use='complete.obs')
	Cb <- ccc_fit$C.b # bias correction factor; cor*Cb=CCC
	CCCinfo0 <- c(ccc, s_shift, l_shift, ccc_lo, ccc_hi, Cb, corr)
	names(CCCinfo0) <- c('ccc', 's_shift', 'l_shift', 'ccc_lo', 'ccc_hi', 'Cb', 'corr') 
	#res <- list(CCCinfo=res, mean_x=mean_x, mean_y=mean_y, sd_x=sd_x, sd_y=sd_y)
	aux <- c(mean_x=mean_x, mean_y=mean_y, sd_x=sd_x, sd_y=sd_y)
	CCCinfo <- c(CCCinfo0, aux)
	res <- CCCinfo
	#browser()
	if(plot){
		## when x is a mat, this xlab is the element of the mat and really long!
		if(length(xlab)>2) xlab <- 'x'
		if(length(ylab)>2) ylab <- 'y'
		#browser()
		if(is.na(main)) main <- sprintf('Concordance plot: CCC=%.3f, Corr=%.3f', ccc, corr)
		rr <- range(c(x, y), na.rm=TRUE) # using largest range possible so that we have a symmetric view
		if(!is.na(ylim[1])){
			ylim <- ylim
		} else {
			ylim <- rr
		}
		if(!is.na(xlim[1])){
			xlim <- xlim
		} else {
			xlim <- rr
		}
		#browser()
		plot(x, y, xlab=xlab, ylab=ylab, main=str_c(tag, '\n', main, collapse=''), pch=pch, cex=cex, xlim=xlim, ylim=ylim)
		#abline(0, 1, col=80, lty=2, lwd=2)
		abline(a = 0, b = 1, lty = 2, col=4, lwd=3) # identical line
		#browser()
		try(abline(lm(y~x), lty = 1, lwd=2), silent=TRUE) # lm line
		legend('bottomright', legend = c("Perfect concordance", "Observed linear trend"), lty = c(2,1), bty = "n", col=c(4, 1), lwd=c(3, 2), text.col=col_legend)
		labCCC <- sprintf("-->CCC: %.3f, (95%% CI %.2f~%.2f)", res['ccc'], res['ccc_lo'], res['ccc_hi'])
		lab_s_shift <- sprintf("-->Scale shift: %.2f", res['s_shift'])
		lab_l_shift <- sprintf("-->Location shift: %.2f", res['l_shift'])
		legend('topleft', legend = c(labCCC, lab_s_shift, lab_l_shift), bty = "n", text.col=col_legend)
	}
	if(plotOutlier){
		#gtitle <- str_c(tag, sprintf('%d disconcordant outliers at significance level of %.3f', sum(dat_outlier$isOutlier, na.rm=TRUE), alpha_outlier), collapse='')
		#p <- ggplot(dat_outlier, aes(x=avg, y=dif))+geom_point(shape=1, size=6)+
		#	geom_text(data=subset(dat_outlier, isOutlier==TRUE), aes(label=sampleName))+
		#	geom_hline(yintercept = extremes, linetype=2, color='red')+
		#	geom_hline(yintercept = diff_mean, linetype=2, color='black')+
		#	xlab('Average measurement')+ylab('Difference between the two measurements')+ggtitle(gtitle)
		#print(p)
		#browser()
		gtitle <- str_c(tag, sprintf('%d disconcordant outliers at significance level of %.3f', sum(dat_outlier$isOutlier, na.rm=TRUE), alpha_outlier), collapse='')
		ylim <- c(min(c(extremes, dat_outlier$dif), na.rm=TRUE), max(c(extremes, dat_outlier$dif), na.rm=TRUE))
		with(dat_outlier, plot(avg, dif, cex=2, pch=1, ylim=ylim, xlab='Average measurement', ylab='Difference between the two measurements', main=gtitle))
		abline(h=extremes, lty=2, col='red')
		abline(h=diff_mean, lty=2, col='black')
		tt <- (data=subset(dat_outlier, isOutlier==TRUE))
		if(nrow(tt)>0){
			with(tt, textxy(avg, dif, labs=sampleName, cex=1.1, offset=0))
		}
		res <- list(CCCinfo=CCCinfo, dat_outlier=dat_outlier, diff_mean=diff_mean, diff_sd=diff_sd, indexOutlier=which(dat_outlier$isOutlier), indexNonOutlier=which(dat_outlier$isOutlier==FALSE))
	}	
	res
}
#plotCCC(scoreEMT0_raw, scoreEMT0_scaled, tag='74 gene signature')


### this function is designed to split the data as in this project. It does not work on the drexplorer template data
### major limit: only allows split by one column: thus it only works for data (3 columns: Cellline, Dose, OD) with one drug and multiple cell lines. 
#' Split a drug-response data  
#' 
#' This split the data into a list so that each element is a list. This is useful if the data contains multiple cell lines or multiple drugs for the same cell line
#' Currently it only allows splitting by one variable (one column). If data contains both multiple cell lines and multiple drugs, a customized function needs to be written for data formatting.
#' 
#' @param dat input data with multiple cell lines for one drug
#' @param colNames colNames to name selected columns
#' @param by string for the colName (just one column allowed) used for split
#' @param minReplicate if a subset does not have nrow exceeding this, throw an error
#' @param verbose whether to produce warning message
#' @return  a list of the split data frames
#' @export
datSplit <- function(dat, colNames, by, minReplicate=3, verbose=TRUE){
	#browser()
	dat <- dat[, 1:length(colNames)] # in case there is extra column at the end, i.e. X, all NA 
	if(!identical(colnames(dat), colNames) & verbose){
		cat('Warning: Data colname is probably inconsistent!\n')
		cat('-----------------------\n')
		warning(cat(sprintf('Original colnames: \n%s\nNow is forced to: \n%s\n', 
			str_c(colnames(dat), collapse=', '), str_c(colNames, collapse=', '))))
		cat('-----------------------\n')
		colnames(dat) <- colNames #
	}
	datsplit <- dlply(dat, by, function(x) x) # a list, each cell line per element
	#dfsplit_labels <- unname(unlist(attributes(datsplit)$split_labels)) # 
	dfsplit_labels <- names(datsplit) # 
	nMeasure <- sapply(datsplit, nrow) # number of rows for each drug-CL
	if(any(nMeasure<minReplicate)){
		indWrong <- which(nMeasure<minReplicate)
		msg <- sprintf('Cell line \n  %s\nhas too few measurements:\n  %s\n 
			Make sure you have the right data!\n', str_c(dfsplit_labels[indWrong], collapse=', '), str_c(nMeasure[indWrong], collapse=', '))
		cat('Error: Data is probably wrong!\n')
		cat('-----------------------\n')
		stop(cat(msg))
		cat('-----------------------\n')
	}
	datsplit
	#browser()
}

############### wrapper
# read data; if extension is csv, treat it as csv; otherwise, treat it as tab delimilated file through read.delim
loadData <- function(filepath){
	ext <- tolower(tools::file_ext(filepath))
	if(ext=='csv'){
		dat <- read.csv(filepath)
	} else {
		dat <- read.delim(filepath)
	}
	dat
}
############### gui
drexplorer <- function(){
	library(drexplorer)
	
}
