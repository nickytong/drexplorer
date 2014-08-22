if(FALSE){
# cd /data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/


#library(roxygen2)
#library(roxygen) # not working
#roxygenize("drexplorer")

source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/drexplorer.R'))
source(file.path('/data/bioinfo2/ptong1/Projects/Coombes/IC50Package/Package/drexplorer/R/drexplorerAdded.R'))

library(devtools)
build('drexplorer')
install('drexplorer')

detach("package:drexplorer", unload=TRUE)
library(drexplorer)


load_all('drexplorer')
}


### this function is designed to split the data as in this project. It does not work on the drexplorer template data
### major limit: only allows split by one column: thus it only works for data (3 columns: Cellline, Dose, OD) with one drug and multiple cell lines. 
#' Split a drug-response data for the template from Suk Young
#' 
#' This split the data into a list so that each element is a drug-Cell Line combination for further analysis
#' 
#' @param dat input data with multiple cell lines for one drug
#' @param colNames colNames to name selected columns
#' @param by string for the colName used for split
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
#' read data; if extension is csv, treat it as csv; otherwise, treat it as tab delimilated file through read.delim
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
