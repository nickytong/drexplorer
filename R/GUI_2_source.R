# setup output dir
set_dirOutput2_ <- function(){
  # set output dir according to input data
  dirOutput <- file.path(getwd(), 'drexplorer_results')
  if(!file.exists(dirOutput)){
  dir.create(dirOutput, recursive=TRUE)
  }
  guiSet("dirOutput_", dirOutput)
  dirOutput
}


GUI_2_main <- function(datFile, exampleDat_UMSCC22B, exampleDat_nl22B2, IAI, IAImode, medianEffect, doseResponseCurve, contour, outputRes){
	dd <- guiGetSafe("PERSONAL_dataset_2")
	if(class(dd)[1] != "data.frame") stop("Data must be loaded.")
	res_design <- detect_ray_design(d1=dd[, 1], d2=dd[, 2], e=dd[, 3])
	#browser()
	## get aux text
	text_design <- ifelse(res_design$isRayDesign==TRUE, sprintf('This is a fixed ratio design.'), sprintf('This is NOT a fixed ratio design.\t'))
	txDim <- sprintf('%d rows; %d columns in the data; %d effective records (0~1)', nrow(res_design$drMat), ncol(res_design$drMat), sum(!res_design$ind_eAbnorm))
	txRcr <- sprintf('drug 1 only has %d records; drug 2 only has %d records, mixture has %d records.', length(res_design$d1), length(res_design$d2), length(res_design$d12))
	tx_d2d1 <- sprintf('Observed d2.d1: %s; d2.d1 to be used: %.2f', str_c(res_design$d2.d1_avail, collapse=', '), res_design$d2.d1)
	texttt <- sprintf('>>>>>%s\n-------->%s\n%s\n%s\n\n', txDim, text_design, txRcr, tx_d2d1)
	#browser()
	#lmgui2Callback("genDRfigure") # not updating
	#resL$ICmat	
	#cat(texttt) # wrong!
	texttt
}
GUI_2_main_callback <- function(arg){
  if( arg == "datFile" ) {
    datFile_press2()
  }else if( arg == "exampleDat_UMSCC22B" ) {
    exampleDat_UMSCC22B_press()
  }else if( arg == "exampleDat_nl22B2" ) {
    exampleDat_nl22B2_press()
  }else if( arg == "IAI" ) {
    IAI_press()
  }else if( arg == "medianEffect" ) {
    medianEffect_press()
  }else if( arg == "doseResponseCurve" ) {
    doseResponseCurve_press()
  }else if( arg == "contour" ) {
    contour_press()
  }else if( arg == "outputRes" ) {
    outputRes_press2()
  }
}


outputRes_press2 <- function(){
	dirOutput <- set_dirOutput2_()
	pdf(file=file.path(dirOutput, sprintf('IAIplot_for_%s.pdf', guiGetSafe("FILE_NAME"))))
	dd <- guiGetSafe("PERSONAL_dataset_2")
	fit <- fitIAI(d1=dd[, 1], d2=dd[, 2], e=dd[, 3], name1=colnames(dd)[1], name2=colnames(dd)[2])
	plotIAI(fit,  type='medianEffect')
	plotIAI(fit,  type='doseResponseCurve')
	plotIAI(fit,  type='contour')
	plotIAI(fit,  type='IAI', mode='both')
	plotIAI(fit,  type='IAI', mode='response')
	plotIAI(fit,  type='IAI', mode='dose')
	dev.off()
	write.csv(fit$CI, file=file.path(dirOutput, sprintf('IAI_CI_%s.csv', guiGetSafe("FILE_NAME"))))
	
}

## fig generation engine
fig_gen <- function(type='IAI', mode='both', return=FALSE){
	dd <- guiGetSafe("PERSONAL_dataset_2")
	fit <- fitIAI(d1=dd[, 1], d2=dd[, 2], e=dd[, 3], name1=colnames(dd)[1], name2=colnames(dd)[2])
	plotIAI(fit, type=type, mode=mode)
	if(return){
		return(fit)
	}
}

IAI_press <- function(){
	#browser()
	fig_gen(type='IAI', mode=guiGetValue("IAImode"))
}

medianEffect_press <- function(){
	#browser()
	fig_gen(type='medianEffect')
}

doseResponseCurve_press <- function(){
	#browser()
	fig_gen(type='doseResponseCurve')	
}

contour_press <- function(){
	#browser()
	fit <- fig_gen(type='contour', return=TRUE)
	if(fit$meta$isRayDesign==TRUE){
		tkmessageBox( message="This is a fixed ratio design; No contour plot can be constructed!", title="Warning" )
	}
}

datFile_press2 <- function()
{
  f_name <- guiGetValue("datFile")
  f_ext <- tolower(tools::file_ext(f_name))
  baseName <- str_split(basename(f_name), '\\.')[[1]][1]
  guiSet("FILE_NAME", baseName)
  if(f_ext=='csv'){
	data <- read.csv(f_name)
  } else {
	data <- read.delim(f_name)
  }
  guiSet("PERSONAL_dataset_2", data)
  tt <- set_dirOutput2_()
}

exampleDat_UMSCC22B_press <- function() {
  # set output dir according to input data
  dirOutput <- set_dirOutput2_()
  data(UMSCC22B)
  ff <- file.path(dirOutput, 'UMSCC22B.csv')
  write.csv(UMSCC22B, file=ff, row.names=F)
  guiSetValue("datFile", ff)
  GUI_2_main_callback("datFile")  
}

exampleDat_nl22B2_press <- function() {
  # set output dir according to input data
  dirOutput <- set_dirOutput2_()
  data(nl22B2)
  ff <- file.path(dirOutput, 'nl22B2.csv')
  write.csv(nl22B2, file=ff, row.names=F)
  guiSetValue("datFile", ff)
  GUI_2_main_callback("datFile")  
}

#library(fgui)
#' Launch Graphical User Interface (GUI) for drug interaction (IAI) analysis
#' 
#' This function will launch GUI for computing Interaction Index (IAI) in fixed ratio design. 
#' The GUI works across different platforms, but the appearance would be slightly different.
#'
#' @export
drexplorerGUI_2 <- function() {
  gui(GUI_2_main,
       argFilename=list(datFile=NULL),
       argCommand = list(exampleDat_UMSCC22B=NULL, exampleDat_nl22B2=NULL, 
			IAI=NULL, medianEffect=NULL, doseResponseCurve=NULL, contour=NULL, outputRes=NULL),
	   argList =list(exampleDat = c('UMSCC22B', 'nl22B2')),
	   argOption=list(IAImode=c('both', 'response', 'dose')), 
       #argGridOrder = c(1, 1, 1, 2, 2, 3, 4, 5, 6),
       argGridOrder = c(1, 2, 2, 3, 3, 4, 5, 6, 7),
	   callback = GUI_2_main_callback,
	   helps=list(datFile="A tab-delimilated file or csv file with column 1 being the drug 1 concentration and column 2 being drug 2 concentration and column 3 being the cell viability (between 0 and 1)",
				exampleDat_UMSCC22B="Select the UMSCC22B data as an example",
				exampleDat_nl22B2="Select the nl22B2 data as an example",
				IAI="Generate Interaction Index plot",
				IAImode="The quantity to be plotted against the Interaction Index",
				alpha="Significance level in detecting outlier data points. Choose Ignore to omit outlier removal"),
	   title="drexplorer GUI--Drug Interaction",
       argText=list(datFile="Load File ...              ", exampleDat_UMSCC22B="Use UMSCC22B data", exampleDat_nl22B2="Use nl22B2 data", 
		IAI="Generate IAI plot", 
		IAImode="Quantity against IAI", medianEffect="Generate Median-Effect plot",
		doseResponseCurve="Generate dose-response curve", contour="Generate contour plot", 
		outputRes="Save results"),
       verbose=FALSE )
}
#res <- drexplorerGUI_2()

if(FALSE){
dat <- read.csv('drexplorer/R/ryegrass.csv')
#resL <- fitOneExp(dat, drug='Drug X', cellLine='Cell line A', models=c("LL.3", "LL.3u", "sigEmax", "logistic"), plot=TRUE, transparency=0.95)
#resL <- fitOneExp(dat, drug='Drug X', cellLine='Cell line A', main='', ind2plot='best', cols='blue', models=c("LL.3", "LL.3u", "sigEmax", "logistic", "linlog"), plot=TRUE, transparency=0.95)
resL <- fitOneExp(dat, drug='Drug X', cellLine='Cell line A', tag='', models=c("LL.3", "LL.3u", "sigEmax", "logistic"), plot=TRUE, transparency=0.95)
}


