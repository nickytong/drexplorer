if(FALSE){
# cd ~/Backup/GitHub/; R --vanilla
library(roxygen2)
library(devtools)

roxygenise("drexplorer")
build('drexplorer')
install('drexplorer')

##
detach("package:drexplorer", unload=TRUE)
library(drexplorer)

build_win('drexplorer')
load_all('drexplorer')

##
library(drexplorer)
#res <- drexplorerGUI_2()
res <- drexplorerGUI_1()
}

# common code in fgui
parse_alpha <- function(){
	alpha <- guiGetValue("alpha")
	if(alpha=='Ignore') {
		alpha <- 1
	} else {
		alpha <- as.numeric(alpha)
	}	
	guiSet("alpha_", alpha) # alpha is modified; so save it as a global var
}
# setup output dir
set_dirOutput_ <- function(){
  # set output dir according to input data
  dirOutput <- file.path(getwd(), 'drexplorer_results')
  if(!file.exists(dirOutput)){
  dir.create(dirOutput, recursive=TRUE)
  }
  guiSet("dirOutput_", dirOutput)
  dirOutput
}


GUI_1_main <- function(datFile, exampleDat, models, freqModels, alpha, genDRfigure, outputRes){
	parse_alpha()
	dat <- guiGetSafe("PERSONAL_dataset")
	if(class(dat)[1] != "data.frame") stop("Data must be loaded.")
	models <- guiGetValue("models")
	alpha <- guiGetSafe("alpha_") # already passed in genDRfigure_press 
	resL <- fitOneExp(dat, drug='', cellLine='', tag='', models=models, alpha=alpha, plot=TRUE, transparency=0.95)
	#lmgui2Callback("genDRfigure") # not updating
	#resL$ICmat
	sprintf('%d outliers detected; The best model is: %s, corresponding IC50=%.3f', sum(resL$datWithOutlierStatus$isOutlier), 
		resL$ICx[1, 'Model'], resL$ICx['IC50'])
}
GUI_1_main_callback <- function(arg){
  if( arg == "datFile" ) {
    datFile_press()
  }else if( arg == "exampleDat" ) {
    exampleDat_press()
  }else if( arg == "freqModels" ) {
    freqModels_press()
  }else if( arg == "genDRfigure" ) {
    genDRfigure_press()
  }else if( arg == "outputRes" ) {
    outputRes_press()
  }
}

freqModels_press <- function(){
	recModels <- c("sigEmax", "LL.4", "LL.5", "linear", "LL.3", "LL.3u", "logistic")
	setListElements("models", recModels)
	#guiSet("models_used", recModels)
	#browser()
	#setListElements("models", recModels)
	#guiGetValue("models") # not working
	#guiGetSafe("models") # this is working
	#getSelectedListElements("models")
}

outputRes_press <- function(){
	parse_alpha()
	drug <- ''
	cellLine <- ''
	dirOutput <- set_dirOutput_()
	pdf(file=file.path(dirOutput, sprintf('%s_%s_DoseResponseCurve.pdf', drug, cellLine)))
	resL <- fitOneExp(dat=guiGetSafe("PERSONAL_dataset"), drug='', cellLine='', tag='', models=guiGetValue("models"), alpha=guiGetSafe("alpha_"), plot=TRUE, transparency=0.95)
	dev.off()
	ICmat <- resL$ICmat
	datWithOutlierStatus <- resL$datWithOutlierStatus
	write.csv(ICmat, file=file.path(dirOutput, sprintf('%s_%s_ICtable.csv', drug, cellLine)))
	write.csv(datWithOutlierStatus, file=file.path(dirOutput, sprintf('%s_%s_datWithOutlierStatus.csv', drug, cellLine)))
	
}


genDRfigure_press <- function(){
	#browser()
	parse_alpha()
	resL <- fitOneExp(dat=guiGetSafe("PERSONAL_dataset"), drug='', cellLine='', tag='', models=guiGetValue("models"), alpha=guiGetSafe("alpha_"), plot=TRUE, transparency=0.95)
}


datFile_press <- function()
{
  f_name <- guiGetValue("datFile")
  f_ext <- tolower(tools::file_ext(f_name))
  if(f_ext=='csv'){
	data <- read.csv(f_name)
  } else {
	data <- read.delim(f_name)
  }
  guiSet("PERSONAL_dataset", data)
  tt <- set_dirOutput_()
}

exampleDat_press <- function() {
  # set output dir according to input data
  dirOutput <- set_dirOutput_()
  # LazyData true #data(ryegrass)
  dose <- ryegrass[, 2]
  response <- ryegrass[, 1]
  ryegrass_dat <- data.frame(dose=dose, response=response)
  ff <- file.path(dirOutput, 'ryegrass.csv')
  write.csv(ryegrass_dat, file=ff, row.names=F)
  guiSetValue("datFile", ff)
  GUI_1_main_callback("datFile")
  
}

#library(fgui)
#' Launch Graphical User Interface (GUI) for dose response curve fitting
#' 
#' This function will launch GUI for curve fitting. The GUI works across different platforms, but the appearance would be slightly different.
#'
#' @export
drexplorerGUI_1 <- function() {
  gui(GUI_1_main,
       argFilename=list(datFile=NULL),
       argCommand = list(exampleDat = NULL, freqModels=NULL, genDRfigure=NULL, outputRes=NULL),
	   #argList =list(models = c("LL.3", "LL.3u", "sigEmax", "logistic")),
	   argList =list(models = c(model.drc, model.DoseFinding)),
	   argOption=list(alpha=c("0.05","0.01","Ignore")), 
       argGridOrder = c(1, 1, 2, 2, 3, 4, 5, 6),
	   callback = GUI_1_main_callback,
	   helps=list(datFile="A tab-delimilated file or csv file with column 1 being the drug concentration and column 2 being the cell viability",
				exampleDat="Use example data from drexplorer package",
				models="A drop-down list of models to be fitted",
				alpha="Significance level in detecting outlier data points. Choose Ignore to omit outlier removal"),
	   title="drexplorer GUI-Curve Fitting",
       argText=list(datFile="Load File ...", exampleDat="Use example data", models="Select models", 
		freqModels="Restrict to frequently used models", alpha="Significance level of outlier detection",
		genDRfigure="Generate dose-response curve", outputRes="Save results"),
       verbose=FALSE )
}
#res <- drexplorerGUI_1()

if(FALSE){
dat <- read.csv('drexplorer/R/ryegrass.csv')
#resL <- fitOneExp(dat, drug='Drug X', cellLine='Cell line A', models=c("LL.3", "LL.3u", "sigEmax", "logistic"), plot=TRUE, transparency=0.95)
#resL <- fitOneExp(dat, drug='Drug X', cellLine='Cell line A', main='', ind2plot='best', cols='blue', models=c("LL.3", "LL.3u", "sigEmax", "logistic", "linlog"), plot=TRUE, transparency=0.95)
resL <- fitOneExp(dat, drug='Drug X', cellLine='Cell line A', tag='', models=c("LL.3", "LL.3u", "sigEmax", "logistic"), plot=TRUE, transparency=0.95)
}


