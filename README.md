drexplorer
==========

The drexplorer R package is developed to facilitate the analysis of dose-response data. It can be used to:
* assess the reproducibility of replicates, 
* detect outlier data points
* fit different models
* select the best model
* estimate IC values at different percentiles
* evaluate drug-drug interactions

Instructions for installation and usage:

    #install package
    library(devtools)
    install_github("drexplorer", "nickytong")
    
	# load the package
	library(drexplorer)
	
	# pull out vignette
	vignette('drexplorer')

	# GUI for dose-response curve fitting
	drexplorerGUI_1()
	
	# GUI for drug-drug interaction
	drexplorerGUI_2()
