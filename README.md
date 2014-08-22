drexplorer
==========

To run the app in your local R session:

    #install app
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