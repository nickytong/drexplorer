drexplorer
==========

To run the app in your local R session:

    #install app
    library(devtools)
    install_github("drexplorer", "nickytong")
    
    #load in opencpu
    library(opencpu)
    opencpu$browse("/library/drexplorer/www")

	# run with a public server from openCPU
	browseURL("https://public.opencpu.org/ocpu/github/nickytong/drexplorer/www/")
