# drexplorer

[![Build Status](https://travis-ci.org/nickytong/drexplorer.svg?branch=master)](https://travis-ci.org/nickytong/drexplorer)

The **drexplorer** R package is developed to facilitate the analysis of dose-response data. It can be used to:
* assess the reproducibility of replicates, 
* detect outlier data points
* fit different models
* select the best model
* estimate IC values at different percentiles
* evaluate drug-drug interactions using interaction index
* implement NCI 60 method (estimation of GI50, TGI, LD50) ![new](https://www.google.com/url?sa=i&rct=j&q=&esrc=s&source=images&cd=&cad=rja&uact=8&ved=0CAcQjRw&url=http%3A%2F%2Fwww.ebay.it%2Fitm%2FNEW-MEMORY-128MB-SMARTMEDIA-CARD-Memoria-Smart-Media-KORG-Music-Olympus-Agfa-%2F140961939560&ei=rl69VMX3DpKTyATbp4JA&psig=AFQjCNHG4nfz1dAzXV6rXV6-3bSlOd380Q&ust=1421783066203410)
* implement Hill equation (EC50, Emax, Einf, E0 and Hill slope)

## Installation

The [**devtools** package](http://cran.r-project.org/web/packages/devtools/index.html) is used to install R packages hosted on Github. To install **drexplorer**, type the following commands in the R console:

```r
    library(devtools)
    install_github("nickytong/drexplorer")
```

## Usage
```r
# load the package
library(drexplorer)
	
# pull out vignette: vignette does not pull out correctly if installed on Linux with non-root (R>3.0)
# instead, check vignette at: http://htmlpreview.github.com/?https://github.com/nickytong/drexplorer/blob/master/inst/doc/drexplorer.html
vignette('drexplorer')
	
# GUI for dose-response curve fitting
drexplorerGUI_1()
	
# GUI for drug-drug interaction
drexplorerGUI_2()
```    

## User Mannual
A [mannual](http://htmlpreview.github.com/?https://github.com/nickytong/drexplorer/blob/master/inst/doc/drexplorer.html) is available. This is the same document as shown by typing the following command in R console.

```r
vignette('drexplorer')
```
## Cite Us
The **drexplorer** package has been published in Bioinformatics which can be cited as:

Tong, Pan, Coombes, Kevin R., Johnson, Faye M., Byers, Lauren A., Diao, Lixia, Liu, Diane D., . . . Wang, Jing. (2015). drexplorer: A tool to explore dose–response relationships and drug–drug interactions. Bioinformatics. 


## Build Status
[Travis
CI](https://github.com/craigcitro/r-travis) is used to check and build the latest version of **drexplorer**, similar to what CRAN does. The current status icon can be found under the title of this page.  
