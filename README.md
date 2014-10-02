# drexplorer

[![Build Status](https://travis-ci.org/nickytong/drexplorer.svg?branch=master)](https://travis-ci.org/nickytong/drexplorer)

The **drexplorer** R package is developed to facilitate the analysis of dose-response data. It can be used to:
* assess the reproducibility of replicates, 
* detect outlier data points
* fit different models
* select the best model
* estimate IC values at different percentiles
* evaluate drug-drug interactions using interaction index

## Installation

The **devtools** package is used to install R packages hosted on Github. To install **drexplorer**, type the following commands in the R console:

```r
    library(devtools)
    install_github("drexplorer", "nickytong")
```

## Usage
```r
# load the package
library(drexplorer)
	
# pull out vignette
vignette('drexplorer')
	
# GUI for dose-response curve fitting
drexplorerGUI_1()
	
# GUI for drug-drug interaction
drexplorerGUI_2()
```    

## User Mannual
A [mannual](https://github.com/nickytong/drexplorer/blob/master/vignettes/drexplorer.pdf) is available. This is the same document as shown by typing the following command in R console.

```r
vignette('drexplorer')
```

## Build Status
[Travis
CI](http://yihui.name/en/2013/04/travis-ci-general-purpose/) is used to check and build the latest version of **drexplorer**, similar to what CRAN does. The current status icon can be found on the top of this page.  
