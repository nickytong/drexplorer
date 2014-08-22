\name{plotIAI}
\alias{plotIAI}
\title{Visualize interaction index}
\usage{
  plotIAI(fit,
    type = c("IAI", "medianEffect", "doseResponseCurve", "contour"),
    contour.level = (1:9)/10, mode = "both")
}
\arguments{
  \item{fit}{the fitted result from fitIAI()}

  \item{type}{type of plot from c('IAI', 'medianEffect',
  'doseResponseCurve', 'contour')}

  \item{contour.level}{contour level. only effective if
  type='contour'}

  \item{mode}{specify if to plot against response, dose or
  both. only effective if type=='IAI'. can be either
  'response', 'dose', or 'both'}
}
\description{
  Visualize interaction index
}

