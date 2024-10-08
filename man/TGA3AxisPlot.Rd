% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TGA functions.R
\name{TGA3AxisPlot}
\alias{TGA3AxisPlot}
\title{Three axis plot maker}
\usage{
TGA3AxisPlot(
  dataframe,
  title = "Placeholder title",
  colorpalette = "Dark2",
  x = "Time (min)",
  y = list("Temperature (°C)", "Weight (mg)", "Heat Flow (mW)")
)
}
\arguments{
\item{dataframe}{Dataframe containing data.}

\item{title}{String to be used as plot title}

\item{colorpalette}{Palette of colors to be used for the first, second, and third lines in the plot. See ?palette for examples. Your own palette can be defined using c("color1","color2","color3")}

\item{x}{String name of column in the dataframe containing the independent variable.}

\item{y}{Vector of strings containing the names of columns to be used as the 3 dependent variables.}
}
\description{
This function generates a plot with three scatterplots of data sharing an independent variable. These plots are commonly made by TA Universal Analysis but may be more generally applicable.
}
