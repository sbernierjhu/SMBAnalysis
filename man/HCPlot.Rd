% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HC functions.R
\name{HCPlot}
\alias{HCPlot}
\title{Heat Capacity Plotter - single dataset}
\usage{
HCPlot(
  dataframe,
  tempvals = Temperature,
  yvals = Cv,
  title = "Heat capacity",
  isCv = TRUE,
  perT3 = TRUE,
  ptcolor = "darkgreen"
)
}
\arguments{
\item{dataframe}{Dataframe containing the data you want to plot.}

\item{tempvals}{Name of the column containing temperature data in Kelvin.}

\item{yvals}{Name of the column containing Cp or Cv data in J/K.}

\item{title}{String to be main title of your plot.}

\item{isCv}{Boolean - when TRUE, units come out for constant volume heat capacity. When FALSE, units come out for constant pressure.}

\item{perT3}{Boolean - when TRUE, units are given for C/T^3. When FALSE, units are given for C.}

\item{ptcolor}{Color for the plot}
}
\description{
This function makes a pretty picture of heat capacity data with labelled axes.
}
