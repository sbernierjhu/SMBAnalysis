% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HC functions.R
\name{CvCpOverlayPlot}
\alias{CvCpOverlayPlot}
\title{Heat Capacity Plotter - overlaid Cp and Cv}
\usage{
CvCpOverlayPlot(
  CvTemps,
  CvVals,
  CvCol = "#1D3732",
  CpTemps,
  CpVals,
  CpCol = "#397F72",
  CpCoeff = 1,
  xmax = 300,
  maintitle = "Heat capacity overlay",
  perT3 = TRUE
)
}
\arguments{
\item{CvTemps}{Column in a dataframe containing temperature data for Cv plot in Kelvin, as in dataframe$column}

\item{CvVals}{Column in a dataframe containing Cv data in J/K, as in dataframe$column}

\item{CvCol}{Color for Cv data and axis}

\item{CpTemps}{Column in a dataframe containing temperature data for Cp plot in Kelvin, as in dataframe$column}

\item{CpVals}{Column in a dataframe containing Cp data in J/K, as in dataframe$column}

\item{CpCol}{Color for Cp data and axis}

\item{CpCoeff}{Scaling coefficient for the Cp data - guess and check until your data shows up with appropriate magnitude}

\item{maintitle}{String to be main title of your plot.}

\item{perT3}{Boolean - when TRUE, units are given for C/T^3. When FALSE, units are given for C.}
}
\description{
This function makes a pretty picture of experimental AND calculated heat capacity data overlaid with color-coded axes.
}
