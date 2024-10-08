% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPMSfunctions.R
\name{PlotWithInset}
\alias{PlotWithInset}
\title{Inset plot generator}
\usage{
PlotWithInset(
  dataframe,
  xcolumn = "TemperatureK",
  ycolumn = "susc",
  xstart = 0,
  xend = 300,
  maintitle = "",
  subtitle = "",
  xtitle = "Temperature (K)",
  ytitle = "Magnetization (emu)",
  color = "steelblue"
)
}
\arguments{
\item{dataframe}{Name of dataframe containing columns of data.}

\item{xcolumn}{String name of the column in the dataframe to be your independent variable.}

\item{ycolumn}{String name of the column in the dataframe to be your dependent variable.}

\item{xstart}{First value in your independent variable list to plot. This is the actual number for the x-axis minimum, not an index.}

\item{xend}{Last value in your independent variable list to plot. This is the actual number for the x-axis maximum, not an index.}

\item{maintitle}{String to be the title of your main plot.}

\item{subtitle}{String to be the title of your inset plot.}

\item{xtitle}{String to be the title of your x-axis.}

\item{ytitle}{String to be the title of your y-axis.}

\item{color}{Color to plot the data in.}
}
\description{
This function makes a plot with an inset. The specifics are tailored to viewing magnetic susceptibility data, but it may be more generally applicable.
}
