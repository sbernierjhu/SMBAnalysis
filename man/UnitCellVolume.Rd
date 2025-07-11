% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cubicity metric.R
\name{UnitCellVolume}
\alias{UnitCellVolume}
\title{Unit cell volume calculator}
\usage{
UnitCellVolume(a, b, c, al, be, ga, sigfigs = 4)
}
\arguments{
\item{a}{a lattice parameter, typically in angstroms}

\item{b}{b lattice parameter, typically in angstroms}

\item{c}{c lattice parameter, typically in angstroms}

\item{al}{alpha lattice parameter (angle between b and c sides) in degrees}

\item{be}{beta lattice parameter (angle between a and c sides) in degrees}

\item{ga}{gamma lattice parameter (angle between a and b sides) in degrees}

\item{sigfigs}{Number of sig figs to round answer to. Default is 4.}
}
\value{
Volume of the unit cell in (units of a/b/c)^3
}
\description{
This function calculates the volume of the unit cell
}
\examples{
library(tidyverse)
UnitCellVolume(10,10,10,90,90,90,decimals=3)
UnitCellVolume(10,8,9,63,70,63,decimals=3)
df <- data.frame(a=c(10,10),b=c(10,8),c=c(10,9),al=c(90,63),be=c(90,70),ga=c(90,63))
df <- df \%>\% mutate(Volume=UnitCellVolume(a,b,c,al,be,ga,4))
head(df)
}
