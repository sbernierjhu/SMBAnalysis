
#'Cubic Deviation Metric calculator (scalar; old)
#'
#'This function calculates the value of the cubicity metric described in http://dx.doi.org/10.1038/s41535-025-00782-3. This function is deprecated but preserved for the purposes of interpreting the data in that paper.
#'
#'The returned value of the cubicity metric will be zero if and only if the lattice parameters input correspond to a cubic unit cell. The returned value will be unitless and normalized such that it can be compared to values for other unit cells. The closer to zero, the "more cubic" your unit cell. It is important to note that multiple unit cells may produce the same metric value.
#'
#' @param a a lattice parameter, typically in angstroms
#' @param b b lattice parameter, typically in angstroms
#' @param c c lattice parameter, typically in angstroms
#' @param al alpha lattice parameter (angle between b and c sides) in degrees
#' @param be beta lattice parameter (angle between a and c sides) in degrees
#' @param ga gamma lattice parameter (angle between a and b sides) in degrees
#' @param sigfigs Number of sig figs to round answer to. Default is 4.
#' @returns Value of cubicity metric
CubicityMetric <- function(a,b,c,al,be,ga,sigfigs=4) signif(abs((a/sqrt(a**2+b**2-2*a*b*cos(ga*(pi/180))))-1/sqrt(2))+abs((b/sqrt(a**2+b**2-2*a*b*cos(ga*(pi/180))))-1/sqrt(2))+abs((c/sqrt(a**2+c**2-2*a*c*cos(be*(pi/180))))-1/sqrt(2))+abs((a/sqrt(a**2+c**2-2*a*c*cos(be*(pi/180))))-1/sqrt(2))+abs((b/sqrt(b**2+c**2-2*b*c*cos(al*(pi/180))))-1/sqrt(2))+abs((c/sqrt(b**2+c**2-2*b*c*cos(al*(pi/180))))-1/sqrt(2)),sigfigs)

#'Cubic Deviation Metric calculator (newest)
#'
#'This function calculates the value of the cubicity metric described in xxx.
#'
#'The returned value of the cubicity metric will be zero if and only if the lattice parameters input correspond to a cubic unit cell. The returned value will be unitless and normalized such that it can be compared to values for other unit cells, with a possible maximum of 1. The closer to zero, the "more cubic" your unit cell; the closer to one, the "less cubic". It is important to note that multiple unit cells may produce the same metric value.
#'
#' @param a a lattice parameter, typically in angstroms
#' @param b b lattice parameter, typically in angstroms
#' @param c c lattice parameter, typically in angstroms
#' @param al alpha lattice parameter (angle between b and c sides) in degrees
#' @param be beta lattice parameter (angle between a and c sides) in degrees
#' @param ga gamma lattice parameter (angle between a and b sides) in degrees
#' @param decimals Number of decimal places to round answer to. Default is 4.
#' @examples
#' library(tidyverse)
#' CubicDeviationMetric(10,10,10,90,90,90,decimals=3)
#' CubicDeviationMetric(10,8,9,63,70,63,decimals=3)
#' df <- data.frame(a=c(10,10),b=c(10,8),c=c(10,9),al=c(90,63),be=c(90,70),ga=c(90,63))
#' df <- df %>% mutate(CDM=CubicDeviationMetric(a,b,c,al,be,ga,4))
#' head(df)
#' @returns Value of cubicity metric
CubicDeviationMetric <- function(a,b,c,al,be,ga,decimals=4)
  {
    round((abs((a/sqrt(a**2+b**2-2*a*b*cos(ifelse(ga<90,180-ga,ga)*(pi/180))))-1/sqrt(2))+abs((b/sqrt(a**2+b**2-2*a*b*cos(ifelse(ga<90,180-ga,ga)*(pi/180))))-1/sqrt(2))+abs((c/sqrt(a**2+c**2-2*a*c*cos(ifelse(be<90,180-be,be)*(pi/180))))-1/sqrt(2))+abs((a/sqrt(a**2+c**2-2*a*c*cos(ifelse(be<90,180-be,be)*(pi/180))))-1/sqrt(2))+abs((b/sqrt(b**2+c**2-2*b*c*cos(ifelse(al<90,180-al,al)*(pi/180))))-1/sqrt(2))+abs((c/sqrt(b**2+c**2-2*b*c*cos(ifelse(al<90,180-al,al)*(pi/180))))-1/sqrt(2)))/3,decimals)
  }

#'Unit cell volume calculator
#'
#'This function calculates the volume of the unit cell
#'
#' @param a a lattice parameter, typically in angstroms
#' @param b b lattice parameter, typically in angstroms
#' @param c c lattice parameter, typically in angstroms
#' @param al alpha lattice parameter (angle between b and c sides) in degrees
#' @param be beta lattice parameter (angle between a and c sides) in degrees
#' @param ga gamma lattice parameter (angle between a and b sides) in degrees
#' @param sigfigs Number of sig figs to round answer to. Default is 4.
#' @examples
#' library(tidyverse)
#' UnitCellVolume(10,10,10,90,90,90,decimals=3)
#' UnitCellVolume(10,8,9,63,70,63,decimals=3)
#' df <- data.frame(a=c(10,10),b=c(10,8),c=c(10,9),al=c(90,63),be=c(90,70),ga=c(90,63))
#' df <- df %>% mutate(Volume=UnitCellVolume(a,b,c,al,be,ga,4))
#' head(df)
#' @returns Volume of the unit cell in (units of a/b/c)^3
UnitCellVolume <- function(a,b,c,al,be,ga,sigfigs=4) signif(a*b*c*sqrt(1+2*(cos(al*(pi/180))*cos(be*(pi/180))*cos(ga*(pi/180)))-(cos(al*(pi/180))^2-cos(be*(pi/180))^2-cos(ga*(pi/180))^2)),sigfigs)
