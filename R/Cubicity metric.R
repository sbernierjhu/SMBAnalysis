
#'Cubic Deviation Metric calculator
#'
#'This function calculates the value of the cubicity metric described in xxx.
#'
#'The returned value of the cubicity metric will be zero if and only if the lattice parameters input correspond to a cubic unit cell. The returned value will be unitless and normalized such that it can be compared to values for other unit cells. The closer to zero, the "more cubic" your unit cell.
#'
#' @param a a lattice parameter, typically in angstroms
#' @param b b lattice parameter, typically in angstroms
#' @param c c lattice parameter, typically in angstroms
#' @param al alpha lattice parameter (angle between b and c sides) in degrees
#' @param be beta lattice parameter (angle between a and c sides) in degrees
#' @param ga gamma lattice parameter (angle between a and b sides) in degrees
#' @param sigfigs Number of sig figs to round answer to. Default is 4.
#' @returns Value of cubicity metric
CubicityMetric <- function(a,b,c,al,be,ga,sigfigs=4) signif(abs((a/sqrt(a**2+b**2-2*a*b*cos(ga*(pi/180))))-1/sqrt(2))+abs((b/sqrt(a**2+b**2-2*a*b*cos(ga*(pi/180))))-1/sqrt(2))+abs((c/sqrt(a**2+c**2-2*a*c*cos(be*(pi/180))))-1/sqrt(2))+abs((a/sqrt(a**2+c**2-2*a*c*cos(be*(pi/180))))-1/sqrt(2))+abs((b/sqrt(b**2+c**2-2*b*c*cos(al*(pi/180))))-1/sqrt(2))+abs((c/sqrt(b**2+c**2-2*b*c*cos(al)*(pi/180)))-1/sqrt(2)),sigfigs)

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
#' @returns Volume of the unit cell in (units of a/b/c)^3
UnitCellVolume <- function(a,b,c,al,be,ga,sigfigs=4) signif(a*b*c*sqrt(1+2*(cos(al*(pi/180))*cos(be*(pi/180))*cos(ga*(pi/180)))-(cos(al*(pi/180))^2-cos(be*(pi/180))^2-cos(ga*(pi/180))^2)),sigfigs)
