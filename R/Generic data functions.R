#'Numeric integrator (constant return)
#'
#'This function performs a numeric integration. In other words, it calculates area under curve defined by series of (x,y) points. The method used is a trapezoidal Riemann sum.
#'
#' @param x list of independent values
#' @param y list of dependent values
#' @param const.step Boolean descriptor for whether (x,y) pairs are evenly spaced. If you do not know, set this to FALSE. In the case the pairs are evenly spaced, the function will use a slightly faster integration method, but the standard method will work in both cases.
#' @examples
#' x<-seq(1,10,1)
#' y<-x^2
#' curve.area(x,y,const.step=TRUE)
#' curve.area(x,y,const.step=FALSE)
#' @returns Integral (aka area under the curve) of input data.
curve.area <- function(x, y, const.step=FALSE){
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  if(const.step==TRUE){
  #sum all values of integrand then subtract off 1/2 first and last point
  result <- sum(y)
  result <- result - 0.5*(y[1])-0.5*(y[length(y)])
  #multiply by stepsize, which is assumed to be constant so computed from first interval
  area <- abs(x[[1]]-x[[2]])*result}
  else{
    n <- length(y)
    results <- vector(length = n) #initialize vector to hold results
    #multiply sum of ith and previous y value by the difference between the ith and previous x value
    for (i in 2:n) {
      results[i] <- 0.5*(y[i] + y[i-1]) * (x[i] - x[i-1])
      area <- sum(results)
  }
}
  return(area)
}

#'Numeric integrator (dataset return)
#'
#'This function performs a numeric integration (quadrature) of a series of (x,y) points. The method used is a trapezoidal Riemann sum.
#'
#' @param x list of independent values
#' @param y list of dependent values
#' @returns Integral of input data as a vector of new integral values. First point in list will be NA.
Riemann.sum <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  n <- length(x)
  results <- vector(length = n)
  results[1] <- NA
  for (i in 2:n) {
    results[i] <- 0.5*(y[i] + y[i-1]) * (x[i] - x[i-1])
  }
  return(results)
}

#'Numeric differentiator (dataset return)
#'
#'This function performs a numeric differentiation. In other words, it calculates the numeric derivative of the input data.
#'
#' @param x list of independent values
#' @param y list of dependent values
#' @returns Integral of input data as a vector of new derivative values. First point in list will be NA.
numeric.derivative <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  # Initialize a vector of length n to enter the derivative approximations
  n <- length(x)
  results <- vector(length = n)
  # Iterate through the values
  #For each datapoint up to the nth one, divide the difference between y1 and y0 by
  #the difference between x1 and x0 (producing dy/dx)
  results[1] <- NA
  for (i in 2:n) {
    results[i] <- (y[i] - y[i-1]) / (x[i] - x[i-1])
  }
  return(results)
}
