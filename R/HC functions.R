#'Cv Calculator
#'
#'This function calculates constant volume heat capacity from phonon density of states as a function of frequency (in Hz) for a single temperature. Please see http://lampz.tugraz.at/~hadley/ss1/dbr/dos2cv.html for functional form.
#'
#' @param f list of frequencies; typically input as column in dataframe. Should be in units of Hz or s^-1. *Please note* that this list is assumed to have a constant spacing between points.
#' @param D list of intensities; typically input as column in dataframe. Should be in units of s.
#' @param T single temperature in K
#' @returns Constant volume heat capacity (Cv) in units of J/K
PhononDOStoCv <- function(f, D, T){
  #abbreviation for hbar/Kb*T
  s <- (1.054571817*10^-34)/(1.380649*10^-23*T)
  #sum all values of integrand
  C <- sum((D*(s*f)^2*exp(s*f))/((-1+exp(s*f))^2))
  #subtract off 1/2 first and last point
  C <- C - 0.5*((D[1]*(s*f[1])^2*exp(s*f[1]))/((-1+exp(s*f[1]))^2))-0.5*((D[length(D)]*(s*f[length(f)])^2*exp(s*f[length(f)]))/((-1+exp(s*f[length(f)]))^2))
  #multiply by stepsize, which is assumed to be constant so computed from first interval
  C <- abs(f[[1]]-f[[2]])*C
  #multiply by k
  C <- C*(1.380649*10^-23)
  return(C)
}

#'Cv vs. T calculator
#'
#'This function calculates constant volume heat capacity from phonon density of states as a function of frequency (in Hz) for a range of temperatures. It also calculates Cv/T^3. Please see http://lampz.tugraz.at/~hadley/ss1/dbr/dos2cv.html for functional form.
#'
#' @param f list of frequencies; typically input as column in dataframe. Should be in units of Hz or s^-1. *Please note* that this list is assumed to have a constant spacing between points.
#' @param D list of intensities; typically input as column in dataframe. Should be in units of s.
#' @param Tmin lowest temperature to calculate in Kelvin
#' @param Tmax highest temperature to calculate in Kelvin
#' @param Tstep stepsize between each temperature to calculate. For example, with Tmin=1, Tmax=5, Tstep=1... the function will return Cv's for T=1,2,3,4,5
#' @returns Dataframe containing column for Temperature in units of K, column for Constant volume heat capacity (Cv) in units of J/K, and column for Constant volume heat capacity divided by temperature cubed (CvperT3) in units of J/K^4
Cv.Tlist <- function(f, D, Tmin=2, Tmax=300, Tstep=1){
  Tlist <- seq(Tmin,Tmax,Tstep)
  Cv <- vector()
  for (T in Tlist) {
    Cv <- c(Cv,PhononDOStoCv(f,D,T))
  }
  Temperature <- c(Tlist)
  CvperT3 <- Cv/(Temperature^3)
  return(data.frame(Temperature,Cv,CvperT3))
}

#'Heat Capacity Plotter - single dataset
#'
#'This function makes a pretty picture of heat capacity data with labelled axes.
#'
#' @param dataframe Dataframe containing the data you want to plot.
#' @param tempvals Name of the column containing temperature data in Kelvin.
#' @param yvals Name of the column containing Cp or Cv data in J/K.
#' @param title String to be main title of your plot.
#' @param isCv Boolean - when TRUE, units come out for constant volume heat capacity. When FALSE, units come out for constant pressure.
#' @param perT3 Boolean - when TRUE, units are given for C/T^3. When FALSE, units are given for C.
#' @param ptcolor Color for the plot
HCPlot <- function(dataframe,tempvals=Temperature,yvals=Cv,title="Heat capacity",isCv=TRUE,perT3=TRUE,ptcolor="darkgreen"){
  ggplot(dataframe, aes(.data[[tempvals]],.data[[yvals]])) + geom_point(shape=1,size=3,color=ptcolor) +
    scale_x_continuous(
      trans = "log10",
      breaks = scales::log_breaks(n = 5)
    ) +
    labs(title = title,
         y=if(isCv==TRUE){if(perT3 == TRUE){
           expression(paste("Cv/", T^3," (J/", K^4,")" ) )}else{
             expression(paste("Cv"," (J/", K,")" ) )
           }}else{if(perT3 == TRUE){
         expression(paste("Cp/", T^3," (J/", K^4,")" ) )}else{
           expression(paste("Cp"," (J/", K,")" ) )
         }},
         x="Temperature (K)") +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5, size=18),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text = element_text(size = 14),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "gray")) +
    annotation_logticks(sides="b")
}

#'Heat Capacity Plotter - overlaid Cp and Cv
#'
#'This function makes a pretty picture of heat capacity data with labelled axes.
#'
#' @param CvTemps Column in a dataframe containing temperature data for Cv plot in Kelvin, as in dataframe$column
#' @param CvVals Column in a dataframe containing Cv data in J/K, as in dataframe$column
#' @param CvCol Color for Cv data and axis
#' @param CpTemps Column in a dataframe containing temperature data for Cp plot in Kelvin, as in dataframe$column
#' @param CpVals Column in a dataframe containing Cp data in J/K, as in dataframe$column
#' @param CpCol Color for Cp data and axis
#' @param CpCoeff Scaling coefficient for the Cp data - guess and check until your data shows up with appropriate magnitude
#' @param maintitle String to be main title of your plot.
#' @param perT3 Boolean - when TRUE, units are given for C/T^3. When FALSE, units are given for C.
CvCpOverlayPlot <- function(CvTemps,CvVals,CvCol="#1D3732",CpTemps,CpVals,CpCol="#397F72",CpCoeff=1,maintitle="Heat capacity overlay",perT3=TRUE){
  ggplot() + geom_line(aes(CvTemps, CvVals),color=CvCol,size=2)+
    geom_point(aes(CpTemps, CpVals*CpCoeff),shape=1,size=3,color=CpCol) +
    scale_x_continuous(
      trans = "log10",
      breaks = scales::log_breaks(n = 5)
    ) +
    scale_y_continuous(

      # Features of the first axis
      name = if(perT3==TRUE){expression(paste("Cv/", T^3," from DOS calculation (J/ ", K^4,")" ) )}else{expression(paste("Cv/", T^3," from DOS calculation (J/ ", K,")" ) )},

      # Add a second axis and specify its features
      sec.axis = sec_axis( transform=~.*1/CpCoeff,
                           name= if(perT3==TRUE){expression(paste(" Cp/", T^3,"from experiment (J/ ", K^4,")" ) )}else{expression(paste(" Cp/", T^3,"from experiment (J/ ", K,")" ) )})
    ) +
    labs(title = maintitle,
         x="Temperature (K)") +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5, size=18),
          axis.title.y = element_text(size = 16),
          axis.title.y.right = element_text(color = CpCol,size = 16),
          axis.line.y.right = element_line(color= CpCol),
          axis.text.y.right= element_text(color = CpCol),
          axis.line.y.left = element_line(color=CvCol),
          axis.title.y.left = element_text(color = CvCol,size = 16),
          axis.text.y.left= element_text(color = CvCol),
          axis.title.x = element_text(color="black",size = 16),
          axis.line.x.bottom = element_line(color="black"),
          axis.text = element_text(size = 14,color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    annotation_logticks(sides="b")
}
