# Author ---------------------------------
# Jordan Brown
# Date: 4/5/2021

# Packaging ------------------------------

library("pracma")

# Constants ------------------------------

# Gravitational Constant 
G <- 6.67430 * 10^-20         # km^3 kg^-1 s^-2
# Mass of the Earth
m_E <- 5.972 * 10^24          # kg
# Radius of the Earth
r_E <- 6378.136               # km
# Standard Gravitational Parameter 
mu_E <- 3.986004356 * 10^5    # km^3 s^-3

# sidereal rotation period of the Earth relative to the stars
p_srE <- 86164.10035 #s (1 day)

# Spherical Harmonics --------------------
J0 <- 1
J1 <- 0
J2 <- 0.0010826359
J3 <- -0.00000254
J4 <- -0.00000161

# Functions ------------------------------

# Eccentric anomaly equation has no closed form solution
# M = E - e*sin(E)
# This function iteratively solves for E
# e = eccentricity
# M = mean anomaly
# E_n = eccentric anomaly guess
# thres = threshold
# Returns: Eccentric anomaly 
NewtonRaphson <- function(e, M, E_n = M, thres = 10^-6){
  E_m <- 2*pi
  while(abs(E_m - E_n) < thres){
    E_m <- E_n - (1 - e*sin(E_n) - M)/(1 - e*cos(E_n))
  }
  return(E_m)
}

# Orbital Functions ----------------------

# a = semi-major axis
period <- function(a){
  return(2*pi*sqrt(a*a*a/mu_E))
}

# a = semi-major axis
# e = eccentricity
apogee <- function(a, e){
  return(a*(1+e))
}

# a = semi-major axis
# e = eccentricity
perigee <- function(a, e){
  return(a*(1-e))
}

# Determines the true anomaly 
# a = semi-major axis
# e = eccentricity
# t = time
# t_e = time of epoch
# t_pg = time of perigee
# Returns: True Anomaly (rad)
true_anomaly <- function(a, e, t, t_pg=0, t_e=0){
  # period
  p <- period(a) # s
  # average rate
  n <- 2*pi/p # rad/s
  # mean anomaly at epoch
  M0 <- n*(t_e-t_pg) # rad
  # mean anomaly 
  M <- M0 + n*(t - t_e) # rad
  # eccentric anomaly
  E <- NewtonRaphson(e, M)
  # true anomaly
  v <- acos((cos(E_t) - e)/(1-e*cos(E_t)))
  return(v)
}

# Determines the altitude required for a repeating groundtrack
# n_p = number of orbit periods
# n_d = number of days
# e = eccentricity
# i = inclination (rad)
# Returns: altitude (km)
ground_track_repeat <- function(n_p, n_d, e, i){
  k1 <- mu_E^(1/3)*(2*pi/p_srE)^(-2/3) # km
  k2 <- 0.75*(360/2/pi)*J2*mu_E^(1/2)*r_E^2 #deg/s
  # initial altitude estimate
  H  <- k1*(n_p/n_d)^(-2/3) - r_E
  H0 <- 0
  # initial semi-major axis estimate
  a  <- H + r_E
  # continues looping until difference is less than threshold
  thres <- 0.000001 
  while(abs(H0 - H) > thres){
    # rotation rate of the Earth
    L_dot <- 360/p_srE #deg/s
    # rate of change of the ascending node
    an_dot<- -2*k2*a^(-7/2)*cos(i)*(1-e*e)^-2 #deg/s
    # rate of change of the perigee
    peri_dot <- k2*a^(-7/2)*(5*cos(i)*cos(i)-1)*(1-e*e)^-2 #deg/s
    # rate of change of the mean anomaly
    M_dot    <- k2*a^(-7/2)*(3*cos(i)*cos(i)-1)*(1-e*e)^(-3/2) #deg/s
    # mean angular motion
    n <- (n_p/n_d)*(L_dot-an_dot)-(peri_dot + M_dot) #deg/s
    # altitude estimate
    H0 <- H
    H  <- mu_E^(1/3)*(n*pi/180)^(-2/3) - r_E
    # semi-major axis estimate
    a <- H + r_E
  }
  return(H)
}

# Earth Imaging --------------------------

# H = altitude
# Returns: max earth central angle (rad)
max_earth_central_ang <- function(H){return(acos(r_E/(r_E + H)))}
# H = altitude
# nadir = nadir angle; angle between line of sight and nadir point
# tea = target elevation angle; angle between line of sight and ground
# Returns: earth central angle (rad)
earth_central_ang <- function(tea, nadir){return(pi/2 - tea - nadir)}
# H = altitude
# Returns: angular radius of the Earth (rad)
ang_rad <- function(H){asin(r_E/(r_E + H))}
# H = altitude
# lambda = central earth angle; earth angle btw nadir and vector from earth center to target
# tea = target elevation angle; angle between line of sight and ground
# Returns: angular radius of the Earth (rad)
nadir_ang <- function(H = NULL, lambda = NULL, tea = NULL){
  if(is.null(H)){return(pi/2-lambda-tea)}
  a_rad <- ang_rad(H)
  if(!is.null(tea)){return(asin(cos(tea)*sin(a_rad)))}
  if(!is.null(lambda)){return((sin(a_rad)*sin(lambda))/(1-sin(a_rad)*sin(lambda)))}
  stop("insufficient arguments")
}
# H = altitude
# lambda = central earth angle; earth angle btw nadir and vector from earth center to target
# nadir  = nadir angle; angle btw line of sight and nadir
# Returns: target elevation angle (rad)
tgt_elev_ang <- function(H = NULL, lambda = NULL, nadir = NULL){
  if(is.null(H)){return(pi/2-lambda-nadir)}
  a_rad <- ang_rad(H)
  if(!is.null(nadir)){return(acos(sin(nadir)/sin(a_rad)))}
  stop("insufficient arguments")
}
# lambda = central earth angle; earth angle btw nadir and vector from earth center to target
# nadir  = nadir angle; angle between line of sight and nadir
# Returns: slant range; distance btw satellite and ground
slant_range <- function(lamda, nadir){return(r_E*sin(lamda)/sin(nadir))} 

# Testing --------------------------------

if(T){
  test <- "test"
  p <- 14
  d <- 1
  e <- .12
  i <- deg2rad(28)
  data <- ground_track_repeat(p, d, e, i)
  print(data)
}
