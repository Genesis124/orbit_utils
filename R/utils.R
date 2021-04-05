# Author ---------------------------------
# Jordan Brown
# Date: 4/5/2021

# Packaging ------------------------------

library()

# Constants ------------------------------

# Gravitational Constant 
G <- 6.67430 * 10^-11         # m^3 kg^-1 s^-2
# Mass of the Earth
m_E <- 5.972 * 10^24          # kg
# Radius of the Earth
r_E <- 6371                  # km
# Standard Gravitational Parameter 
mu_E <- 3.986004356 * 10^14   # m^3 s^-3

# sidereal rotation period of the Earth relative to the stars
p_srE <- 86164.10035 #s (1 day)

# Spherical Harmonics --------------------
J0 <- 1
J1 <- 0
J2 <- 0.0010826359
J3 <- -0.00000254
J4 <- -0.00000161

# Functions ------------------------------

# e = eccentricity
# M = mean anomaly
# E_n = eccentric anomaly guess
# thres = threshold
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

perigee <- function(a, e){
  return(a*(1-e))
}

# a = semi-major axis
# e = eccentricity
# t = time
# t_e = time of epoch
# t_pg = time of perigee
true_anomaly <- function(a, e, t_pg=0, t_e=0){
  # period
  T <- period(a) # s
  # average rate
  n <- 2*pi/T # rad/s
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

# n_p = number of orbit periods
# n_d = number of days
# e = eccentricity
# i = inclination (rad)
ground_track_repeat <- function(n_p, n_d, e, i){
  k1 <- mu_E^(1/3)*(2*pi/p_srE)^(-2/3)
  k2 <- 0.75*(360/2/pi)*J2*mu_E^(1/2)*r_E^2
  # initial altitude estimate
  H  <- k1*(n_p/n_d)^(-2/3) - r_E
  H0 <- 0
  # initial semi-major axis estimate
  a  <- H0 + r_E
  thres <- 0.000001
  while(abs(H0 - H) > thres){
    # rotation rate of the Earth
    L_dot <- 2*pi/p_srE #rad/s
    # rate of change of the ascending node
    an_dot<- -2*k2*a^(-7/2)*cos(i)*(1-e*e)^-2 #rad/s
    # rate of change of the perigee
    peri_dot <- k2*a^(-7/2)*(5*cos(i)*cos(i)-1)*(1-e*e)^-2 #rad/s
    # rate of change of the mean anomaly
    M_dot    <- k2*a^(-7/2)*(3*cos(i)*cos(i)-1)*(1-e*e)^(-3/2) #rad/s
    # mean angular motion
    n <- (n_p/n_d)*(L_dot-an_dot)-(peri_dot + M_dot) #rad/s
    # altitude estimate
    H0 <- H
    H  <- mu_E^(1/3)*(n*pi/180/p_srE)^(-2/3) - r_E
    # semi-major axis estimate
    a <- H + r_E
  }
  return(H)
}

# Testing --------------------------------

if(T){
  test <- "test"
}
