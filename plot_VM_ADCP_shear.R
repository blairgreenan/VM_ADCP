# plot_VM_ADCP_shear.R
# Blair Greenan
# Fisheries and Oceans Canada
# 7 June 2023
#
# Description: this script generates a plot of the RV NB Palmer vessel-mounted ADCP
#
# load libraries
library(oce)
library(R.matlab)
library(R.oo)
library(tidyverse)
library(lubridate)
library(cmocean)
library(interp)
library(patchwork)
library(ggquiver)
library(plot3D)


# load to processed data file
# https://www.nodc.noaa.gov/archive/arc0039/0082184/1.1/data/0-data/KN201L01/adcp/programs/adcp_doc/UHDAS_scidoc/Processing/adcpsect_output.html
ADCP_uv <- readMat("contour_uv.mat")
ADCP_xy <- readMat("contour_xy.mat")

# Extract the u & v components from the ADCP_uv matrix the alternates u & v 
# throughout the matrix (1:2:end is u, 2:2:end is v)
ADCP_u <- ADCP_uv$uv[,seq(from=1,to=24820,by=2)]
ADCP_v <- ADCP_uv$uv[,seq(from=2,to=24820,by=2)]

# Find the data that corresponds to the days of the Ross Bank VPR Survey
#Start_time <- 383.7 # 18 Jan 2012
#End_time <- 386.0 # 21 Jan 2012
Start_time <- 383.75 # 18 Jan 2012
End_time <- 384.8 # 21 Jan 2012

# Find the indices in the time vector that match this criteria
VPR_survey_time <- which(ADCP_xy$xyt[3,]>=Start_time & ADCP_xy$xyt[3,]<=End_time)

# Extract the latitude, longitude, u and v that correspond to the survey time
# Limited this to the upper 200m, so bins from approx 20-200m - this is the top 36 bins (5m bins)
# u_survey <- ADCP_u[1:36,VPR_survey_time]
# v_survey <- ADCP_v[1:36,VPR_survey_time]
u_survey <- ADCP_u[,VPR_survey_time]
v_survey <- ADCP_v[,VPR_survey_time]
lon_survey <- ADCP_xy$xyt[1,VPR_survey_time]
lat_survey <- ADCP_xy$xyt[2,VPR_survey_time]
time_survey <- ADCP_xy$xyt[3,VPR_survey_time]
# The time in the VM_ADCP data file is set by a base year (2011) and then fractional
# days are added to this.  

# Compute the depth-average ADCP velocity components
u_mean <- colMeans(u_survey,na.rm=TRUE)
v_mean <- colMeans(v_survey,na.rm=TRUE)
# Compute the speed and direction
mean_speed <- sqrt(u_mean^2 + v_mean^2)
mean_direction <- atan2(v_mean,u_mean)

# bin depths for shear - note the the u & v are at zc depths, so the I have assigned
# the shear to be in the middle of those bin depths, whic corresponds to the z values.
z <- ADCP_xy$z[2:79]
u_diff <- diff(u_survey)
v_diff <- diff(v_survey)
shear <- sqrt(u_diff^2 + v_diff^2)/5  # dividing by the bin depth of 5 m

########## Create a tidy data set
# First step is to create a data frame from the 2D matrix
shear_df <- data.frame(shear)
# Convert the data frame to a tibble
shear_tibble <- as_tibble(shear_df)
# Use bind_cols to add a depth vector to the first column of the tibble
shear_tibble_with_depth <- bind_cols(z,shear_tibble)
# Add a name to the depth column
names(shear_tibble_with_depth)[1] <- "Depth"
# User pivot_longer to create a tibble with a shear value for each depth and time 
pivot_shear <- pivot_longer(shear_tibble_with_depth, cols = 2:303, names_to = "Time", values_to = "Shear")
# Convert time from a string to number 
pivot_shear$Time <- as.numeric(pivot_shear$Time)

# Need to figure out how to get Lat and Lon into this tibble.  No recycling with a tibble.

