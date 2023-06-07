# plot_VM_ADCP.R
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
u_survey <- ADCP_u[1:36,VPR_survey_time]
v_survey <- ADCP_v[1:36,VPR_survey_time]
lon_survey <- ADCP_xy$xyt[1,VPR_survey_time]
lat_survey <- ADCP_xy$xyt[2,VPR_survey_time]

# Compute the depth-average velocity components
u_mean <- colMeans(u_survey,na.rm=TRUE)
v_mean <- colMeans(v_survey,na.rm=TRUE)

# Compute the speed and direction
mean_speed <- sqrt(u_mean^2 + v_mean^2)
mean_direction <- atan2(v_mean,u_mean)

ADCP_df <- data.frame(lon_survey,lat_survey,u_mean,v_mean, mean_speed, mean_direction)

# plot the data
# Note that the length of the arrows on the plot can be controlled by dividing the mean_speed below
quiver_plot <- ggplot(ADCP_df, aes(x=lon_survey,y=lat_survey)) +
  geom_point() + 
  geom_spoke(angle=mean_direction,radius=mean_speed/3, show.legend = TRUE, na.rm = TRUE, arrow = arrow(length = unit(.05, 'inches'))) +
  xlab("Longitude") +
  ylab("Latitude")



