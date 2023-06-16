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
# Start_time <- 383.75 # 18 Jan 2012
# End_time <- 384.8 # 21 Jan 2012
# Find the data that corresponds to the days of the Eddy 3 Survey
# Start_time <- 382 # 17 Jan 2012
# End_time <- 383.5 # 18 Jan 2012
# Find the data that corresponds to the days of the Ice Edge Eddy Survey
 Start_time <- 395.5 # 31 Jan 2012
 End_time <- 396.5 # 01 Feb 2012
 

# Find the indices in the time vector that match this criteria
VPR_survey_time <- which(ADCP_xy$xyt[3,]>=Start_time & ADCP_xy$xyt[3,]<=End_time)

# Extract the latitude, longitude, u and v that correspond to the survey time
# Limited this to the upper 200m, so bins from approx 20-200m - this is the top 36 bins (5m bins)
# u_survey <- ADCP_u[1:36,VPR_survey_time]
# v_survey <- ADCP_v[1:36,VPR_survey_time]
u_survey <- ADCP_u[,VPR_survey_time]
u_survey[is.nan(u_survey)] <- NA
v_survey <- ADCP_v[,VPR_survey_time]
v_survey[is.nan(v_survey)] <- NA
lon_survey <- ADCP_xy$xyt[1,VPR_survey_time]
lon_survey[is.nan(lon_survey)] <- NA
lat_survey <- ADCP_xy$xyt[2,VPR_survey_time]
lat_survey[is.nan(lat_survey)] <- NA
time_survey <- ADCP_xy$xyt[3,VPR_survey_time]
time_survey[is.nan(time_survey)] <- NA
# The time in the VM_ADCP data file is set by a base year (2011) and then fractional
# days are added to this.  

# Compute the depth-average ADCP velocity components
u_mean <- colMeans(u_survey,na.rm=TRUE)
v_mean <- colMeans(v_survey,na.rm=TRUE)
# Compute the speed and direction
mean_speed <- sqrt(u_mean^2 + v_mean^2)
mean_direction <- atan2(v_mean,u_mean)

# bin depths for shear - note the the u & v are at zc depths, so the I have assigned
# the shear to be in the middle of those bin depths, which corresponds to the z values.
z <- ADCP_xy$z[2:79]
u_diff <- diff(u_survey)
v_diff <- diff(v_survey)
shear <- sqrt(u_diff^2 + v_diff^2)/5  # dividing by the bin depth of 5 m

########## Create a tidy data set
# First step is to create a data frame from the 2D matrix
shear_df <- data.frame(shear)
# Convert the data frame to a tibble
shear_tibble <- as_tibble(shear_df)
# Add time as names of the columns
names(shear_tibble) <- as.character(time_survey)
# Use bind_cols to add a depth vector to the first column of the tibble
shear_tibble_with_depth <- bind_cols(z,shear_tibble)
# Add a name to the depth column
names(shear_tibble_with_depth)[1] <- "Depth"
# User pivot_longer to create a tibble with a shear value for each depth and time 
#pivot_shear <- pivot_longer(shear_tibble_with_depth, cols = 2:303, names_to = "Time", values_to = "Shear")
len_time_survey <- length(time_survey)+1
pivot_shear <- pivot_longer(shear_tibble_with_depth, cols = 2:len_time_survey, names_to = "Time", values_to = "Shear")
# Convert time from a string to number 
pivot_shear$Time <- as.numeric(pivot_shear$Time)
# Use the append function to create a latitude vector of the same length as others in pivot_shear since there is no recycling in a tibble
Latitude <- lat_survey
for (i in 2:length(z)) {Latitude <- append(Latitude, lat_survey)}
pivot_shear$Latitude <- Latitude
# Use the append function to create a longitude vector of the same length as others in pivot_shear since there is no recycling in a tibble
Longitude <- lon_survey
for (i in 2:length(z)) {Longitude <- append(Longitude, lon_survey)}
pivot_shear$Longitude <- Longitude

# Plot the results
library(viridisLite)
library(plotly)
fig <- plot_ly(pivot_shear, x = ~Longitude, y = ~Latitude, z = ~Depth*-1, color = ~Shear)
fig <- fig %>% add_markers()
fig
# Need to type fig at command prompt to get this to display


dev.new()
scatter3D(pivot_shear$Longitude, pivot_shear$Latitude, -1*pivot_shear$Depth, colvar = pivot_shear$Shear,
phi = 45, theta = 45, col = viridis(256), pch = 19, cex = 0.75, cex.main = 2,
cex.axis = 0.75, cex.lab = 0.75, xlab = "Longitude", ylab = "Latitude", zlab = "Depth (m)",
main = "Shear", ticktype = "detailed")

# Plot a map of the locations of the VM_ADCP samples 
dev.new()
plot(pivot_shear$Longitude, pivot_shear$Latitude)
# Plot the mean 20-200m currents for context of the shear plots
# create data frames
ADCP_df <- data.frame(lon_survey,lat_survey,u_mean,v_mean, mean_speed, mean_direction)
# plot the data
# Note that the length of the arrows on the plot can be controlled by dividing the mean_speed below
ADCP_quiver_plot <- ggplot(ADCP_df, aes(x=lon_survey,y=lat_survey)) +
  geom_point() + 
  geom_spoke(angle=mean_direction,radius=mean_speed/3, show.legend = TRUE, na.rm = TRUE, arrow = arrow(length = unit(.05, 'inches'))) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("VM_ADCP Depth Average 20-200m")
# Need to type dev.new() and ADCP_quiver_plot at command prompt to get this to display
