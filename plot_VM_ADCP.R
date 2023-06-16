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
u_survey <- ADCP_u[1:36,VPR_survey_time]
v_survey <- ADCP_v[1:36,VPR_survey_time]
lon_survey <- ADCP_xy$xyt[1,VPR_survey_time]
lat_survey <- ADCP_xy$xyt[2,VPR_survey_time]
time_survey <- ADCP_xy$xyt[3,VPR_survey_time]
# The time in the VM_ADCP data file is set by a base year (2011) and then fractional
# days are added to this.  Need to adjust this some that it is a date/time value
# that can be handled by R and Matlab (to run the Padman tidal model)
DateTime_survey <- date_decimal(2012 + ((time_survey-366)/366))  # scratched my head quite a bit to figure out that this needs to be 366 not 365
DateTime_survey = as.POSIXlt(DateTime_survey, tz = "UTC") #convert to UTC time zone
DateTime_survey = as.numeric(DateTime_survey) / 86400 #convert to days
DateTime_survey = DateTime_survey + 719529 #convert to MATLAB datenum

# write output to a .mat file
writeMat("Survey_metadat.mat", lon_survey=lon_survey, lat_survey=lat_survey, DateTime_survey=DateTime_survey)

###### Start and instance of Matlab and run the TMD (Tidal Model Driver v2.5, Padman et al, 2002)

# Use the R.matlab package functionality to execute the tidal model in Matlab
# First copy the metadata file to the TMD directory
file.copy("Survey_metadat.mat", "C:/Users/greenanb/Documents/Science Projects/Current/Ross Sea/Data/TMD_Matlab_Toolbox_v2_5/TMD/Survey_metadat.mat", overwrite = TRUE)
# Start the instance of Matlab
Matlab$startServer()
# Create a MATLAB client object used to communicate with MATLAB
matlab <- Matlab()
# Check status of MATLAB connection (not yet connected)
print(matlab)
# If you experience any problems, ask for detailed outputs
#     by uncommenting the next line
setVerbose(matlab, -2)
# Connect to the MATLAB server.
isOpen <- open(matlab)
# Confirm that the MATLAB server is open, and running
if (!isOpen)
  throw("MATLAB server is not running: waited 30 seconds.")
# Check status of MATLAB connection (now connected)
print(matlab)
# Change the working directory of Matlab
evaluate(matlab, "cd ('C:/Users/greenanb/Documents/Science Projects/Current/Ross Sea/Data/TMD_Matlab_Toolbox_v2_5/TMD')")
# Load the metadata file in Matlab
evaluate(matlab, "load Survey_metadat.mat")
# Execuate the tidal model functions
evaluate(matlab, "[u_tide,ConList]=tmd_tide_pred('./DATA/CATS2008/Model_CATS2008',DateTime_survey,lat_survey,lon_survey,'u');")
evaluate(matlab, "[v_tide,ConList]=tmd_tide_pred('./DATA/CATS2008/Model_CATS2008',DateTime_survey,lat_survey,lon_survey,'v');")
# Save the results to a .mat file
evaluate(matlab, "save survey_tidal_model.mat DateTime_survey lat_survey lon_survey u_tide v_tide")
# Close the MATLAB client, which will also shutdown
# the MATLAB server and the connection to it.
close(matlab)
# Check status of MATLAB connection (now disconnected)
print(matlab)

########## Return to R and complete the analysis

# Copy tidal model estimates back to the VM_ADCP directory
file.copy("C:/Users/greenanb/Documents/Science Projects/Current/Ross Sea/Data/TMD_Matlab_Toolbox_v2_5/TMD/survey_tidal_model.mat", "C:/Users/greenanb/Documents/Science Projects/Current/Ross Sea/Documents/Papers/Ross Bank/Figures/VM_ADCP/survey_tidal_model.mat", overwrite = TRUE)

# Load the tidal model results
tidal_model <- readMat("survey_tidal_model.mat")
tidal_u <- (tidal_model$u.tide[,1:length(tidal_model$u.tide)])/100
tidal_v <- (tidal_model$v.tide[,1:length(tidal_model$v.tide)])/100

# Compute the depth-average ADCP velocity components
u_mean <- colMeans(u_survey,na.rm=TRUE)
v_mean <- colMeans(v_survey,na.rm=TRUE)
# Compute the speed and direction
mean_speed <- sqrt(u_mean^2 + v_mean^2)
mean_direction <- atan2(v_mean,u_mean)

# Compute the speed and direction from tidal model
tide_speed <- sqrt(tidal_u^2 + tidal_v^2)
tide_direction <- atan2(tidal_v,tidal_u)

# Sub-tidal currents
subtidal_u <- u_mean - tidal_u
subtidal_v <- v_mean - tidal_v
subtidal_speed <- sqrt(subtidal_u^2 + subtidal_v^2)
subtidal_direction <- atan2(subtidal_v,subtidal_u)

# create data frames
ADCP_df <- data.frame(lon_survey,lat_survey,u_mean,v_mean, mean_speed, mean_direction)
tide_df <- data.frame(lon_survey,lat_survey,tidal_u,tidal_v, tide_speed, tide_direction)
subtidal_df <- data.frame(lon_survey,lat_survey,subtidal_u,subtidal_v, subtidal_speed, subtidal_direction)

# plot the data
# Note that the length of the arrows on the plot can be controlled by dividing the mean_speed below
ADCP_quiver_plot <- ggplot(ADCP_df, aes(x=lon_survey,y=lat_survey)) +
  geom_point() + 
  geom_spoke(angle=mean_direction,radius=mean_speed/3, show.legend = TRUE, na.rm = TRUE, arrow = arrow(length = unit(.05, 'inches'))) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("VM_ADCP Depth Average 20-200m")

# Note that the length of the arrows on the plot can be controlled by dividing the mean_speed below
tide_quiver_plot <- ggplot(tide_df, aes(x=lon_survey,y=lat_survey)) +
  geom_point() + 
  geom_spoke(angle=tide_direction,radius=tide_speed/3, show.legend = TRUE, na.rm = TRUE, arrow = arrow(length = unit(.05, 'inches'))) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Tidal Model")

# Note that the length of the arrows on the plot can be controlled by dividing the mean_speed below
subtidal_quiver_plot <- ggplot(subtidal_df, aes(x=lon_survey,y=lat_survey)) +
  geom_point() + 
  geom_spoke(angle=subtidal_direction,radius=subtidal_speed/3, show.legend = TRUE, na.rm = TRUE, arrow = arrow(length = unit(.05, 'inches'))) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Subtidal Currents")


# Plot the results
# dev.new()
# Use Patchwork package to create 3-panel plot
# ADCP_quiver_plot + tide_quiver_plot + subtidal_quiver_plot
# ggsave("VMADCP_tides.png", device = "png", width = 7, height = 3.5, units = "in", dpi = 1200, scale = 2)


# Need to look at the vertical shear estimates from the VM_ADCP



