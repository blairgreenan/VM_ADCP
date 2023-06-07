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
ADCP <- readMat("contour_uv.mat")

