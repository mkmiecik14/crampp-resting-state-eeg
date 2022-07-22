# R preparation script
# Matt Kmiecik
# Started 22 JULY 2022

# Purpose: to prepare the R workspace
# Any packages, functions, variables, etc. that are used across scripts are
# coded in here to reduce redundancy of code

# Libraries
library(tidyverse)  # general use
library(Hmisc)      # general use
library(R.matlab)   # general use (preprocessing scripts)
library(readxl)     # general use

library(akima)        # topo_tools.R
library(scales)       # topo_tools.R
library(mgcv)         # topo_tools.R
library(RColorBrewer) # topo_tools.R

library(lme4)         # stats
library(broomExtra)   # stats
library(lmerTest)     # stats
library(performance)  # stats
library(effectsize)   # stats
library(MuMIn)        # stats
library(r2glmm)       # stats

# For the API configuration
library(httr)
library(jsonlite)

# Plotting Tools
# use geom_flat_violin()
# see: https://neuroconscience.wordpress.com/2018/03/15/introducing-raincloud-plots/
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
library(patchwork)