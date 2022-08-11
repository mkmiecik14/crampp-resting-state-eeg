# CRAMPP Resting State EEG

This repository contains all the code for processing, visualizing, and analysis of the CRAMPP resting state EEG data.

## Processing Pipeline

### All-Purpose Scripts

* `workspace_prep_rs.m` - This MATLAB script prepares the workspace by preparing directory information and loading EEGLAB. You can find it run at the top of most MATLAB processing scripts in this project.

* `r-prep.R` - This R script prepares the R workspace in R scripts. It loads libraries, prepares plotting functions, etc. You can find it run at the top of most R processing scripts in this project.

### Step 1 - Preprocessing EEG data

1) `prepro_resting_state.m`

This MATLAB script preforms some initial preprocessing of EEG data. The order of steps are:
* re-referencing (averaged mastoid reference)
* downsampling (to 256Hz)
* mean centering each channel (DC offset)
* Highpass filter at 1Hz (-6dB @ 1Hz, 425 point highpass, 2Hz transition band width)
* line noise removal using Cleanline plugin

2) Preprocessed data were visually inspected for bad channels and noisy sections of EEG

If noisy sections of EEG were marked for rejection, then these data were manually saved out with the suffix "-prepro-vis-rej.set"

3) `ica_resting_state.m`

This MATLAB script imports visually inspected data, interpolates any bad channels, and performs independent components analysis (ICA).

4) `epochs_spect_resting_state_noise.m`

This MATLAB script performs the following steps
* Labels, flags ICs, and removes ICs (if blink probability > 80%)
* surface Laplacian filtering (i.e., current source density)
* computes broadband power spectral density of each block
* computes peak alpha frequency via peak picking and center of gravity estimates
* estimates and corrects for pink and white noise
* saves out estimates for further analysis in R

### Step 2 - Analyzing EEG Data

1) `prepro-spec-res.R`

This R script reads in:
* all spectral results saved from MATLAB: "../output/*spec-res.mat"
* channel location info: read_xlsx("../data/ss-info.xlsx", sheet = "elec")

and processes these data to prepare for analysis in R. Frequency bandwidths are defined.



