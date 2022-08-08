% Workspace Preparation
% This script is meant to be run at the beginning of each script in this
% project to prepare MATLAB with paths and other code that is redundant in
% each script
% Matt Kmiecik
% Started 15 July 2022

% Sets working directory ----
this_pwd = 'C:\Analysis\crampp-resting-state\';
cd(this_pwd); % goes there

% Initializations ----
raw_data_path = fullfile(this_pwd, 'data\'); % path to raw data
chan_loc_path = fullfile(this_pwd, 'data\chan-32-TP9.bvef'); % path to chan locs
outpath = fullfile(pwd, 'output\'); % path to output data
ss_info = fullfile(raw_data_path, 'ss-info.xlsx');

% Starts EEGLAB ----
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Loads in participant information ----
[NUM,TXT,RAW] = xlsread(ss_info);