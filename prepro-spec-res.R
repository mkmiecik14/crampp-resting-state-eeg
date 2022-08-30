# Preprocessing script to organize the spectral data into R
# Matt Kmiecik
# Started 22 JULY 2022

source("r-prep.R") # preps R workspace 

# Loading various data - - - -

# Loading in spectral results from matlab (.mat)
spec_files <- list.files("../output/", pattern = "*spec-res.mat")

# Loading in channel locations
chan_locs <- 
  read_xlsx("../data/ss-info.xlsx", sheet = "elec") %>%
  mutate(labels = gsub("'", "", labels)) # removes apostrophe

# saving out channel locations
save(chan_locs, file = "../output/chan-locs.rda")
write_csv(chan_locs, file = "../output/chan-locs.csv")

# Gathering subject numbers
subjs <- 
  tibble(fnames = spec_files) %>% 
  separate(fnames, into = c("task", "visit", "ss", "type1", "type2", "ext")) %>%
  mutate(
    ss_i = 1:n(), # to help join below / for specificity 
    ss = as.numeric(ss)
    )

# trigger definitions
triggers <-
  tibble(
    trigger = c(111, 102, 103, 114, 105, 116, 117, 108),
    eyes = c(
      "open", "closed", "closed", "open", "closed", "open", "open", "closed"
    ),
    block = 1:8
  )

# Reading in and unpacking spectral results ----
spec_res <- 
  spec_files %>%
  map(~readMat(file.path("../output/", .x))) %>%
  map("spec.res") %>%
  map(~as.matrix(.x))

# Pulling out various spectral results - - - -
# resting-state blocks 1-8 along z dim in numerical order
stim_spectra  <- spec_res %>% map(1) # broadband PSD (in dB)
stim_freqs <- spec_res %>% map(2) # frequencies
freqs <- stim_freqs[[1]][,,1] # vector of frequencies
stim_paf <- spec_res %>% map(3) # peak alpha frequency
stim_cog <- spec_res %>% map(4) # center of gravity
stim_iaf <- spec_res %>% map(5) # individual alpha freq
# broadband PSD (converted to uV^2/Hz + surface Laplacian)
stim_psd <- spec_res %>% map(6) 
# broadband PSD (converted to uV^2/Hz + surface Laplacian) corrected for 
# pink and white noise
stim_psd_cor <- spec_res %>% map(7)
stim_freqvec <- spec_res %>% map(8) # frequencies for noise correction
pwn_freqs <-  stim_freqvec[[1]][,1] # frequency vector for noise correction

# Cleanup
rm(spec_res) # removes from memory 
gc() # garbage collection

#####################################
#                                   #
# BROADBAND SPECTRAL RESULTS in dB  #
#                                   #
#####################################

# Extracting and tidying spectral results

# frequency resolution is determined in MATLAB spectral decomposition scripts
# so see those scripts for these numbers:
eeg_srate <- 256 # data sampled at 256 Hz
fft_window <- 4 # in seconds
n_points <- eeg_srate*fft_window # number of data points in FFT window
freq_bins <- eeg_srate/n_points # frequency bins

# defining frequency bands (for later)
delta <-  seq(.5, 4, freq_bins)
theta <-  seq(4, 7.5, freq_bins)
alpha <-  seq(7.5, 13, freq_bins)
beta <-   seq(13, 30, freq_bins)

# Each participant has:
# 30 rows (electrodes w/o ref) * n (depends on frequency bins) * 8 (rs blocks)
dB_res <-
  stim_spectra %>%
  map_dfr(
    ~bind_rows(apply(.x, 3, as_tibble), .id = "block") %>% 
      mutate(elec = rep(chan_locs$labels, 8)) %>%
      relocate(elec, .after = block) %>%
      pivot_longer(c(-block, -elec)) %>%
      mutate(name = rep(freqs, 8*length(chan_locs$labels))),
    .id = "ss_i"
  ) %>%
  # converts to numeric
  mutate(
    ss_i = as.numeric(ss_i),
    block = as.numeric(block)
  ) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>% # joins with subject numbers
  left_join(., triggers %>% select(block, eyes), by = "block") %>%
  rename(freq = name, psd = value) %>% # renaming
  select(ss, block, eyes, elec, freq, psd) %>% # reordering + deselecting (ss_i)
  # injects frequency band labels here
  mutate(
    band = case_when(
      freq %in% delta ~ "delta",
      freq %in% theta ~ "theta",
      freq %in% alpha ~ "alpha",
      freq %in% beta ~ "beta",
      TRUE ~ "outside"
    )
  )

# check for missing blocks
# dB_res_blocks <- 
#   dB_res %>%
#   filter(complete.cases(psd)) %>%
#   select(ss:elec) %>% 
#   distinct() %>% 
#   count(ss, block, eyes) %>% 
#   count(ss, eyes)
# 
# dB_res_blocks %>% filter(n<4)
  
# Saving out broadband spectral results
save(dB_res, file = "../output/dB-res.rda")
write_csv(dB_res, file = "../output/dB-res.csv")
rm(dB_res) # removes from memory
gc() # garbage collection

########################
#                      #
# PEAK ALPHA FREQUENCY #
#                      #
########################

# Cleaning up peak alpha frequency
paf_res <- 
  stim_paf %>%
  map_df(~as_tibble(.x) %>% mutate(elec = chan_locs$labels), .id = "ss_i") %>%
  mutate(ss_i = as.numeric(ss_i)) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>%
  select(-ss_i) %>%
  pivot_longer(cols = c(-ss, -elec), names_to = "block", values_to = "paf") %>%
  mutate(block = as.numeric(regmatches(block, regexpr("\\d", block)))) %>%
  left_join(., triggers %>% select(block, eyes), by = "block") %>%
  select(ss, elec, block, eyes, paf)

# check for missing blocks
# paf_res_blocks <- 
#   paf_res %>%
#   filter(complete.cases(paf)) %>%
#   select(ss:eyes) %>% 
#   distinct() %>% 
#   count(ss, block, eyes) %>% 
#   count(ss, eyes)
# 
# paf_res_blocks %>% filter(n<4)

# Saving out paf results
save(paf_res, file = "../output/paf-res.rda")
write_csv(paf_res, file = "../output/paf-res.csv")
rm(paf_res) # removes from memory
gc() # garbage collection

#####################
#                   #
# CENTER OF GRAVITY #
#                   #
#####################

# Cleaning up cog results
cog_res <-
  stim_cog %>%
  map_df(~as_tibble(.x) %>% mutate(elec = chan_locs$labels), .id = "ss_i") %>%
  mutate(ss_i = as.numeric(ss_i)) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>%
  select(-ss_i) %>%
  pivot_longer(cols = c(-ss, -elec), names_to = "block", values_to = "cog") %>%
  mutate(block = as.numeric(regmatches(block, regexpr("\\d", block)))) %>%
  left_join(., triggers %>% select(block, eyes), by = "block") %>%
  select(ss, elec, block, eyes, cog)

# check for missing blocks
# cog_res_blocks <-
#   cog_res %>%
#   filter(complete.cases(cog)) %>%
#   select(ss:eyes) %>%
#   distinct() %>%
#   count(ss, block, eyes) %>%
#   count(ss, eyes)
# 
# cog_res_blocks %>% filter(n<4)

# Saving out cog results
save(cog_res, file = "../output/cog-res.rda")
write_csv(cog_res, file = "../output/cog-res.csv")
rm(cog_res) # removes from memory
gc() # garbage collection

##############################################
#                                            #
# INDIVIDUAL ALPHA FREQUENCY - GRANDAVERAGES #
#                                            #
##############################################

# Cleaning up iaf results
iaf_res <- 
  stim_iaf %>%
  map_df(
    ~as_tibble(.x) %>% 
      rename(paf = V1, cog = V2) %>% 
      mutate(block = 1:n()),
    .id = "ss_i"
    ) %>%
  mutate(ss_i = as.numeric(ss_i)) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>%
  left_join(., triggers %>% select(block, eyes), by = "block") %>%
  select(ss, block, eyes, paf, cog)

# check for missing blocks
# iaf_res_blocks <-
#   iaf_res %>%
#   filter(complete.cases(.)) %>%
#   count(ss, block, eyes) %>%
#   count(ss, eyes)
# 
# iaf_res_blocks %>% filter(n<4)

# Saving out iaf results
save(iaf_res, file = "../output/iaf-res.rda")
write_csv(iaf_res, file = "../output/iaf-res.csv")
rm(iaf_res) # removes from memory
gc() # garbage collection

#######################################
#                                     #
# BROADBAND SPECTRAL RESULTS in PSD   #
#             uV^2/Hz                 #
#                                     #
#######################################

# Technically, the signals were filtered via surface Laplacian prior to FFT. 
# The surface Laplacian puts the units as uV/cm^2 then an FFT is performed, 
# putting the final units as ((uV/cm^2)^2)/Hz.....I think

# Extracting and tidying spectral results
psd_res <-
  stim_psd %>%
  map_dfr(
    ~bind_rows(apply(.x, 3, as_tibble), .id = "block") %>% 
      mutate(elec = rep(chan_locs$labels, 8)) %>%
      relocate(elec, .after = block) %>%
      pivot_longer(c(-block, -elec)) %>%
      mutate(name = rep(freqs, 8*length(chan_locs$labels))),
    .id = "ss_i"
  ) %>%
  # converts to numeric
  mutate(
    ss_i = as.numeric(ss_i),
    block = as.numeric(block)
  ) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>% # joins with subject numbers
  left_join(., triggers %>% select(block, eyes), by = "block") %>%
  rename(freq = name, psd = value) %>% # renaming
  select(ss, block, eyes, elec, freq, psd) %>% # reordering + deselecting (ss_i)
  # injects frequency band labels here
  mutate(
    band = case_when(
      freq %in% delta ~ "delta",
      freq %in% theta ~ "theta",
      freq %in% alpha ~ "alpha",
      freq %in% beta ~ "beta",
      TRUE ~ "outside"
    )
  )

# check for missing blocks
# psd_res_blocks <- 
#   psd_res %>%
#   filter(complete.cases(psd)) %>%
#   select(ss:elec) %>% 
#   distinct() %>% 
#   count(ss, block, eyes) %>% 
#   count(ss, eyes)
# 
# psd_res_blocks %>% filter(n<4)
  
# Saving out broadband spectral results
save(psd_res, file = "../output/psd-res.rda")
write_csv(psd_res, file = "../output/psd-res.csv")
rm(psd_res) # removes from memory
gc() # garbage collection

################################################################
#                                                              #
# SPECTRAL RESULTS in PSD CORRECTED FOR PINK AND WHITE NOISE   #
#                       units:  uV^2/Hz                        #
#                                                              #
################################################################

# Technically, the signals were filtered via surface Laplacian prior to FFT. 
# The surface Laplacian puts the units as uV/cm^2 then an FFT is performed, 
# putting the final units as ((uV/cm^2)^2)/Hz.....I think

psd_cor_res <-
  stim_psd_cor %>%
  map_dfr(
    ~bind_rows(apply(.x, 3, as_tibble), .id = "block") %>%
      mutate(freq = rep(pwn_freqs, 8)) %>%
      relocate(freq, .after = block) %>%
      pivot_longer(c(-block, -freq), names_to = "elec", values_to = "psd_cor") %>%
      mutate(elec = rep(chan_locs$labels, 8*length(pwn_freqs))),
    .id = "ss_i"
    ) %>%
  # converts to numeric
  mutate(
    ss_i = as.numeric(ss_i),
    block = as.numeric(block)
  ) %>% 
  left_join(., subjs %>% select(ss_i, ss), by = "ss_i") %>% # joins with subject numbers
  left_join(., triggers %>% select(block, eyes), by = "block") %>%
  select(ss, block, eyes, elec, freq, psd_cor) %>% # reordering + deselecting (ss_i)
  # injects frequency band labels here
  mutate(
    band = case_when(
      freq %in% delta ~ "delta",
      freq %in% theta ~ "theta",
      freq %in% alpha ~ "alpha",
      freq %in% beta ~ "beta",
      TRUE ~ "outside"
    )
  )

# check for missing blocks
# psd_cor_res_blocks <- 
#   psd_cor_res %>%
#   filter(complete.cases(psd_cor)) %>%
#   select(ss:elec) %>% 
#   distinct() %>% 
#   count(ss, block, eyes) %>% 
#   count(ss, eyes)
# 
# psd_res_blocks %>% filter(n<4)

# Saving out pink&white noise corrected spectral results
save(psd_cor_res, file = "../output/psd-cor-res.rda")
write_csv(psd_cor_res, file = "../output/psd-cor-res.csv")
rm(psd_cor_res) # removes from memory
gc() # garbage collection

# Cleans work space ----
rm(
  chan_locs,
  stim_freqs,
  stim_spectra,
  triggers,
  freqs,
  spec_files,
  subjs,
  alpha,
  beta,
  delta,
  theta,
  stim_cog,
  stim_iaf,
  stim_paf,
  eeg_srate,
  fft_window,
  freq_bins,
  n_points,
  pwn_freqs,
  stim_freqvec,
  stim_psd,
  stim_psd_cor
)

# play sound to know that the script is finished
beep(2)
