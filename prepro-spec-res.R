# Preprocessing script to organize the spectral data into R
# Matt Kmiecik
# Started 22 JULY 2022

source("r-prep.R") # preps R workspace 

# Loading various data ----

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
  mutate(ss_i = 1:n()) # to help join below / for specificity 

# Reading in and unpacking spectral results ----
spec_res <- 
  spec_files %>%
  map(~readMat(file.path("../output/", .x))) %>%
  map("spec.res") %>%
  map(~as.matrix(.x))

# Pulling out various spectral results ----
# resting-state blocks 1-8 along z dim in numerical order
stim_spectra  <- spec_res %>% map(1) # broadband PSD
stim_freqs <- spec_res %>% map(2) # frequencies
freqs <- stim_freqs[[1]][,,1] # vector of frequencies
stim_paf <- spec_res %>% map(3) # peak alpha frequency
stim_cog <- spec_res %>% map(4) # center of gravity
stim_iaf <- spec_res %>% map(5) # individual alpha freq


# trigger definitions
triggers <-
  tibble(
    trigger = c(111, 102, 103, 114, 105, 116, 117, 108),
    eyes = c(
      "open", "closed", "closed", "open", "closed", "open", "open", "closed"
      ),
    block = 1:8
  )

##############################
#                            #
# BROADBAND SPECTRAL RESULTS #
#                            #
##############################

# Extracting and tidying spectral results

# frequency resolution is determined in MATLAB spectral decomposition scripts
# so see those scripts for these numbers:
eeg_srate <- 256 # data sampled at 256 Hz
fft_window <- 8 # in seconds
n_points <- eeg_srate*fft_window # number of data points in FFT window
freq_bins <- eeg_srate/n_points # frequency bins

# defining frequency bands (for later)
delta <-  seq(.5, 4, freq_bins)
theta <-  seq(4, 7.5, freq_bins)
alpha <-  seq(7.5, 13, freq_bins)
beta <-   seq(13, 30, freq_bins)

# Each participant has:
# 30 rows (electrodes w/o ref) * 257 (frequencies) * 8 (rs blocks)
psd_res <-
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
  
# Saving out broadband spectral results
save(psd_res, file = "../output/psd-res.rda")
write_csv(psd_res, file = "../output/psd-res.csv")

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

# Saving out paf results
save(paf_res, file = "../output/paf-res.rda")
write_csv(paf_res, file = "../output/paf-res.csv")

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

# Saving out cog results
save(cog_res, file = "../output/cog-res.rda")
write_csv(cog_res, file = "../output/cog-res.csv")

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

# Saving out iaf results
save(iaf_res, file = "../output/iaf-res.rda")
write_csv(iaf_res, file = "../output/iaf-res.csv")


# Cleans work space ----
rm(
  chan_locs,
  psd_res,
  spec_res,
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
  cog_res,
  iaf_res,
  paf_res,
  stim_cog,
  stim_iaf,
  stim_paf,
  eeg_srate,
  fft_window,
  freq_bins,
  n_points
)
