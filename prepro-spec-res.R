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
stim_spectra  <- spec_res %>% map(1) # PSD
stim_freqs <- spec_res %>% map(2) # frequencies
freqs <- stim_freqs[[1]][,,1] # vector of frequencies

# trigger definitions
triggers <-
  tibble(
    trigger = c(111, 102, 103, 114, 105, 116, 117, 108),
    eyes = c(
      "open", "closed", "closed", "open", "closed", "open", "open", "closed"
      ),
    block = 1:8
  )

# Extracting and tidying spectral results ----

# defining frequency bands (for later)
delta <- seq(.5, 4, .5)
theta <- seq(4, 8, .5)
alpha <- seq(8, 13, .5)
beta <- seq(13, 30, .5)

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
  
# Saving out data ----
# Spectral results
save(psd_res, file = "../output/psd-res.rda")
write_csv(psd_res, file = "../output/psd-res.csv")

# channel locations
save(chan_locs, file = "../output/chan-locs.rda")
write_csv(chan_locs, file = "../output/chan-locs.csv")

# Cleans workspace ----
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
  theta
)
