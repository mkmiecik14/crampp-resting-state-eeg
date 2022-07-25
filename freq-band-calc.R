# Frequency Band Calculations
# Matt Kmiecik
# 25 July 2022

# Purpose: calculates frequency bands here so as not to recalculate them 
# for every analysis script

source("r-prep.R") # preps R workspace 

# Loads in preprocessed spectral results
load("../output/psd-res.rda")

# defining frequency bands
delta <- seq(.5, 4, .5)
theta <- seq(4, 8, .5)
alpha <- seq(8, 13, .5)
beta <- seq(13, 30, .5)


psd_res %>% 
  mutate(
    band = case_when(
      freq %in% delta ~ "delta",
      freq %in% theta ~ "theta",
      freq %in% alpha ~ "alpha",
      freq %in% beta ~ "beta",
      TRUE ~ "outside"
      )
    )

# EYES OPEN - - - -
# Summary collapsed across blocks
eyes_open <-
  psd_res %>% 
  filter(eyes == "open") %>%
  group_by(ss, elec, freq) %>%
  summarise(m = mean(psd), n = n()) %>%
  ungroup()

# inserting band descriptors after collapse
eyes_open_bands_ss <- 
  eyes_open %>%
  mutate(
    band = case_when(
      freq %in% delta ~ "delta",
      freq %in% theta ~ "theta",
      freq %in% alpha ~ "alpha",
      freq %in% beta ~ "beta",
      TRUE ~ "outside"
    )
    )

# calculates band-wise summary
eyes_open_bands <- 
  eyes_open_bands_ss %>%
  group_by(ss, elec, band) %>%
  summarise(M = mean(m), N = n()) %>%
  ungroup()
  