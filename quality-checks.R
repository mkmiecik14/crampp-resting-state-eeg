# Early visualization script
# Matt Kmiecik
# 25 July 2022

# Purpose: making sure that data are looking good

source("r-prep.R") # preps R workspace

# Loads data - - - -
load("../output/psd-res.rda")

# collpses across block (subject-wise summary)
eyes_open_closed_ss <- 
  psd_res %>%
  group_by(ss, eyes, elec, freq, band) %>%
  summarise(m = mean(psd), n = n()) %>%
  ungroup()

# collapses across block (group-wise summary)
eyes_open_closed_sum <- 
  eyes_open_closed_ss %>%
  group_by(eyes, elec, freq, band) %>%
  summarise(
    M = mean(m),
    SD = sd(m),
    N = n(),
    SEM = SD/sqrt(N)
  ) %>%
  ungroup()

# Eyes open vs. closed plot
ggplot(eyes_open_closed_sum, aes(freq, M, group = eyes, color = eyes)) +
  geom_line() +
  geom_vline(xintercept = 12, color = rdgy_pal[8], linetype = 2) +
  coord_cartesian(xlim = c(0, 30), ylim = c(-40, -10)) +
  theme_bw() +
  scale_color_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  labs(x = "Frequency", y = "Power Spectral Density (dB)") +
  theme(legend.position = "bottom") +
  facet_wrap(~elec)
