# Early visualization script
# Matt Kmiecik
# 25 July 2022

# Purpose: making sure that data are looking good

source("r-prep.R") # preps R workspace
source("topo_tools.R") # for topo plotting tools

# parallel processing settings 
no_cores <- availableCores() - 1 # leave one core for background processing
plan(multicore, workers = no_cores) # plans the cores

###########################
#                         #
# BROADBAND SPECTRA in dB #
#                         #
###########################

# Loads data
load("../output/dB-res.rda")

dB_eyes_open_closed_ss <- 
  dB_res %>% 
  split(.$ss) %>%
  future_map_dfr(
    ~filter(.x, complete.cases(.)) %>%
      group_by(ss, eyes, elec, freq, band) %>%
      summarise(m = mean(psd), n = n()) %>%
      ungroup()
    )

# Removing objects not used anymore for speed
rm(dB_res); gc();

# collapses across block (group-wise summary)
dB_eyes_open_closed_sum <- 
  dB_eyes_open_closed_ss %>%
  filter(complete.cases(m)) %>%
  group_by(eyes, elec, freq, band) %>%
  summarise(
    M = mean(m),
    SD = sd(m),
    N = n(),
    SEM = SD/sqrt(N)
  ) %>%
  ungroup()

# Eyes open vs. closed plot
ggplot(dB_eyes_open_closed_sum, aes(freq, M, group = eyes, color = eyes, fill = eyes)) +
  geom_ribbon(
    aes(ymin = M-SEM, ymax = M+SEM, fill = eyes), color = NA, alpha = 1/2
  ) +
  geom_line() +
  geom_vline(xintercept = 12, color = rdgy_pal[8], linetype = 2) +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(-40, -10)) +
  theme_bw() +
  scale_color_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  scale_fill_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  labs(
    x = "Frequency", 
    y = "Power Spectral Density (dB)", 
    caption = "SEM error bars."
    ) +
  theme(legend.position = "bottom") +
  facet_wrap(~elec)

#############################
#                           #
# BROADBAND SPECTRA in PSD  #
#                           #
#############################

# Loads data
load("../output/psd-res.rda")

# collapses across block (subject-wise summary)
psd_eyes_open_closed_ss <- 
  psd_res %>%
  filter(complete.cases(psd)) %>%
  group_by(ss, eyes, elec, freq, band) %>%
  summarise(m = mean(psd), n = n()) %>%
  ungroup()

# Removing objects not used anymore for speed
rm(psd_res); gc();

# collapses across block (group-wise summary)
psd_eyes_open_closed_sum <- 
  psd_eyes_open_closed_ss %>%
  filter(complete.cases(.)) %>%
  group_by(eyes, elec, freq, band) %>%
  summarise(
    M = mean(m),
    SD = sd(m),
    N = n(),
    SEM = SD/sqrt(N)
  ) %>%
  ungroup()

# Eyes open vs. closed plot (FIX THIS PLOT! REMOVE COORD_CARTE)
ggplot(
  psd_eyes_open_closed_sum, 
  aes(freq, M, group = eyes, color = eyes, fill = eyes)
  ) +
  geom_ribbon(
    aes(ymin = M-SEM, ymax = M+SEM, fill = eyes), color = NA, alpha = 1/2
  ) +
  geom_line() +
  geom_vline(xintercept = 12, color = rdgy_pal[8], linetype = 2) +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, .5)) +
  theme_bw() +
  scale_color_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  scale_fill_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  labs(x = "Frequency", y = "Power Spectral Density [(uV/cm^2))^2/Hz]") +
  theme(legend.position = "bottom") +
  facet_wrap(~elec)

#############################################################
#                                                           #
# BROADBAND SPECTRA in PSD CORRECTED FOR PINK&WHITE NOISE   #
#                                                           #
#############################################################

# Loads data
load("../output/psd-cor-res.rda")

# collapses across block (subject-wise summary)
psd_cor_eyes_open_closed_ss <- 
  psd_cor_res %>%
  filter(complete.cases(psd_cor)) %>%
  group_by(ss, eyes, elec, freq, band) %>%
  summarise(m = mean(psd_cor), n = n()) %>%
  ungroup()

# collapses across block (group-wise summary)
psd_cor_eyes_open_closed_sum <- 
  psd_cor_eyes_open_closed_ss %>%
  group_by(eyes, elec, freq, band) %>%
  summarise(
    M = mean(m),
    SD = sd(m),
    N = n(),
    SEM = SD/sqrt(N)
  ) %>%
  ungroup()

# Eyes open vs. closed plot
ggplot(psd_cor_eyes_open_closed_sum, aes(freq, M, group = eyes, color = eyes)) +
  geom_ribbon(
    aes(ymin = M-SEM, ymax = M+SEM, fill = eyes), color = NA, alpha = 1/2
  ) +
  geom_line() +
  geom_vline(xintercept = 12, color = rdgy_pal[8], linetype = 2) +
  scale_x_continuous(limits = c(0, 25)) +
  scale_y_continuous(limits = c(0, .15)) + # may need adjusting
  theme_bw() +
  scale_color_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  scale_fill_manual(values = c(rdgy_pal[3], rdgy_pal[11])) +
  labs(
    x = "Frequency", 
    y = "Power Spectral Density [(uV/cm^2))^2/Hz]",
    title = "Corrected for Pink & White Noise"
    ) +
  theme(legend.position = "bottom") +
  facet_wrap(~elec)

#########################################
#                                       #
# COMPARING PSD vs. PSD NOISE CORRECTED #
#                                       #
#########################################

# combines both summary dfs from above
psd_both_eyes_open_closed_sum <-
  bind_rows(
    psd_eyes_open_closed_sum %>% mutate(noise = "uncorrected"),
    psd_cor_eyes_open_closed_sum %>% mutate(noise = "corrected")
  )

# Eyes open
ggplot(
  psd_both_eyes_open_closed_sum %>% filter(eyes == "open"), 
  aes(freq, M, group = noise, color = noise)
  ) +
  geom_line() +
  geom_vline(xintercept = 12, color = rdgy_pal[8], linetype = 2) +
  coord_cartesian(xlim = c(0, 15)) +
  theme_bw() +
  scale_color_manual(values = c(rdgy_pal[11], ghibli_palettes$PonyoMedium[3])) +
  labs(
    x = "Frequency", 
    y = "Power Spectral Density [(uV/cm^2))^2/Hz]",
    title = "Eyes Open"
  ) +
  theme(legend.position = "bottom") +
  facet_wrap(~elec)

# Eyes closed
ggplot(
  psd_both_eyes_open_closed_sum %>% filter(eyes == "closed"), 
  aes(freq, M, group = noise, color = noise)
) +
  geom_line() +
  geom_vline(xintercept = 12, color = rdgy_pal[8], linetype = 2) +
  coord_cartesian(xlim = c(0, 15)) +
  theme_bw() +
  scale_color_manual(values = c(rdgy_pal[11], ghibli_palettes$PonyoMedium[3])) +
  labs(
    x = "Frequency", 
    y = "Power Spectral Density [(uV/cm^2))^2/Hz]",
    title = "Eyes Closed"
  ) +
  theme(legend.position = "bottom") +
  facet_wrap(~elec)

########################
#                      #
# PEAK ALPHA FREQUENCY #
#                      #
########################

load("../output/paf-res.rda") # loads data

ggplot(paf_res, aes(paf, group = block, color = eyes)) +
  geom_density() +
  theme_bw() +
  facet_wrap(~elec)

# subject_wise summary
paf_ss <- 
  paf_res %>%
  filter(complete.cases(.)) %>%
  group_by(ss, elec, eyes) %>%
  summarise(m = mean(paf), n = n()) %>%
  ungroup()

# group-wise summary
paf_sum <- 
  paf_ss %>%
  group_by(elec, eyes) %>%
  summarise(M = mean(m), SD = sd(m), N = n(), SEM = SD/sqrt(N)) %>%
  ungroup()

pn <- position_nudge(x = .3)
pj <- position_jitter(width = .1, height = 0)
ggplot(paf_sum, aes(eyes, M, fill = eyes)) +
  geom_point(data = paf_ss, aes(y=m, color = eyes), alpha = 1/2, position = pj) +
  geom_bar(stat = "identity", width = .2, color = "black", position = pn) +
  geom_errorbar(aes(ymin=M-SEM, ymax = M+SEM), width = .1, position = pn) +
  scale_fill_manual(
    values = c(
      ghibli_palettes$MononokeLight[4], ghibli_palettes$MononokeLight[6]
      )
    ) +
  scale_color_manual(
    values = c(
      ghibli_palettes$MononokeLight[4], ghibli_palettes$MononokeLight[6]
    )
  ) +
  coord_cartesian(ylim = c(7.5, 13)) +
  labs(x = "Condition", y = "Peak Alpha Frequency", caption = "SEM error bars.") +
  theme_bw() +
  facet_wrap(~elec) +
  theme(legend.position = "none")

#####################
#                   #
# CENTER OF GRAVITY #
#                   #
#####################

load("../output/cog-res.rda") # loads data

ggplot(cog_res, aes(cog, group = block, color = eyes)) +
  geom_density() +
  theme_bw() +
  coord_cartesian(xlim = c(7.5, 13)) +
  facet_wrap(~elec)

# subject_wise summary
cog_ss <- 
  cog_res %>%
  filter(complete.cases(.)) %>%
  group_by(ss, elec, eyes) %>%
  summarise(m = mean(cog), n = n()) %>%
  ungroup()

# group-wise summary
cog_sum <- 
  cog_ss %>%
  group_by(elec, eyes) %>%
  summarise(M = mean(m), SD = sd(m), N = n(), SEM = SD/sqrt(N)) %>%
  ungroup()

pn <- position_nudge(x = .3)
pj <- position_jitter(width = .1, height = 0)
ggplot(cog_sum, aes(eyes, M, fill = eyes)) +
  geom_point(data = cog_ss, aes(y=m, color = eyes), alpha = 1/2, position = pj) +
  geom_bar(stat = "identity", width = .2, color = "black", position = pn) +
  geom_errorbar(aes(ymin=M-SEM, ymax = M+SEM), width = .1, position = pn) +
  scale_fill_manual(
    values = c(
      ghibli_palettes$MononokeLight[4], ghibli_palettes$MononokeLight[6]
    )
  ) +
  scale_color_manual(
    values = c(
      ghibli_palettes$MononokeLight[4], ghibli_palettes$MononokeLight[6]
    )
  ) +
  coord_cartesian(ylim = c(7.5, 13)) +
  labs(x = "Condition", y = "Center of Gravity", caption = "SEM error bars.") +
  theme_bw() +
  facet_wrap(~elec) +
  theme(legend.position = "none")

# Comparing PAF & COG

# combines PAF and COG data
paf_cog_sum <- 
  bind_rows(paf_sum %>% mutate(meas = "paf"), cog_sum %>% mutate(meas = "cog"))

pj <- position_jitter(width = .1, height = 0)
ggplot(paf_cog_sum, aes(eyes, M, fill = meas)) +
  geom_bar(stat = "identity", width = .4, color = "black", position = position_dodge(.5)) +
  geom_errorbar(aes(ymin=M-SEM, ymax = M+SEM), width = .1, position = position_dodge(.5)) +
  scale_fill_manual(
    values = c(
      ghibli_palettes$MononokeLight[3], ghibli_palettes$MononokeLight[5]
    )
  ) +
  scale_color_manual(
    values = c(
      ghibli_palettes$MononokeLight[3], ghibli_palettes$MononokeLight[5]
    )
  ) +
  coord_cartesian(ylim = c(7.5, 13)) +
  labs(x = "Condition", y = "Alpha Frequency", caption = "SEM error bars.") +
  theme_bw() +
  facet_wrap(~elec) +
  theme(legend.position = "bottom")


##############################
#                            #
# INDIVIDUAL ALPHA FREQUENCY #
#                            #
##############################

# IAF averages across electrodes that were able to generate a quality peak alpha
# frequency and center of gravity estimates
load("../output/iaf-res.rda")

# creates long format
iaf_res_long <- iaf_res %>% pivot_longer(c(-ss, -block, -eyes))

ggplot(iaf_res_long, aes(value, group = block, color = eyes)) +
  geom_density() +
  theme_bw() +
  coord_cartesian(xlim = c(7.5, 13)) +
  facet_wrap(~name)

# subject-wise summary
iaf_res_long_ss <- 
  iaf_res_long %>% 
  group_by(ss, eyes, name) %>% 
  summarise(m = mean(value), n = n()) %>%
  ungroup()

# sample-wise summary
iaf_res_long_sum <-
  iaf_res_long_ss %>%
  group_by(eyes, name) %>%
  summarise(M = mean(m), SD = sd(m), N = n(), SEM = SD/sqrt(N)) %>%
  ungroup()

# plot
pn <- position_nudge(x = .3)
pj <- position_jitter(width = .1, height = 0)
ggplot(iaf_res_long_sum, aes(eyes, M, fill = name)) +
  geom_bar(stat = "identity", width = .4, color = "black", position = position_dodge(.5)) +
  geom_errorbar(aes(ymin=M-SEM, ymax = M+SEM), width = .1, position = position_dodge(.5)) +
  scale_fill_manual(
    values = c(
      ghibli_palettes$MononokeLight[3], ghibli_palettes$MononokeLight[5]
    )
  ) +
  scale_color_manual(
    values = c(
      ghibli_palettes$MononokeLight[3], ghibli_palettes$MononokeLight[5]
    )
  ) +
  coord_cartesian(ylim = c(7.5, 13)) +
  labs(x = "Condition", y = "Individual Alpha Frequency", caption = "SEM error bars.") +
  theme_bw() +
  theme(legend.position = "bottom")

