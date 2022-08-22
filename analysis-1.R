# Analysis 1 - Analyzing Individual Alpha Frequency (IAF) via global peak#
# picking and center of gravity estimates
# Matt Kmiecik
# Started 22 AUGUST 2022

source("r-prep.R")

# IAF averages across electrodes that were able to generate a quality peak alpha
# frequency and center of gravity estimates
load("../output/iaf-res.rda") # loads data

# creates long format
iaf_res_long <- 
  iaf_res %>% 
  pivot_longer(c(-ss, -block, -eyes)) %>%
  mutate(eyes = factor(eyes))

# Sets contrast (open > closed)
contrasts(iaf_res_long$eyes) <- cbind(open_vs_closed = c(-.5, .5))

# Linear mixed effects modeling

# PEAK ALPHA FREQUENCY (via peak picking)

# Minimal model
mod_1_paf <- 
  lmer(
    value ~ 1 + block + eyes + (1 | ss), 
    data = iaf_res_long %>% filter(name == "paf"),
    REML = TRUE
  )
summary(mod_1_paf)
check_model(mod_1_paf)

mod_2_paf <-
  lmer(
    value ~ 1 + block + eyes + (1 + block + eyes | ss), 
    data = iaf_res_long %>% filter(name == "paf"),
    REML = TRUE
  )
summary(mod_2_paf)
check_model(mod_2_paf)

# Model comparison
anova(mod_1_paf, mod_2_paf)

# CENTER OF GRAVITY (essentially a weighted average)

# Minimal model
mod_1_cog <- 
  lmer(
    value ~ 1 + block + eyes + (1 | ss), 
    data = iaf_res_long %>% filter(name == "cog"),
    REML = TRUE
  )
summary(mod_1_cog)
check_model(mod_1_cog)

mod_2_cog <-
  lmer(
    value ~ 1 + block + eyes + (1 + block + eyes | ss), 
    data = iaf_res_long %>% filter(name == "cog"),
    REML = TRUE
  )
summary(mod_2_cog)
check_model(mod_2_cog)

# Model comparison
anova(mod_1_cog, mod_2_cog)


# Visualization
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
pn <- position_nudge(x = .4)
pj <- position_jitter(width = .1, height = 0)
ggplot(iaf_res_long_sum, aes(eyes, M, fill = name)) +
  geom_bar(stat = "identity", width = .3, color = "black", position = pn) +
  geom_errorbar(aes(ymin=M-SEM, ymax = M+SEM), width = .1, position = pn) +
  geom_point(
    data = iaf_res_long_ss, 
    aes(y = m, color = name), 
    alpha = 1/2,
    position = pj,
    shape = 16
    ) +
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
  theme(legend.position = "none") +
  facet_wrap(~name)
