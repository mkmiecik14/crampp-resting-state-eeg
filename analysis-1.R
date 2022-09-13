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
  mutate(
    eyes = factor(eyes), # turns into a factor 
    block_mc = block - 4.5 # mean centers
  )

# proof that it is mean centered
iaf_res_long %>% 
  group_by(ss, name) %>% 
  summarise(m = mean(block_mc)) %>% 
  ungroup() %>% 
  filter(m != 0)

# summary: subject-wise
iaf_res_long_ss <- 
  iaf_res_long %>% 
  filter(complete.cases(value)) %>%
  group_by(ss, eyes, name) %>%
  summarise(m = mean(value), n = n()) %>%
  ungroup()

# summary: sample-wise
iaf_res_long_sum <- 
  iaf_res_long_ss %>%
  group_by(eyes, name) %>%
  summarise(
    M = mean(m), 
    SD = sd(m), 
    N = n(), 
    SEM = SD/sqrt(N), 
    LL = quantile(m, .025), 
    UL = quantile(m, .975)
    ) %>%
  ungroup()

# visualization
pn <- position_nudge(x = .3)
pj <- position_jitter(width = .1, height = 0)
ggplot(iaf_res_long_sum, aes(eyes, M, color = name, fill = name)) +
  geom_point(
    data = iaf_res_long_ss, 
    aes(y = m), 
    position = pj, 
    alpha = 1/3,
    shape = 16) +
  geom_bar(stat = "identity", width = .25, color = "black", position = pn) +
  geom_errorbar(
    aes(ymin=M-SEM, ymax = M+SEM), 
    width = .1, 
    color = "black", 
    position = pn
    ) +
  coord_cartesian(ylim = c(8, 13)) +
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
  labs(x = "Eyes", y = "Frequency (Hz)", caption = "SEM error bars.") +
  theme_bw() + 
  facet_wrap(~name) +
  theme(legend.position = "none")

# Linear mixed effects modeling - - - -

# Sets contrast (open > closed)
contrasts(iaf_res_long$eyes) <- cbind(open_vs_closed = c(-.5, .5))


# PEAK ALPHA FREQUENCY (via peak picking)

# Minimal model
mod_1_paf <- 
  lmer(
    value ~ 1 + block_mc*eyes + (1 | ss), 
    data = iaf_res_long %>% filter(name == "paf"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_1_paf)
check_model(mod_1_paf)

mod_2_paf <-
  lmer(
    value ~ 1 + block_mc*eyes + (1 + block_mc + eyes | ss), 
    data = iaf_res_long %>% filter(name == "paf"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_2_paf)
check_model(mod_2_paf)

# original optimizer was not converging; switched to bobyqa
# mod_2_paf_allfit <- allFit(mod_2_paf)
# mod_2__paf_allfit_OK <- mod_2_paf_allfit[sapply(mod_2_paf_allfit, is, "merMod")]
# lapply(mod_2__paf_allfit_OK, function(x) x@optinfo$conv$lme4$messages)

mod_3_paf <-
  lmer(
    value ~ 1 + block_mc*eyes + (1 + block_mc*eyes | ss), 
    data = iaf_res_long %>% filter(name == "paf"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_3_paf)
check_model(mod_3_paf)

# Plots interaction
interactions::interact_plot(
  mod_3_paf, 
  pred = block_mc, 
  modx = eyes, 
  plot.points = TRUE
)


# Model comparison
anova(mod_1_paf, mod_2_paf, mod_3_paf)

# CENTER OF GRAVITY (essentially a weighted average)

# Minimal model
mod_1_cog <- 
  lmer(
    value ~ 1 + block_mc*eyes + (1 | ss), 
    data = iaf_res_long %>% filter(name == "cog"),
    REML = TRUE
  )
summary(mod_1_cog)
check_model(mod_1_cog)

# block_mc does not have sufficient variability
mod_2_cog <-
  lmer(
    value ~ 1 + block_mc*eyes + (1 + eyes | ss), 
    data = iaf_res_long %>% filter(name == "cog"),
    REML = TRUE
  )
summary(mod_2_cog)
check_model(mod_2_cog)

# Model comparison
anova(mod_1_cog, mod_2_cog)


