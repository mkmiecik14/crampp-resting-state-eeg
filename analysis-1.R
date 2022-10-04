# Analysis 1 - Analyzing Individual Alpha Frequency (IAF) via global peak
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

#######
#     #
# MMH #
#     #
#######

# Looking at the relationship between peak alpha and MMH

# Loads in the principal components from the MMH paper
load("../data/mmh-res.RData")

# simplifies row-wise factor scores
fi <- 
  as_tibble(pca_res$Fixed.Data$ExPosition.Data$fi, rownames = "ss") %>%
  mutate(ss = as.numeric(ss))

# collapses across block
iaf_res_long_ss <- 
  iaf_res_long %>%
  group_by(ss, eyes, name) %>%
  summarise(m = mean(value), n = n()) %>%
  ungroup()

# wide format + pcs
iaf_res_long_ss_pcs <-
  iaf_res_long_ss %>% 
  pivot_wider(id_cols = ss, names_from = c(name, eyes), values_from = m) %>%
  left_join(., fi, by = "ss") %>%
  select(ss:V3) %>%
  rename(mmh = V1, ppt_sr = V2, bladder_hyper = V3)

iaf_res_long_ss_pcs %>% filter(complete.cases(.))

# Zero order correlations
# Computes bootstrapped correlations (seed was set earlier in script)
iaf_cor_res <-
  psych::corr.test(
    iaf_res_long_ss_pcs %>% select(-ss),
    use = "pairwise",
    method = "pearson", 
    adjust = "none",
    ci = TRUE,
    minlength = 100 # extends the abbreviations
  )
iaf_cor_res$ci # the results

# correlations into df
iaf_cor_res_ci <-
  as_tibble(iaf_cor_res$ci, rownames = "meas") %>%
  separate(meas, into = c("V1", "V2"), sep = "-")

# simplifies for plotting
iaf_cor_res_ci_simple <- 
  iaf_cor_res_ci %>% 
  filter(V2 %in% c("mmh", "ppt_sr", "bladder_hyper")) %>%
  filter(V1 %nin% c("mmh", "ppt_sr", "bladder_hyper")) %>%
  mutate(V2 = fct_relevel(V2, c("mmh", "ppt_sr", "bladder_hyper")))

# correlation plot
ggplot(iaf_cor_res_ci_simple, aes(V2, r)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_y_continuous(breaks = seq(-1, 1, .2), minor_breaks = NULL) +
  labs(x = "PCs", caption = "Error bars are bootstrapped 95% CIs.") +
  theme_bw() +
  facet_wrap(~V1, nrow = 1)

# Linear Mixed Modeling - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# block-wise data with principal components
iaf_res_long_pcs <- 
  iaf_res_long %>% 
  left_join(., fi, by = "ss") %>%
  rename(mmh = V1, ppt_sr = V2, bladder_hyper = V3)

# contrasts are already set
contrasts(iaf_res_long_pcs$eyes)

# PEAK ALPHA FREQUENCY (via peak picking) - - - - -

# Minimal model
mod_1_paf_pcs <- 
  lmer(
    value ~ 1 + mmh*block_mc*eyes + ppt_sr*block_mc*eyes + bladder_hyper*block_mc*eyes + (1 | ss), 
    data = iaf_res_long_pcs %>% filter(name == "paf"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_1_paf_pcs)
check_model(mod_1_paf_pcs)

# Plots interaction
interactions::interact_plot(
  mod_1_paf_pcs, 
  pred = mmh, 
  modx = eyes, 
  plot.points = TRUE
)

mod_2_paf_pcs <- 
  lmer(
    value ~ 1 + mmh*block_mc*eyes + ppt_sr*block_mc*eyes + bladder_hyper*block_mc*eyes + (1 + block_mc + eyes | ss), 
    data = iaf_res_long_pcs %>% filter(name == "paf"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_2_paf_pcs)
check_model(mod_2_paf_pcs)

mod_3_paf_pcs <- 
  lmer(
    value ~ 1 + mmh*block_mc*eyes + ppt_sr*block_mc*eyes + bladder_hyper*block_mc*eyes + (1 + block_mc*eyes | ss), 
    data = iaf_res_long_pcs %>% filter(name == "paf"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_3_paf_pcs)
check_model(mod_3_paf_pcs)

# Plots interaction
interactions::interact_plot(
  mod_1_paf_pcs, 
  pred = block_mc, 
  modx = ppt_sr, 
  plot.points = TRUE
)

# anova of models
anova(mod_1_paf_pcs, mod_2_paf_pcs, mod_3_paf_pcs)

# CENTER OF GRAVITY - - - - - 

# Minimal model
mod_1_cog_pcs <- 
  lmer(
    value ~ 1 + mmh*block_mc*eyes + ppt_sr*block_mc*eyes + bladder_hyper*block_mc*eyes + (1 | ss), 
    data = iaf_res_long_pcs %>% filter(name == "cog"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_1_cog_pcs)
check_model(mod_1_cog_pcs)

mod_2_cog_pcs <- 
  lmer(
    value ~ 1 + mmh*block_mc*eyes + ppt_sr*block_mc*eyes + bladder_hyper*block_mc*eyes + (1 + block_mc + eyes | ss), 
    data = iaf_res_long_pcs %>% filter(name == "cog"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_2_cog_pcs)
check_model(mod_2_cog_pcs)

mod_3_cog_pcs <- 
  lmer(
    value ~ 1 + mmh*block_mc*eyes + ppt_sr*block_mc*eyes + bladder_hyper*block_mc*eyes + (1 + block_mc*eyes | ss), 
    data = iaf_res_long_pcs %>% filter(name == "cog"),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa")
  )
summary(mod_3_cog_pcs)
check_model(mod_3_cog_pcs)

# Plots interaction
interactions::interact_plot(
  mod_1_cog_pcs, 
  pred = block_mc, 
  modx = ppt_sr, 
  plot.points = TRUE
)

# anova of models
anova(mod_1_cog_pcs, mod_2_cog_pcs, mod_3_cog_pcs)
