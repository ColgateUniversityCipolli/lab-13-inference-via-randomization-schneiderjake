library(pwr)
library(tidyverse)
library(patchwork)
library(effectsize)
library(e1071)
library(boot)
library(boot.pval)

finch.dat <- read_csv("zebrafinches.csv")


#######################################################################
#######################################################################
#task 2

boot.mean <- function(d, i){
  mean(d[i])
}


################################################################################
# Bootstrap Test Confidence Interval Closer (T-statistic version)

R <- 5000
n <- nrow(finch.dat)
s <- sd(finch.dat$closer)
resamples <- tibble(tstats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$closer,
                          size = n,
                          replace = TRUE)
  
  resamples$tstats[i] <- (mean(curr.resample) - 0) / (s / sqrt(n))
}

# No longer a confidence interval for  this is for the T-statistic distribution
# Shift so H0 is true (centered at 0)
(delta <- mean(resamples$tstats))

resamples <- resamples |>
  mutate(tstats.shifted = tstats - delta)

resamples.null.closer <- resamples$tstats.shifted

mean(resamples.null.closer) #this is close to zero

# Use observed t-stat from original data
obs_t <- (mean(finch.dat$closer) - 0) / (s / sqrt(n))
low <- -abs(obs_t)
high <- abs(obs_t)

resamples |>
  summarize(mean = mean(tstats.shifted),
            p.low = mean(tstats.shifted <= low),
            p.high = mean(tstats.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()

# Plot the Bootstrap Test
ggdat.obs <- tibble(xbar = obs_t,
                    y = 0) |>
  mutate(mirror = -xbar)  # reflection around 0 for two-sided

crit_upper <- quantile(resamples$tstats.shifted, 0.95)
crit_upper

bootstrap.plot.closer.tstat <- ggplot() +
  # Plot Resampling Distribution
  stat_density(data = resamples, aes(x = tstats), geom = "line", color = "lightgrey") +
  stat_density(data = resamples, aes(x = tstats.shifted), geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  # Plot Observation
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"), shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name = bquote(t),
    limits = range(c(resamples$tstats, resamples$tstats.shifted)) * 1.1,
    breaks = pretty(c(resamples$tstats, resamples$tstats.shifted), n = 5)
  ) +
  ylab("Density") +
  ggtitle("Bootstrap Hypothesis Test (T-Statistic)")


boots <- boot(data = finch.dat$closer,
              statistic = boot.mean,
              R = R)
ci.closer  <- boot.ci(boots, type = "bca")$bca[4:5] 
ci.closer <- sprintf("[%.3f, %.3f]", ci.closer[1], ci.closer[2])

################################################################################
# Bootstrap Test Confidence Interval Further (T-statistic version)

R <- 5000
n <- nrow(finch.dat)
s <- sd(finch.dat$further)
resamples <- tibble(tstats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$further,
                          size = n,
                          replace = TRUE)
  
  resamples$tstats[i] <- (mean(curr.resample) - 0) / (s / sqrt(n))
}

# Shift so H0 is true (centered at 0)
(delta <- mean(resamples$tstats))

resamples <- resamples |>
  mutate(tstats.shifted = tstats - delta)

resamples.null.further <- resamples$tstats.shifted

mean(resamples.null.further) #this is close to zero

# Use observed t-stat from original data
obs_t <- (mean(finch.dat$further) - 0) / (s / sqrt(n))
low <- -abs(obs_t)
high <- abs(obs_t)

resamples |>
  summarize(mean = mean(tstats.shifted),
            p.low = mean(tstats.shifted <= low),
            p.high = mean(tstats.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()

# Plot
ggdat.obs <- tibble(xbar = obs_t,
                    y = 0) |>
  mutate(mirror = -xbar)

crit_lower <- quantile(resamples$tstats.shifted, 0.05)
crit_lower

bootstrap.plot.further.tstat <- ggplot() +
  stat_density(data = resamples, aes(x = tstats), geom = "line", color = "lightgrey") +
  stat_density(data = resamples, aes(x = tstats.shifted), geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  geom_point(data = ggdat.obs, aes(x = xbar, y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"), shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name = bquote(t),
    limits = range(c(resamples$tstats, resamples$tstats.shifted)) * 1.1,
    breaks = pretty(c(resamples$tstats, resamples$tstats.shifted), n = 5)
  ) +
  ylab("Density") +
  ggtitle("Bootstrap Hypothesis Test (T-Statistic, Further)")

boots <- boot(data = finch.dat$further,
              statistic = boot.mean,
              R = R)
ci.further <- boot.ci(boots, type = "bca")$bca[4:5]
ci.further <- sprintf("[%.3f, %.3f]", ci.further[1], ci.further[2])

################################################################################
# Bootstrap Test Confidence Interval Diff (T-statistic version)

R <- 5000
n <- nrow(finch.dat)
s <- sd(finch.dat$diff)
resamples <- tibble(tstats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$diff,
                          size = n,
                          replace = TRUE)
  
  resamples$tstats[i] <- (mean(curr.resample) - 0) / (s / sqrt(n))
}

# Shift so H0 is true (centered at 0)
(delta <- mean(resamples$tstats))

resamples <- resamples |>
  mutate(tstats.shifted = tstats - delta)

resamples.null.diff <- resamples$tstats.shifted

mean(resamples.null.diff) #this is close to zero

# Use observed t-stat from original data
obs_t <- (mean(finch.dat$diff) - 0) / (s / sqrt(n))
low <- -abs(obs_t)
high <- abs(obs_t)

resamples |>
  summarize(mean = mean(tstats.shifted),
            p.low = mean(tstats.shifted <= low),
            p.high = mean(tstats.shifted >= high))|>
  mutate(p = p.low + p.high) 


# Plot
ggdat.obs <- tibble(xbar = obs_t,
                    y = 0) |>
  mutate(mirror = -xbar)

## 2.5 % and 97.5 % quantiles of the null (shifted) bootstrap distribution
crit_two_sided <- quantile(resamples$tstats.shifted, c(0.025, 0.975))
crit_two_sided

bootstrap.plot.diff.tstat <- ggplot() +
  stat_density(data = resamples, aes(x = tstats), geom = "line", color = "lightgrey") +
  stat_density(data = resamples, aes(x = tstats.shifted), geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  geom_point(data = ggdat.obs, aes(x = xbar, y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"), shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name = bquote(t),
    limits = range(c(resamples$tstats, resamples$tstats.shifted)) * 1.1,
    breaks = pretty(c(resamples$tstats, resamples$tstats.shifted), n = 5)
  ) +
  ylab("Density") +
  ggtitle("Bootstrap Hypothesis Test (T-Statistic, Diff)")



boots <- boot(data = finch.dat$diff,
              statistic = boot.mean,
              R = R)
ci.diff    <- boot.ci(boots, type = "bca")$bca[4:5]   
ci.diff <- sprintf("[%.3f, %.3f]", ci.diff[1], ci.diff[2])

ci.table <- tibble( 
  Variable = c("Closer", "Further", "Diff"),
  `Bootstrap 95% CI(BCa)` = c(ci.closer, ci.further, ci.diff)
)



p_ttest_closer  <- t.test(finch.dat$closer, mu = 0)$p.value
p_ttest_further <- t.test(finch.dat$further, mu = 0)$p.value
p_ttest_diff    <- t.test(finch.dat$diff, mu = 0)$p.value

# -----------------------------
# Bootstrap p-values (two-sided)
# -----------------------------
p_boot_closer  <- mean(resamples.null.closer <= low) + mean(resamples.null.closer >= high)
p_boot_further <- mean(resamples.null.further <= low) + mean(resamples.null.further >= high)
p_boot_diff    <- mean(resamples.null.diff <= low) + mean(resamples.null.diff >= high)

# -----------------------------
# Combine into a comparison table
# -----------------------------
pval.comparison <- tibble(
  Variable = c("Closer", "Further", "Diff"),
  `Bootstrap p-value` = round(c(p_boot_closer, p_boot_further, p_boot_diff), 4),
  `T-Test p-value` = round(c(p_ttest_closer, p_ttest_further, p_ttest_diff), 4)
)

##########################
#2c 
df        <- nrow(finch.dat) - 1
alpha     <- 0.05

# Helper: obs t and null-shifted bootstrap t’s
boot_tstats <- function(x, R = 5000, mu0 = 0) {
  n  <- length(x)
  s  <- sd(x)
  t_stats <- replicate(R, {
    xs <- sample(x, n, TRUE)
    (mean(xs) - mu0)/(s/sqrt(n))
  })
  t_shift <- t_stats - mean(t_stats)
  t_obs   <- (mean(x) - mu0)/(s/sqrt(n))
  list(
    t_obs = t_obs,
    qs    = quantile(t_shift, c(.025, .05, .95, .975), names = FALSE)
  )
}

# Compute for each variable
closer  <- boot_tstats(finch.dat$closer)
further <- boot_tstats(finch.dat$further)
diff    <- boot_tstats(finch.dat$diff)

# Theoretical cut-offs
t_left05    <- qt(alpha,    df)    # 5% left-tail
t_right95   <- qt(1-alpha,  df)    # 95% right-tail
t_two_low   <- qt(alpha/2,  df)    # 2.5% two-tail
t_two_high  <- qt(1-alpha/2,df)    # 97.5% two-tail

compare <- tibble(
  Variable           = c("Closer", "Further", "Diff"),
  `Oberserved t`              = c(closer$t_obs,    further$t_obs,    diff$t_obs),
  `Boot 5th percentile`        = c(closer$qs[2],    further$qs[2],    diff$qs[2]),
  `Boot 95th percentile`      = c(closer$qs[3],    further$qs[3],    diff$qs[3]),
  `Boot 2.5th percentile` = c(closer$qs[1],    further$qs[1],    diff$qs[1]),
  `Boot 97.5th percentile` = c(closer$qs[4],    further$qs[4],    diff$qs[4]),
  `Ttest 5th percentile`       = t_left05,
  `Ttest 95th percentile`     = t_right95,
  `Ttest 2.5th percentile`= t_two_low,
  `Ttest 97.5th percentile`= t_two_high
)

compare |> 
  mutate(across(-Variable, ~round(.x, 3)))
# Left‐tailed test
# Left‐tailed test
left_tail <- compare|>
  select(
    Variable,
    `Oberserved t`,
    `Boot 5th percentile`,
    `Ttest 5th percentile`
  ) 

# Right‐tailed test
right_tail <- compare|>
  select(
    Variable,
    `Oberserved t`,
    `Boot 95th percentile`,
    `Ttest 95th percentile`
  )

# Two‐tailed test
two_tail <- compare |>
  select(
    Variable,
    `Oberserved t`,
    `Boot 2.5th percentile`,
    `Boot 97.5th percentile`,
    `Ttest 2.5th percentile`,
    `Ttest 97.5th percentile`
  ) 
  


# View each
view(left_tail)
view(right_tail)
view(two_tail)

