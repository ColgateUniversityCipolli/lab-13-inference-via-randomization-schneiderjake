library(pwr)
library(tidyverse)
library(patchwork)
library(effectsize)
library(e1071)
library(boot)
library(boot.pval)

finch.dat <- read_csv("zebrafinches.csv")





#######################################################################
#USING BOOT PACKAGE
boot.tstat <- function(d, ind) {
  x    <- d[ind]
  n    <- length(x)
  (mean(x) - 0) / (sd(d) / sqrt(n))
}
################################# 
#Closer Data
closer0 <- finch.dat$closer - mean(finch.dat$closer)

boots <- boot(
  data      = closer0,
  statistic = boot.tstat,
  R         = 10000
)

boot.ci(boots, type = "bca")

boot.pval(
  boot_res   = boots,
  theta_null = 0
)

################################# 
#Further Data
further0 <- finch.dat$further - mean(finch.dat$further)
boots <- boot(
  data      = further0,
  statistic = boot.tstat,
  R         = 10000
)

boot.ci(boots, type = "bca")

boot.pval(
  boot_res   = boots,
  theta_null = 0
)

################################# 
#Difference Data
diff0 <- finch.dat$diff - mean(finch.dat$diff)
boots <- boot(
  data      = diff0,
  statistic = boot.tstat,
  R         = 10000
)

boot.ci(boots, type = "bca")

boot.pval(
  boot_res   = boots,
  theta_null = 0
)



















#######################################################################
#######################################################################
#task 2

################################################################################
# Bootstrap Test Confidence Interval Closer (T-statistic version)

R <- 100000
n <- nrow(finch.dat)
s <- sd(finch.dat$closer)
resamples <- tibble(tstats = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$closer,
                          size = n,
                          replace = TRUE)
  
  resamples$tstats[i] <- (mean(curr.resample) - 0) / (s / sqrt(n))
}

# No longer a confidence interval for x̄ — this is for the T-statistic distribution
# Shift so H0 is true (centered at 0)
(delta <- mean(resamples$tstats))

resamples <- resamples |>
  mutate(tstats.shifted = tstats - delta)

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







################################################################################
# Bootstrap Test Confidence Interval Further (T-statistic version)

R <- 100000
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



################################################################################
# Bootstrap Test Confidence Interval Diff (T-statistic version)

R <- 100000
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

# Use observed t-stat from original data
obs_t <- (mean(finch.dat$diff) - 0) / (s / sqrt(n))
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


